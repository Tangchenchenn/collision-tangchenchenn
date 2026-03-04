% length_vs_r.m
% 目的：扫描不同轮毂半径，保持 L_rope : R_hub = 5:1
%       观察绳子在1000 rpm旋转下的稳定扫掠半径

clc; clearvars; close all; % 使用 clearvars 避免性能下降

projectRoot = fileparts(mfilename('fullpath')); 
cd(projectRoot);
addpath(genpath(projectRoot));

%% 1. 参数设置
RPM_target     = 1000;
omega_target   = RPM_target * 2 * pi / 60;     % rad/s

ratio_L_R      = 5;
% 【修复笔误】将 0.5 改为了 0.050
R_hub_list     = [0.010, 0.060];
n_cases        = length(R_hub_list);

stable_radius  = nan(n_cases,1);
tip_radius_std = nan(n_cases,1);               

%% 2. 批量模拟
for ic = 1:n_cases

    R_hub = R_hub_list(ic);
    fprintf('\n=== 案例 %2d / %2d   R_hub = %.4f m ===\n', ic, n_cases, R_hub);

    % 设置轮毂半径
    geom.hub_radius = R_hub;
    
    % 先加载底层配置文件 (这会载入默认的材料、几何结构)
    robot2DescriptionDeicing;           

    % 【核心修复：在配置文件加载后，再强制覆盖我们这个特定实验需要的参数】
    sim_params.totalTime   = 1.5;   % 只需2秒即可达到稳定
    sim_params.dt          = 1e-4;  % 适当放大步长加快扫描速度
    sim_params.tol         = 5e-4;
    sim_params.maximum_iter= 25;
    sim_params.log_data    = true;
    sim_params.logStep     = 20;
    sim_params.ramp_time   = 0.3;
    sim_params.omega_target= omega_target;

    env.ext_force_list     = ["gravity", "viscous", "aerodynamic", "centrifugal", "coriolis"];
    env.eta                = 0.05;

    % 确保严格关闭接触引擎
    env.selfContact        = false;
    env.floorContact       = false;

    % 保存初始节点
    nodes_original = nodes;
    
    % ...... 下面的 [nodes, edges, rod_edges...] 等创建几何对象的代码保持不变 ......

    % 创建几何对象、弹簧等
    [nodes, edges, rod_edges, shell_edges, rod_shell_joint_edges, rod_shell_joint_total_edges, ...
     face_nodes, face_edges, face_shell_edges, elStretchRod, elStretchShell, elBendRod, ...
     elBendSign, elBendShell, sign_faces, face_unit_norms] ...
        = createGeometry(nodes, edges, face_nodes);

    twist_angles = zeros(size(rod_edges,1) + size(rod_shell_joint_total_edges,1), 1);

    [environment, ~] = createEnvironmentAndIMCStructs(env, geom, material, sim_params);

    softRobot = MultiRod(geom, material, twist_angles, nodes, edges, rod_edges, shell_edges, ...
        rod_shell_joint_edges, rod_shell_joint_total_edges, face_nodes, sign_faces, ...
        face_edges, face_shell_edges, sim_params, environment);
    softRobot.nodes_local = nodes_original;

    % ─────────────── 创建弹簧（使用对象数组，与 timeStepper/getFbJb 等兼容） ───────────────
    n_stretch    = size(elStretchRod,1) + size(elStretchShell,1);
    n_bend_twist = size(elBendRod,1);

    if n_stretch > 0
        for s = 1:n_stretch
            if s <= size(elStretchRod,1)
                stretch_springs(s) = stretchSpring(softRobot.refLen(s), elStretchRod(s,:), softRobot);
            else
                stretch_springs(s) = stretchSpring(softRobot.refLen(s), ...
                    elStretchShell(s-size(elStretchRod,1),:), softRobot, softRobot.ks(s));
            end
        end
    else
        stretch_springs = [];
    end

    if n_bend_twist > 0
        for b = 1:n_bend_twist
            bend_twist_springs(b) = bendTwistSpring(elBendRod(b,:), elBendSign(b,:), [0 0], 0, softRobot);
        end
    else
        bend_twist_springs = [];
    end

    hinge_springs    = [];
    triangle_springs = [];

    % ─────────────── 初始化参考框架 ───────────────
    softRobot = computeSpaceParallel(softRobot);
    theta = softRobot.q0(3*softRobot.n_nodes+1 : 3*softRobot.n_nodes+softRobot.n_edges_dof);
    [softRobot.m1, softRobot.m2] = computeMaterialDirectors(softRobot.a1, softRobot.a2, theta);

    % ─────────────── 调用 setkappa 并计算参考扭转 ───────────────
    if ~isempty(bend_twist_springs)
        bend_twist_springs = setkappa(softRobot, bend_twist_springs);

        softRobot.undef_refTwist = computeRefTwist_bend_twist_spring(...
            bend_twist_springs, softRobot.a1, softRobot.tangent, ...
            zeros(numel(bend_twist_springs),1));

        softRobot.refTwist = computeRefTwist_bend_twist_spring(...
            bend_twist_springs, softRobot.a1, softRobot.tangent, ...
            softRobot.undef_refTwist);
    else
        softRobot.undef_refTwist = [];
        softRobot.refTwist = [];
    end

    % 边界条件
    softRobot.fixed_nodes = fixed_node_indices;
    softRobot.fixed_edges = [];
    [softRobot.fixedDOF, softRobot.freeDOF] = FindFixedFreeDOF(softRobot.fixed_nodes, softRobot.fixed_edges, softRobot.n_DOF, softRobot.n_nodes);

    % ─────────────── 仿真循环 ───────────────
    Nsteps = round(sim_params.totalTime / sim_params.dt);
    ctime = 0;
    tip_r_log = zeros(1, floor(Nsteps/sim_params.logStep)+1);
    log_idx = 0;

    for step = 1:Nsteps

        % 转速爬坡
        if ctime < sim_params.ramp_time
            ratio = ctime / sim_params.ramp_time;
            current_omega = ratio * sim_params.omega_target;
        else
            current_omega = sim_params.omega_target;
        end
        sim_params.omega = [0; 0; current_omega];

        % 时间步进（不传 imc）
        [softRobot, stretch_springs, bend_twist_springs, hinge_springs, ~, ~] = ...
            timeStepper(softRobot, stretch_springs, bend_twist_springs, hinge_springs, ...
                        triangle_springs, [], environment, [], sim_params, ctime);

        ctime = ctime + sim_params.dt;
        softRobot.q0 = softRobot.q;

        % 记录最远端距离
        if mod(step, sim_params.logStep) == 0
            log_idx = log_idx + 1;
            q_pos = reshape(softRobot.q(1:3*softRobot.n_nodes), 3, []);
            r_all = sqrt(sum(q_pos.^2, 1));
            tip_r_log(log_idx) = max(r_all);
        end
    end

    % 后处理：取最后 40% 的数据计算平均与标准差
    n_stable = round(length(tip_r_log)*0.4);
    if n_stable < 1, n_stable = 1; end
    stable_tip_r = tip_r_log(end-n_stable+1 : end);

    stable_radius(ic)  = mean(stable_tip_r);
    tip_radius_std(ic) = std(stable_tip_r);

    fprintf('  稳定扫掠半径 ≈ %.4f ± %.4f m\n', stable_radius(ic), tip_radius_std(ic));
end


%% 4. 结果绘图
figure('Color','w','Position',[300 200 1100 650]);

% ==================== 计算目标点及其对应的各项数值 ====================
target_radius = 0.524 / 2;  % 目标扫掠半径 (Y轴) = 0.262

% 1. 反向插值求出对应的轮毂半径 (X轴)
target_R_hub = interp1(stable_radius, R_hub_list, target_radius, 'linear'); 

% 2. 计算该轮毂半径对应的“理论最大值” (使用原代码的理论公式)
target_theoretical = target_R_hub * ratio_L_R; 

% 3. 插值求出该轮毂半径对应的“波动范围(标准差)”
target_std = interp1(R_hub_list, tip_radius_std, target_R_hub, 'linear');
% ======================================================================

% -------- 左图：扫掠半径 --------
subplot(1,2,1);
plot(R_hub_list, stable_radius, 'b.-', 'LineWidth',2.2, 'MarkerSize',24);
hold on;
plot(R_hub_list, R_hub_list * ratio_L_R, 'k--', 'LineWidth',1.4);

% 标注1：模拟出来的稳定扫掠半径点 (红星)
plot(target_R_hub, target_radius, 'rp', 'MarkerSize', 14, 'MarkerFaceColor', 'r');
text(target_R_hub, target_radius, sprintf('  模拟: (%.4f, %.3f)', target_R_hub, target_radius), ...
    'Color', 'r', 'FontSize', 11, 'FontWeight', 'bold', 'VerticalAlignment', 'bottom');

% 标注2：理论最大值点 (品红星)
plot(target_R_hub, target_theoretical, 'mp', 'MarkerSize', 14, 'MarkerFaceColor', 'm');
text(target_R_hub, target_theoretical, sprintf('  理论: (%.4f, %.3f)', target_R_hub, target_theoretical), ...
    'Color', 'm', 'FontSize', 11, 'FontWeight', 'bold', 'VerticalAlignment', 'top');

xlabel('轮毂半径 R_{hub} (m)');
ylabel('稳定扫掠半径 (m)');
title(sprintf('转速 %d rpm，L:R_{hub}=%d:1', RPM_target, ratio_L_R));
grid on; box on; 
legend('模拟结果','理论最大（直线）','目标模拟点', '对应理论点', 'Location','southeast');

% -------- 右图：波动幅度 --------
subplot(1,2,2);
plot(R_hub_list, tip_radius_std, 'r.-', 'LineWidth',2, 'MarkerSize',20);
hold on;

% 标注3：对应的波动范围标准差点 (蓝星)
plot(target_R_hub, target_std, 'bp', 'MarkerSize', 14, 'MarkerFaceColor', 'b');
text(target_R_hub, target_std, sprintf('  波动: (%.4f, %.5f)', target_R_hub, target_std), ...
    'Color', 'b', 'FontSize', 11, 'FontWeight', 'bold', 'VerticalAlignment', 'bottom');

xlabel('轮毂半径 R_{hub} (m)');
ylabel('扫掠半径波动标准差 (m)');
title('稳定阶段波动幅度');
grid on; box on;
legend('波动标准差', '目标对应波动值', 'Location','northeast');

sgtitle('绳子高速旋转扫掠半径随轮毂尺寸变化');

% 保存图片
saveas(gcf, sprintf('rope_%drpm_ratio%d.png', RPM_target, ratio_L_R));

% 在命令行输出具体数据，方便核对
disp(table(R_hub_list(:), stable_radius(:), tip_radius_std(:), ...
    'VariableNames', {'R_hub_m', '稳定半径_m', '波动_std_m'}));

fprintf('\n================== 目标点分析结果 ==================\n');
fprintf('目标扫掠半径:       %.3f m\n', target_radius);
fprintf('-> 插值求得轮毂半径: %.4f m\n', target_R_hub);
fprintf('-> 对应理论最大半径: %.3f m\n', target_theoretical);
fprintf('-> 对应扫掠波动幅度: %.5f m (标准差)\n', target_std);
fprintf('====================================================\n');