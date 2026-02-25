clc; clear; close all;

% 路径设置
projectRoot = fileparts(mfilename('fullpath')); 
cd(projectRoot);
addpath(genpath(projectRoot));

%% ==========================================
%% 0. 设定要测试的自转速度数组 (参数扫掠 0 - 0.1 r/s)
%% ==========================================
% 定义一系列候选自转速度 (单位: r/s，转每秒)
spin_test_rps = linspace(0, 0.1, 6); % 测试 0, 0.02, 0.04, 0.06, 0.08, 0.1 r/s
% 转换为底层物理引擎需要的角速度 (rad/s)
spin_test_values = spin_test_rps * 2 * pi; 

num_tests = length(spin_test_values);
deicing_rates = zeros(1, num_tests);

fprintf('开始进给速度(自转)寻优，测试范围: 0 ~ 0.1 r/s，共需测试 %d 组...\n', num_tests);

for test_idx = 1:num_tests
    current_spin_rad = spin_test_values(test_idx);
    current_spin_rps = spin_test_rps(test_idx);
    
    fprintf('\n==================================================\n');
    fprintf('>>> 正在测试第 %d/%d 组 | 设定进给速度: %.3f r/s (%.3f rad/s) <<<\n', ...
            test_idx, num_tests, current_spin_rps, current_spin_rad);
    fprintf('==================================================\n');

    %% 1. 加载配置 (每次循环重新加载，确保系统底层状态干净无残留)
    robotDescriptionDeicing; 
    
    % 【核心修改】：通过代码覆盖原来的自转速度
    env.contact_params.omega_spin = current_spin_rad; 

    % 强制参数修正
    if ~isfield(sim_params, 'omega'), sim_params.omega = [0; 0; sim_params.omega_target]; end
    env.ext_force_list = unique([env.ext_force_list, "centrifugal", "coriolis"]);

    % 创建几何与对象
    nodes_original = nodes;
    [nodes, edges, rod_edges, shell_edges, rod_shell_joint_edges, rod_shell_joint_total_edges, face_nodes, face_edges, face_shell_edges, ...
        elStretchRod, elStretchShell, elBendRod, elBendSign, elBendShell, sign_faces, face_unit_norms] ...
        = createGeometry(nodes, edges, face_nodes);

    twist_angles = zeros(size(rod_edges,1)+size(rod_shell_joint_total_edges,1),1);
    [environment, imc] = createEnvironmentAndIMCStructs(env, geom, material, sim_params);

    softRobot = MultiRod(geom, material, twist_angles, nodes, edges, rod_edges, shell_edges, ...
        rod_shell_joint_edges, rod_shell_joint_total_edges, face_nodes, sign_faces, ...
        face_edges, face_shell_edges, sim_params, environment);
    softRobot.nodes_local = nodes_original; 

    % 创建弹簧
    n_stretch = size(elStretchRod,1) + size(elStretchShell,1);
    n_bend_twist = size(elBendRod,1);
    if n_stretch==0, stretch_springs=[]; else
        for s=1:n_stretch
            if s <= size(elStretchRod,1), stretch_springs(s) = stretchSpring(softRobot.refLen(s), elStretchRod(s,:), softRobot);
            else, stretch_springs(s) = stretchSpring(softRobot.refLen(s), elStretchShell(s-size(elStretchRod,1),:), softRobot, softRobot.ks(s)); end
        end
    end
    if n_bend_twist==0, bend_twist_springs=[]; else
        for b=1:n_bend_twist
            bend_twist_springs(b) = bendTwistSpring(elBendRod(b,:), elBendSign(b,:), [0 0], 0, softRobot);
        end
    end
    hinge_springs = []; triangle_springs = [];

    % 初始化系统
    softRobot = computeSpaceParallel(softRobot);
    theta = softRobot.q0(3*softRobot.n_nodes+1 : 3*softRobot.n_nodes+softRobot.n_edges_dof);
    [softRobot.m1, softRobot.m2] = computeMaterialDirectors(softRobot.a1, softRobot.a2, theta);
    bend_twist_springs = setkappa(softRobot, bend_twist_springs);
    softRobot.undef_refTwist = computeRefTwist_bend_twist_spring(bend_twist_springs, softRobot.a1, softRobot.tangent, zeros(n_bend_twist,1));
    softRobot.refTwist = computeRefTwist_bend_twist_spring(bend_twist_springs, softRobot.a1, softRobot.tangent, softRobot.undef_refTwist);

    softRobot.fixed_nodes = fixed_node_indices; 
    softRobot.fixed_edges = [];
    [softRobot.fixedDOF, softRobot.freeDOF] = FindFixedFreeDOF(softRobot.fixed_nodes, softRobot.fixed_edges, softRobot.n_DOF, softRobot.n_nodes);

    %% 2. 物理仿真循环
    Nsteps = round(sim_params.totalTime/sim_params.dt);
    ctime = 0; 
    break_times = inf(1, imc.num_ice);

    for timeStep = 1:Nsteps
        if ctime < sim_params.ramp_time
            current_omega_mag = (ctime / sim_params.ramp_time) * sim_params.omega_target;
        else
            current_omega_mag = sim_params.omega_target;
        end
        
        sim_params.omega = [0; 0; current_omega_mag]; 
        imc.omega_mag = current_omega_mag;            
        imc.theta_accumulated = imc.theta_accumulated + current_omega_mag * sim_params.dt;

        if(sim_params.use_midedge), tau_0 = updatePreComp_without_sign(softRobot.q, softRobot); else, tau_0 = []; end
        
        [softRobot, stretch_springs, bend_twist_springs, hinge_springs, force_now, imc] = ...
            timeStepper(softRobot, stretch_springs, bend_twist_springs, hinge_springs, ...
            triangle_springs, tau_0, environment, imc, sim_params, ctime);
            
        % 捕获断裂事件
        if any(imc.is_broken)
            newly_broken = imc.is_broken & isinf(break_times);
            if any(newly_broken)
                broken_indices = find(newly_broken);
                for idx = broken_indices
                    break_times(idx) = ctime;
                end
            end
        end
        
        ctime = ctime + sim_params.dt;
        softRobot.q0 = softRobot.q; 
    end
    
    % 记录这一组测试的除冰率
    broken_count = sum(imc.is_broken);
    deicing_rates(test_idx) = (broken_count / imc.num_ice) * 100;
    fprintf('>>> 本组测试完成 | 进给速度: %.3f r/s | 击碎: %d/%d | 除冰率: %.1f%%\n', ...
            current_spin_rps, broken_count, imc.num_ice, deicing_rates(test_idx));
end

%% 3. 绘制寻优结果对比曲线
figure('Name', '除冰机器人进给速度寻优', 'Color', 'w');
plot(spin_test_rps, deicing_rates, '-ko', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
title('绕绝缘子进给速度与除冰率的关系', 'FontSize', 14);
xlabel('进给速度 (r/s)', 'FontSize', 12);
ylabel('除冰率 (%)', 'FontSize', 12);
ylim([0 110]);
xlim([0 max(spin_test_rps)]);
grid on;

% 找到除冰率最高且速度最大的“最优作业点”（兼顾除冰率和效率）
best_rate = max(deicing_rates);
% 寻找达到最高除冰率的所有索引，取其中速度最大的一个（更高效）
best_indices = find(deicing_rates == best_rate);
best_idx = best_indices(end); 
best_spin_rps = spin_test_rps(best_idx);

hold on;
plot(best_spin_rps, best_rate, 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
text(best_spin_rps, best_rate - 8, sprintf('推荐作业速度: %.2f r/s\n(除冰率: %.1f%%)', best_spin_rps, best_rate), ...
    'HorizontalAlignment', 'center', 'Color', 'r', 'FontWeight', 'bold');

fprintf('\n====== 寻优结束 ======\n');
fprintf('测试完毕！在达到最高除冰率 %.1f%% 的前提下，最高效的进给速度为: %.3f r/s\n', best_rate, best_spin_rps);