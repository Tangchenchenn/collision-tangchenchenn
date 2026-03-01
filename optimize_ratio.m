clc; clear all; close all;

% 路径设置
projectRoot = fileparts(mfilename('fullpath')); 
cd(projectRoot);
addpath(genpath(projectRoot));

%% 1. 寻优参数设置
total_length = 0.262; % 单侧绳子长度 + 轮毂半径 = 0.262m
ratios = 9:-1:1;      % 刚柔比 (绳长 : 轮毂半径) 9:1 到 1:1
num_cases = length(ratios);

% 网格控制参数：固定每段的物理长度，而不是固定总节点数
target_l_bar = 0.015; % 目标网格单元长度约为 1.5 cm

% 结果记录数组
res_peak_torque = zeros(num_cases, 1);
res_peak_stress = zeros(num_cases, 1); 
res_peak_force  = zeros(num_cases, 1);
res_is_broken   = false(num_cases, 1);
break_threshold = 58.9; % 击碎阈值 (N)

fprintf('==================================================\n');
fprintf('开始刚柔比自动寻优，共 %d 组参数...\n', num_cases);
fprintf('==================================================\n');

%% 2. 遍历每一个刚柔比进行仿真
for idx = 1:num_cases
    ratio = ratios(idx);
    
    % 计算当前的轮毂半径和绳子长度
    R_hub = total_length / (ratio + 1);
    L_rope = total_length - R_hub;
    
    % 动态计算节点数量，确保网格密度(l_bar)基本一致
    n_nodes_per_rod = max(5, round(L_rope / target_l_bar) + 1); 
    l_bar = L_rope / (n_nodes_per_rod - 1); % 重新计算真实的等分长度
    
    fprintf('\n[Case %d/%d] 刚柔比 %d:1 | 轮毂: %.3fm | 绳长: %.3fm | 节点数: %d\n', ...
        idx, num_cases, ratio, R_hub, L_rope, n_nodes_per_rod);
    
    % --------------------------------------------------------
    % 2.1 动态生成几何数据
    % --------------------------------------------------------
    nodes1 = zeros(n_nodes_per_rod, 3);
    nodes1(1, :) = [R_hub, 0, 0]; 
    
    % 生成具有微小初始弯曲的绳索形状（更容易启动旋转）
    for k = 2:n_nodes_per_rod
        sum_val = [0, 0, 0];
        for i = 0:(k-2)
            term_x = cos(pi/4 - (i*pi)/(4*(n_nodes_per_rod-1))) * cos((i*pi)/(n_nodes_per_rod-1));
            term_y = cos(pi/4 - (i*pi)/(4*(n_nodes_per_rod-1))) * sin((i*pi)/(n_nodes_per_rod-1));
            sum_val = sum_val + l_bar * [term_x, term_y, 0];
        end
        nodes1(k, :) = nodes1(1, :) + sum_val;
    end
    
    nodes2 = -nodes1; % 第二根绳子对称
    nodes = [nodes1; nodes2];
    
    edges1 = [(1:n_nodes_per_rod-1)', (2:n_nodes_per_rod)'];
    edges2 = [(n_nodes_per_rod+1:2*n_nodes_per_rod-1)', (n_nodes_per_rod+2:2*n_nodes_per_rod)'];
    edges = [edges1; edges2];
    face_nodes = [];
    fixed_node_indices = [1, n_nodes_per_rod + 1]; % 动态获取根部节点索引
    
    % --------------------------------------------------------
    % 2.2 配置物理与环境参数
    % --------------------------------------------------------
    sim_params.static_sim = false;
    sim_params.TwoDsim = false;
    sim_params.use_lineSearch = true;
    sim_params.use_midedge = false;
    sim_params.showFrames = false;
    
    sim_params.ramp_time = 0.5;
    RPM_target = 1000;
    sim_params.omega_target = (RPM_target * 2 * pi) / 60;
    sim_params.omega = [0; 0; 0];
    sim_params.hub_radius = R_hub; 
    
    % 补全容差参数以修复 timeStepper 报错
    sim_params.dt = 1e-4; 
    sim_params.totalTime = 1.0; 
    sim_params.tol = 1e-3;
    sim_params.ftol = 1e-3; 
    sim_params.dtol = 0.01;  
    sim_params.maximum_iter = 20;
    sim_params.log_data = false; 
    
    geom.rod_r0 = 0.002;
    geom.shell_h = 0;
    r = geom.rod_r0;
    geom.Axs = pi * r^2;
    geom.Ixs1 = (pi * r^4) / 4;
    geom.Ixs2 = (pi * r^4) / 4;
    geom.Jxs = (pi * r^4) / 2;
    
    % 补全材料与环境参数 
    material.youngs_shell = 0;  
    material.poisson_shell = 0; 
    material.density = 1100;
    material.youngs_rod = 1e9; 
    material.poisson_rod = 0.35;
    material.contact_stiffness = 50000;
    material.mu = 0.3;
    material.velTol = 1e-4; 
    
    env.ext_force_list = ["gravity", "viscous", "aerodynamic", "centrifugal", "coriolis", "selfContact"];
    env.g = [0; 0; -9.81];
    env.air_density = 1.225;
    env.Cd = 1.0;
    env.rho = 1.225; 
    env.eta = 0.1;
    
    % === [核心配置] 绕过底层校验机制的单冰柱配置 ===
    env.contact_params.ice_radius = 0.01;
    env.contact_params.num_ice = 1; 
    
    % 占位的假环形阵列参数，以骗过 createEnvironment 检查
    env.contact_params.array_radius = 0; 
    env.contact_params.array_center_dist = 0.18; 
    
    env.contact_params.omega_mag = sim_params.omega_target; 
    env.contact_params.omega_spin = 0; 
    env.contact_params.compute_friction = true; 
    env.contact_params.active_time = 0.5; 
    env.contact_params.sigma_t = 1.5e6;  
    env.contact_params.z_root = 0.07;     
    env.contact_params.rod_radius = geom.rod_r0;
    env.contact_params.is_broken = false(1, 1); 
    env.contact_params.peak_force = zeros(1, 1); 
    
    % --------------------------------------------------------
    % 2.3 初始化与仿真步进
    % --------------------------------------------------------
    [nodes, edges, rod_edges, shell_edges, rod_shell_joint_edges, rod_shell_joint_total_edges, face_nodes, face_edges, face_shell_edges, ...
        elStretchRod, elStretchShell, elBendRod, elBendSign, elBendShell, sign_faces, face_unit_norms] ...
        = createGeometry(nodes, edges, face_nodes);
    
    twist_angles = zeros(size(rod_edges,1)+size(rod_shell_joint_total_edges,1),1);
    
    % 打包生成 environment 和 imc (此时里面还没有 ice_positions)
    [environment, imc] = createEnvironmentAndIMCStructs(env, geom, material, sim_params);
    
    % 手动将真实的离散坐标强制写入 imc，让 IceContact.m 读取
    imc.ice_positions = [0.18, 0.0]; 
    imc.theta_accumulated = 0;
    
    softRobot = MultiRod(geom, material, twist_angles, nodes, edges, rod_edges, shell_edges, ...
        rod_shell_joint_edges, rod_shell_joint_total_edges, face_nodes, sign_faces, ...
        face_edges, face_shell_edges, sim_params, environment);
    softRobot.nodes_local = nodes;
    
    % 【修复对象转换为 double 报错】手动清理上一个循环遗留的对象数组
    clear stretch_springs bend_twist_springs hinge_springs triangle_springs;
    
    n_stretch = size(elStretchRod,1) + size(elStretchShell,1);
    n_bend_twist = size(elBendRod,1);
    
    if n_stretch > 0
        for s=1:n_stretch, stretch_springs(s) = stretchSpring(softRobot.refLen(s), elStretchRod(s,:), softRobot); end
    else, stretch_springs = []; end
    
    if n_bend_twist > 0
        for b=1:n_bend_twist, bend_twist_springs(b) = bendTwistSpring(elBendRod(b,:), elBendSign(b,:), [0 0], 0, softRobot); end
    else, bend_twist_springs = []; end
    
    hinge_springs = []; triangle_springs = [];
    
    softRobot = computeSpaceParallel(softRobot);
    theta = softRobot.q0(3*softRobot.n_nodes+1 : 3*softRobot.n_nodes+softRobot.n_edges_dof);
    [softRobot.m1, softRobot.m2] = computeMaterialDirectors(softRobot.a1, softRobot.a2, theta);
    bend_twist_springs = setkappa(softRobot, bend_twist_springs);
    softRobot.undef_refTwist = computeRefTwist_bend_twist_spring(bend_twist_springs, softRobot.a1, softRobot.tangent, zeros(n_bend_twist,1));
    softRobot.refTwist = computeRefTwist_bend_twist_spring(bend_twist_springs, softRobot.a1, softRobot.tangent, softRobot.undef_refTwist);
    
    softRobot.fixed_nodes = fixed_node_indices; 
    softRobot.fixed_edges = [];
    [softRobot.fixedDOF, softRobot.freeDOF] = FindFixedFreeDOF(softRobot.fixed_nodes, softRobot.fixed_edges, softRobot.n_DOF, softRobot.n_nodes);
    
    Nsteps = round(sim_params.totalTime/sim_params.dt);
    ctime = 0; 
    local_max_torque = 0; local_max_stress = 0;
    
    for timeStep = 1:Nsteps
        if ctime < sim_params.ramp_time
            current_omega_mag = (ctime / sim_params.ramp_time) * sim_params.omega_target;
        else
            current_omega_mag = sim_params.omega_target;
        end
        
        sim_params.omega = [0; 0; current_omega_mag];
        imc.omega_mag = current_omega_mag;            
        imc.theta_accumulated = imc.theta_accumulated + current_omega_mag * sim_params.dt;
        
        [softRobot, stretch_springs, bend_twist_springs, hinge_springs, force_now, imc] = ...
            timeStepper(softRobot, stretch_springs, bend_twist_springs, hinge_springs, ...
            triangle_springs, [], environment, imc, sim_params, ctime);
        
        if isfield(force_now, 'contact')
            F_resistive = force_now.drag + force_now.coriolis + force_now.contact; 
            F_res_vec = reshape(F_resistive(1:3*softRobot.n_nodes), 3, []);
            nodes_pos = reshape(softRobot.q(1:3*softRobot.n_nodes), 3, []);
            node_torques_z = nodes_pos(1, :) .* F_res_vec(2, :) - nodes_pos(2, :) .* F_res_vec(1, :);
            current_torque = abs(sum(node_torques_z));
            if current_torque > local_max_torque, local_max_torque = current_torque; end
        end
        
        % === 【核心修复】动态计算根部弯曲应力 ===
        % 1. 取出与第一根弯曲弹簧相关的边与节点索引
        edge1_idx = bend_twist_springs(1).edges_ind(1);
        edge2_idx = bend_twist_springs(1).edges_ind(2);
        node1_idx = bend_twist_springs(1).nodes_ind(1);
        node2_idx = bend_twist_springs(1).nodes_ind(2);
        node3_idx = bend_twist_springs(1).nodes_ind(3);

        % 2. 从软体机器人的全局状态 q 中提取连续 3 个节点的位置向量，转置为横向量(1x3)
        node1_loc = softRobot.q(3*node1_idx-2 : 3*node1_idx)';
        node2_loc = softRobot.q(3*node2_idx-2 : 3*node2_idx)';
        node3_loc = softRobot.q(3*node3_idx-2 : 3*node3_idx)';

        % 3. 提取对应的材料主轴向量 (考虑符号修正)
        m1e = softRobot.m1(edge1_idx, :);
        m2e = bend_twist_springs(1).sgn(1) * softRobot.m2(edge1_idx, :);
        m1f = softRobot.m1(edge2_idx, :);
        m2f = bend_twist_springs(1).sgn(2) * softRobot.m2(edge2_idx, :);

        % 4. 调用底层函数 computekappa 计算此刻真实曲率矢量
        kappa_root = computekappa(node1_loc, node2_loc, node3_loc, m1e, m2e, m1f, m2f);
        
        % 5. 求曲率模长并换算为 MPa 应力
        root_kappa_mag = norm(kappa_root); 
        current_stress_MPa = (material.youngs_rod * root_kappa_mag * geom.rod_r0) / 1e6; 
        if current_stress_MPa > local_max_stress, local_max_stress = current_stress_MPa; end
        % ======================================
        
        ctime = ctime + sim_params.dt;
        softRobot.q0 = softRobot.q; 
    end
    
    res_peak_torque(idx) = local_max_torque;
    res_peak_stress(idx) = local_max_stress;
    res_peak_force(idx)  = imc.peak_force(1);
    res_is_broken(idx)   = (imc.peak_force(1) >= break_threshold);
    
    fprintf(' -> 扭矩: %.3f N·m | 应力: %.1f MPa | 碰撞力: %.1f N | 击碎: %d\n', ...
        local_max_torque, local_max_stress, res_peak_force(idx), res_is_broken(idx));
end

%% 3. 计算综合代价函数 (Cost Function)
% 目标：在保证击碎冰柱 (约束) 的前提下，寻找扭矩和应力最小的平衡点

norm_torque = res_peak_torque / max(res_peak_torque);
norm_stress = res_peak_stress / max(res_peak_stress);

w_torque = 0.5; 
w_stress = 0.5;

cost_scores = zeros(num_cases, 1);
for i = 1:num_cases
    if ~res_is_broken(i)
        % 硬约束约束：未能击碎冰块的，直接淘汰
        cost_scores(i) = inf;
    else
        % 目标优化：计算综合代价 (越小越好)
        cost_scores(i) = w_torque * norm_torque(i) + w_stress * norm_stress(i);
    end
end

[min_cost, best_idx] = min(cost_scores);
best_ratio = ratios(best_idx);

%% 4. 可视化与综合报告
figure('Color', 'w', 'Position', [100, 100, 1200, 800]);

% 子图 1: 马达扭矩与弯曲应力
subplot(3, 1, 1);
yyaxis left;
plot(ratios, res_peak_torque, '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b');
ylabel('Motor Torque [N·m]');
ylim([0, max(res_peak_torque)*1.2]);
yyaxis right;
plot(ratios, res_peak_stress, '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
ylabel('Root Stress [MPa]');
ylim([0, max(res_peak_stress)*1.2]);
set(gca, 'XDir', 'reverse'); 
xticks(ratios); xticklabels(arrayfun(@(x) sprintf('%d:1', x), ratios, 'UniformOutput', false));
title('Raw Metrics: Motor Torque & Bending Stress'); grid on;

% 子图 2: 破冰碰撞力检测
subplot(3, 1, 2); hold on;
bar(ratios, res_peak_force, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'k');
yline(break_threshold, 'r--', 'LineWidth', 2, 'Label', sprintf('Break Threshold %.1fN', break_threshold));
set(gca, 'XDir', 'reverse');
xticks(ratios); xticklabels(arrayfun(@(x) sprintf('%d:1', x), ratios, 'UniformOutput', false));
ylabel('Collision Force [N]');
title('Constraint Check: Ice Breaking Capability'); grid on;

% 子图 3: 综合代价函数评估曲线
subplot(3, 1, 3); hold on;
plot_costs = cost_scores; plot_costs(isinf(plot_costs)) = NaN; 
plot(ratios, plot_costs, '-ko', 'LineWidth', 2, 'MarkerFaceColor', 'k', 'MarkerSize', 8);
if ~isinf(min_cost)
    plot(ratios(best_idx), min_cost, 'p', 'MarkerSize', 16, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k');
    text(ratios(best_idx), min_cost * 1.05, ' ★ Optimal', 'Color', 'g', 'FontWeight', 'bold', 'FontSize', 12);
end
set(gca, 'XDir', 'reverse');
xticks(ratios); xticklabels(arrayfun(@(x) sprintf('%d:1', x), ratios, 'UniformOutput', false));
xlabel('Rigidity-Flexibility Ratio (L_{rope} : R_{hub})');
ylabel('Cost Score (Lower is Better)');
title(sprintf('Comprehensive Evaluation (Weights: %.1f Torque, %.1f Stress)', w_torque, w_stress));
grid on; hold off;

% 打印总结报告
fprintf('\n================== 寻优结果总结 ==================\n');
fprintf('比例\t节点数\t扭矩(N·m)\t根部应力(MPa)\t碰撞力(N)\t状态\t综合代价\n');
for i = 1:num_cases
    statStr = '淘汰'; costStr = 'INF';
    if res_is_broken(i)
        statStr = '达标'; costStr = sprintf('%.3f', cost_scores(i));
    end
    n_nodes = max(5, round((total_length - total_length/(ratios(i)+1)) / target_l_bar) + 1);
    fprintf('%d:1\t%d\t%.4f\t\t%.2f\t\t%.2f\t\t%s\t%s\n', ...
        ratios(i), n_nodes, res_peak_torque(i), res_peak_stress(i), res_peak_force(i), statStr, costStr);
end
fprintf('==================================================\n');
if isinf(min_cost)
    fprintf('【警告】: 所有刚柔比配置均未能达到 58.9N 的击碎阈值，请考虑提高转速或选用更粗的绳索。\n');
else
    fprintf('>>> 自动寻优完成！推荐的最佳刚柔比为: 【 %d:1 】 <<<\n', best_ratio);
end