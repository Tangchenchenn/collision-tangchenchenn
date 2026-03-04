clc; clear all; close all;

% 路径设置
projectRoot = fileparts(mfilename('fullpath')); 
cd(projectRoot);
addpath(genpath(projectRoot));

%% 1. 寻优参数设置
total_length = 0.262; % 单侧绳子长度 + 轮毂半径 = 0.262m
ratios = 9:-1:1;      % 刚柔比 (绳长 : 轮毂半径) 9:1 到 1:1
num_cases = length(ratios);

% === 【核心修复1：时间参数修正】 ===
% 修复原代码中 totalTime < stable_time 导致接触永远不触发的Bug
sim_params.ramp_time = 0.5;            % 马达加速时间
sim_params.stable_time = 0.8;          % 稳定运行 0.8 秒后再激活冰柱
sim_params.totalTime = 5;            % 仿真总时长 1.5 秒 (留出 0.7s 的接触时间)
% ====================================

% 网格控制参数
target_l_bar = 0.015; % 目标网格单元长度约为 1.5 cm

% === 【核心修复2：评估指标修改为 RMS 扭矩与击碎率】 ===
num_ice_total = 43;          % 环形阵列的冰柱总数
target_break_rate = 0.9;     % 目标击碎率：必须击碎至少 80% 的冰柱才算达标

res_rms_torque  = zeros(num_cases, 1); % 替换原来的峰值扭矩
res_peak_stress = zeros(num_cases, 1); % 根部最大应力（弯曲疲劳考察）
res_break_rate  = zeros(num_cases, 1); % 击碎率
res_is_valid    = false(num_cases, 1); % 是否满足目标击碎率硬约束
% ======================================================

fprintf('==================================================\n');
fprintf('开始刚柔比自动寻优(环形阵列 + RMS均方根评估)，共 %d 组...\n', num_cases);
fprintf('==================================================\n');

%% 2. 遍历每一个刚柔比进行仿真
for idx = 1:num_cases
    ratio = ratios(idx);

    R_hub = total_length / (ratio + 1);
    L_rope = total_length - R_hub;

    n_nodes_per_rod = max(5, round(L_rope / target_l_bar) + 1); 
    l_bar = L_rope / (n_nodes_per_rod - 1); 

    fprintf('\n[Case %d/%d] 刚柔比 %d:1 | 轮毂: %.3fm | 绳长: %.3fm | 节点数: %d\n', ...
        idx, num_cases, ratio, R_hub, L_rope, n_nodes_per_rod);

    % --------------------------------------------------------
    % 2.1 动态生成几何数据 (微小初始弯曲)
    % --------------------------------------------------------
    nodes1 = zeros(n_nodes_per_rod, 3);
    nodes1(1, :) = [R_hub, 0, 0]; 
    for k = 2:n_nodes_per_rod
        sum_val = [0, 0, 0];
        for i = 0:(k-2)
            term_x = cos(pi/4 - (i*pi)/(4*(n_nodes_per_rod-1))) * cos((i*pi)/(n_nodes_per_rod-1));
            term_y = cos(pi/4 - (i*pi)/(4*(n_nodes_per_rod-1))) * sin((i*pi)/(n_nodes_per_rod-1));
            sum_val = sum_val + l_bar * [term_x, term_y, 0];
        end
        nodes1(k, :) = nodes1(1, :) + sum_val;
    end
    nodes2 = -nodes1; 
    nodes = [nodes1; nodes2];

    edges1 = [(1:n_nodes_per_rod-1)', (2:n_nodes_per_rod)'];
    edges2 = [(n_nodes_per_rod+1:2*n_nodes_per_rod-1)', (n_nodes_per_rod+2:2*n_nodes_per_rod)'];
    edges = [edges1; edges2];
    face_nodes = [];
    fixed_node_indices = [1, n_nodes_per_rod + 1]; 

    % --------------------------------------------------------
    % 2.2 配置物理与环境参数
    % --------------------------------------------------------
    sim_params.static_sim = false;
    sim_params.TwoDsim = false;
    sim_params.use_lineSearch = true;
    sim_params.use_midedge = false;
    sim_params.showFrames = false;

    RPM_target = 1000;
    sim_params.omega_target = (RPM_target * 2 * pi) / 60;
    sim_params.omega = [0; 0; 0];
    sim_params.hub_radius = R_hub; 

    sim_params.dt = 1e-4; 
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

    % === 【核心修复3：恢复环形阵列接触配置】 ===
    env.contact_params.ice_radius = 0.007;
    env.contact_params.num_ice = num_ice_total; % 40根冰柱
    env.contact_params.array_radius = 0.145;    % 圆周半径
    env.contact_params.array_center_dist = 0.25;% 圆心距中心距离
    
    env.contact_params.omega_mag = sim_params.omega_target; 
    env.contact_params.omega_spin = 1.5; 
    env.contact_params.compute_friction = true; 
    env.contact_params.active_time = sim_params.stable_time; % 稳定后再检测
    env.contact_params.sigma_t = 1.5e6;  
    env.contact_params.z_root = 0.02;     
    env.contact_params.ice_length = 0.040; 
    env.contact_params.ice_z_min = env.contact_params.z_root - env.contact_params.ice_length; 
    env.contact_params.ice_z_max = env.contact_params.z_root; 
    
    env.contact_params.rod_radius = geom.rod_r0;
    env.contact_params.is_broken = false(1, num_ice_total); 
    env.contact_params.peak_force = zeros(1, num_ice_total); 
    % =================================================

    % --------------------------------------------------------
    % 2.3 初始化与仿真步进
    % --------------------------------------------------------
    [nodes, edges, rod_edges, shell_edges, rod_shell_joint_edges, rod_shell_joint_total_edges, face_nodes, face_edges, face_shell_edges, ...
        elStretchRod, elStretchShell, elBendRod, elBendSign, elBendShell, sign_faces, face_unit_norms] ...
        = createGeometry(nodes, edges, face_nodes);

    twist_angles = zeros(size(rod_edges,1)+size(rod_shell_joint_total_edges,1),1);
    [environment, imc] = createEnvironmentAndIMCStructs(env, geom, material, sim_params);
    
    % 删除原来强制写入 ice_positions 的占位符，让底层自动生成环形坐标
    imc.theta_accumulated = 0;

    softRobot = MultiRod(geom, material, twist_angles, nodes, edges, rod_edges, shell_edges, ...
        rod_shell_joint_edges, rod_shell_joint_total_edges, face_nodes, sign_faces, ...
        face_edges, face_shell_edges, sim_params, environment);
    softRobot.nodes_local = nodes;

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
    
    local_max_stress = 0;
    active_torques = []; % 用于记录激活期间的所有扭矩样本

    for timeStep = 1:Nsteps
        if ctime < sim_params.ramp_time
            current_omega_mag = (ctime / sim_params.ramp_time) * sim_params.omega_target;
        else
            current_omega_mag = sim_params.omega_target;
        end
        sim_params.omega = [0; 0; current_omega_mag];
        imc.omega_mag = current_omega_mag;

        if ctime < sim_params.stable_time
            environment.selfContact = false;   
            imc.active = false;                
            imc.is_active = false;             
        else
            environment.selfContact = true;
            imc.active = true;
            imc.is_active = true;
        end

        imc.theta_accumulated = imc.theta_accumulated + current_omega_mag * sim_params.dt;

        [softRobot, stretch_springs, bend_twist_springs, hinge_springs, force_now, imc] = ...
            timeStepper(softRobot, stretch_springs, bend_twist_springs, hinge_springs, ...
                        triangle_springs, [], environment, imc, sim_params, ctime);

       % === [方案 A] 同步修改：提取根部单元内部扭转力矩 (Transmitted Torque) ===
        % 逻辑：计算由于碰撞导致的应力波传回马达轴心的真实力矩
        n_nodes = softRobot.n_nodes;
        twist_offset = 3 * n_nodes; 
        
        % 动态计算每根绳索的边缘数量（用于处理不同网格密度）
        num_edges_per_rod = softRobot.n_edges_rod_only / 2;
        root_edge_indices = [1, num_edges_per_rod + 1];
        
        total_internal_torque = 0;
        for e_idx = root_edge_indices
            if (twist_offset + e_idx) <= length(softRobot.q) && e_idx <= length(softRobot.refTwist)
                % 1. 提取当前根部扭转角 theta
                theta_current = softRobot.q(twist_offset + e_idx);
                % 2. 计算扭转偏差 (减去马达端固定角 0 和参考扭率)
                delta_theta = theta_current - 0 - softRobot.refTwist(e_idx);
                % 3. M = (GJ / L) * delta_theta (注意：GJ 是标量)
                torque_i = (softRobot.GJ / softRobot.refLen(e_idx)) * delta_theta;
                
                total_internal_torque = total_internal_torque + abs(torque_i);
            end
        end
        current_torque = total_internal_torque;
        
        % 只有当达到稳定时间(冰柱激活)后，才收集扭矩样本以供后续计算 RMS
        if ctime >= sim_params.stable_time
            active_torques = [active_torques, current_torque];
        end

        % --- 计算根部弯曲应力不变 ---
        edge1_idx = bend_twist_springs(1).edges_ind(1);
        edge2_idx = bend_twist_springs(1).edges_ind(2);
        node1_idx = bend_twist_springs(1).nodes_ind(1);
        node2_idx = bend_twist_springs(1).nodes_ind(2);
        node3_idx = bend_twist_springs(1).nodes_ind(3);

        node1_loc = softRobot.q(3*node1_idx-2 : 3*node1_idx)';
        node2_loc = softRobot.q(3*node2_idx-2 : 3*node2_idx)';
        node3_loc = softRobot.q(3*node3_idx-2 : 3*node3_idx)';

        m1e = softRobot.m1(edge1_idx, :);
        m2e = bend_twist_springs(1).sgn(1) * softRobot.m2(edge1_idx, :);
        m1f = softRobot.m1(edge2_idx, :);
        m2f = bend_twist_springs(1).sgn(2) * softRobot.m2(edge2_idx, :);

        kappa_root = computekappa(node1_loc, node2_loc, node3_loc, m1e, m2e, m1f, m2f);
        root_kappa_mag = norm(kappa_root); 
        current_stress_MPa = (material.youngs_rod * root_kappa_mag * geom.rod_r0) / 1e6; 
        if current_stress_MPa > local_max_stress, local_max_stress = current_stress_MPa; end

        ctime = ctime + sim_params.dt;
        softRobot.q0 = softRobot.q; 
    end

    % === 【核心修复5：结算 RMS 扭矩与击碎率】 ===
    if ~isempty(active_torques)
        res_rms_torque(idx) = sqrt(mean(active_torques.^2)); % 计算均方根扭矩
    else
        res_rms_torque(idx) = 0; % 完美闪避，没有数据
    end
    
    res_peak_stress(idx) = local_max_stress;
    
    % 统计击碎了多少根冰柱
    broken_count = sum(imc.is_broken);
    res_break_rate(idx) = broken_count / num_ice_total; 
    res_is_valid(idx)   = (res_break_rate(idx) >= target_break_rate); % 是否达标

    fprintf(' -> 均方根扭矩: %.3f N·m | 最大应力: %.1f MPa | 击碎率: %.1f%% (%d/%d) | 达标: %d\n', ...
        res_rms_torque(idx), res_peak_stress(idx), res_break_rate(idx)*100, broken_count, num_ice_total, res_is_valid(idx));
    % ====================================================
end

%% 3. 计算综合代价函数 (Cost Function)
norm_torque = res_rms_torque / max(res_rms_torque);
norm_stress = res_peak_stress / max(res_peak_stress);

w_torque = 0.5; 
w_stress = 0.5;

cost_scores = zeros(num_cases, 1);
for i = 1:num_cases
    if ~res_is_valid(i)
        % 硬约束：未能达到 80% 击碎率的刚柔比，直接淘汰
        cost_scores(i) = inf;
    else
        % 目标优化：计算综合代价 (扭矩越小越好，应力越小越好)
        cost_scores(i) = w_torque * norm_torque(i) + w_stress * norm_stress(i);
    end
end

[min_cost, best_idx] = min(cost_scores);
best_ratio = ratios(best_idx);

%% 4. 可视化与综合报告
figure('Color', 'w', 'Position', [100, 100, 1200, 800]);

% 子图 1: RMS 扭矩与弯曲应力
subplot(3, 1, 1);
yyaxis left;
plot(ratios, res_rms_torque, '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b');
ylabel('Motor RMS Torque [N·m]');
ylim([0, max(res_rms_torque)*1.2 + 0.01]);
yyaxis right;
plot(ratios, res_peak_stress, '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r');
ylabel('Root Peak Stress [MPa]');
ylim([0, max(res_peak_stress)*1.2 + 0.1]);
set(gca, 'XDir', 'reverse'); 
sorted_ratios = sort(ratios); 
xticks(sorted_ratios); 
xticklabels(arrayfun(@(x) sprintf('%d:1', x), sorted_ratios, 'UniformOutput', false));
title('Raw Metrics: Steady-State RMS Torque & Peak Bending Stress'); grid on;

% 子图 2: 击碎率检测 (Breakage Rate)
subplot(3, 1, 2); hold on;
bar(ratios, res_break_rate * 100, 'FaceColor', [0.3 0.7 0.9], 'EdgeColor', 'k');
yline(target_break_rate * 100, 'r--', 'LineWidth', 2, 'Label', sprintf('Target Rate: %.0f%%', target_break_rate * 100));
set(gca, 'XDir', 'reverse');
xticks(sorted_ratios); 
xticklabels(arrayfun(@(x) sprintf('%d:1', x), sorted_ratios, 'UniformOutput', false));
ylabel('Ice Breakage Rate [%]');
ylim([0, 110]);
title('Constraint Check: Ice Breakage Rate (Ring Array)'); grid on;

% 子图 3: 综合代价函数评估曲线
subplot(3, 1, 3); hold on;
plot_costs = cost_scores; plot_costs(isinf(plot_costs)) = NaN; 
plot(ratios, plot_costs, '-ko', 'LineWidth', 2, 'MarkerFaceColor', 'k', 'MarkerSize', 8);
if ~isinf(min_cost)
    plot(ratios(best_idx), min_cost, 'p', 'MarkerSize', 16, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k');
    text(ratios(best_idx), min_cost * 1.05, ' ★ Optimal', 'Color', 'g', 'FontWeight', 'bold', 'FontSize', 12);
end
set(gca, 'XDir', 'reverse');
xticks(sorted_ratios); 
xticklabels(arrayfun(@(x) sprintf('%d:1', x), sorted_ratios, 'UniformOutput', false));
xlabel('Rigidity-Flexibility Ratio (L_{rope} : R_{hub})');
ylabel('Cost Score (Lower is Better)');
title(sprintf('Comprehensive Evaluation (Weights: %.1f Torque, %.1f Stress)', w_torque, w_stress));
grid on; hold off;

% 打印总结报告
fprintf('\n================== 寻优结果总结 ==================\n');
fprintf('比例\tRMS扭矩(N·m)\t根部应力(MPa)\t击碎率(%%)\t状态\t综合代价\n');
for i = 1:num_cases
    statStr = '淘汰'; costStr = 'INF';
    if res_is_valid(i)
        statStr = '达标'; costStr = sprintf('%.3f', cost_scores(i));
    end
    fprintf('%d:1\t%.4f\t\t%.2f\t\t%.1f\t\t%s\t%s\n', ...
        ratios(i), res_rms_torque(i), res_peak_stress(i), res_break_rate(i)*100, statStr, costStr);
end
fprintf('==================================================\n');
if isinf(min_cost)
    fprintf('【警告】: 所有刚柔比配置均未能达到 %.0f%% 的击碎率约束，请考虑提高转速或选用更粗的绳索。\n', target_break_rate*100);
else
    fprintf('>>> 自动寻优完成！推荐的最佳刚柔比为: 【 %d:1 】 <<<\n', best_ratio);
end