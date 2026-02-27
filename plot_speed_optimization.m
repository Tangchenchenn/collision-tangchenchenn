clc; clear all; close all;

% 路径设置
projectRoot = fileparts(mfilename('fullpath')); 
cd(projectRoot);
addpath(genpath(projectRoot));

%% 1. 定义要扫描的参数范围 (冰柱圈自转速度)
% 设定 5 个测试梯度 (rad/s)
test_spin_speeds = [0.5, 1.0, 1.5, 2.0, 2.5]; 
num_tests = length(test_spin_speeds);

% 用于存储从仿真中提取的真实数据
record_deicing_rate = zeros(1, num_tests);
record_peak_torque  = zeros(1, num_tests);

fprintf('开始批量仿真: 探究最合适自转速度 (共 %d 组)...\n', num_tests);

%% 2. 批量仿真循环
for i = 1:num_tests
    current_spin = test_spin_speeds(i);
    fprintf('\n>>> 开始第 %d 组测试: omega_spin = %.2f rad/s <<<\n', i, current_spin);
    
    % --- 加载配置 ---
    robotDescriptionDeicing; 
    
    % [核心覆盖] 覆盖配置中的自转速度参数
    env.contact_params.omega_spin = current_spin;
    
    % 强制参数修正
    if ~isfield(sim_params, 'omega'), sim_params.omega = [0; 0; sim_params.omega_target]; end
    env.ext_force_list = unique([env.ext_force_list, "centrifugal", "coriolis"]);
    
    % 创建几何与系统环境
    [nodes, edges, rod_edges, shell_edges, rod_shell_joint_edges, rod_shell_joint_total_edges, face_nodes, face_edges, face_shell_edges, ...
        elStretchRod, elStretchShell, elBendRod, elBendSign, elBendShell, sign_faces, face_unit_norms] = createGeometry(nodes, edges, face_nodes);
    twist_angles = zeros(size(rod_edges,1)+size(rod_shell_joint_total_edges,1),1);
    [environment, imc] = createEnvironmentAndIMCStructs(env, geom, material, sim_params);
    
    softRobot = MultiRod(geom, material, twist_angles, nodes, edges, rod_edges, shell_edges, ...
        rod_shell_joint_edges, rod_shell_joint_total_edges, face_nodes, sign_faces, face_edges, face_shell_edges, sim_params, environment);
    
    % 初始化弹簧
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
    
    % 初始化系统约束
    softRobot = computeSpaceParallel(softRobot);
    theta = softRobot.q0(3*softRobot.n_nodes+1 : 3*softRobot.n_nodes+softRobot.n_edges_dof);
    [softRobot.m1, softRobot.m2] = computeMaterialDirectors(softRobot.a1, softRobot.a2, theta);
    bend_twist_springs = setkappa(softRobot, bend_twist_springs);
    softRobot.undef_refTwist = computeRefTwist_bend_twist_spring(bend_twist_springs, softRobot.a1, softRobot.tangent, zeros(n_bend_twist,1));
    softRobot.refTwist = computeRefTwist_bend_twist_spring(bend_twist_springs, softRobot.a1, softRobot.tangent, softRobot.undef_refTwist);
    softRobot.fixed_nodes = fixed_node_indices; softRobot.fixed_edges = [];
    [softRobot.fixedDOF, softRobot.freeDOF] = FindFixedFreeDOF(softRobot.fixed_nodes, softRobot.fixed_edges, softRobot.n_DOF, softRobot.n_nodes);
    
    % --- 物理仿真循环 ---
    Nsteps = round(sim_params.totalTime/sim_params.dt);
    ctime = 0; 
    F_history_motor_torque = zeros(Nsteps,1);
    break_times = inf(1, env.contact_params.num_ice);
    
    for timeStep = 1:Nsteps
        if ctime < sim_params.ramp_time
            current_omega_mag = (ctime / sim_params.ramp_time) * sim_params.omega_target;
        else
            current_omega_mag = sim_params.omega_target;
        end
        sim_params.omega = [0; 0; current_omega_mag]; 
        imc.omega_mag = current_omega_mag;            
        imc.theta_accumulated = imc.theta_accumulated + current_omega_mag * sim_params.dt;
        
        tau_0 = [];
        [softRobot, stretch_springs, bend_twist_springs, hinge_springs, force_now, imc] = ...
            timeStepper(softRobot, stretch_springs, bend_twist_springs, hinge_springs, ...
            triangle_springs, tau_0, environment, imc, sim_params, ctime);
        
        % 记录断裂状态
        if any(imc.is_broken)
            newly_broken = imc.is_broken & isinf(break_times);
            break_times(newly_broken) = ctime;
        end
        
        % 记录马达扭矩
        if isfield(force_now, 'contact')
            F_resistive = force_now.drag + force_now.coriolis + force_now.contact; 
            F_res_vec = reshape(F_resistive(1:3*softRobot.n_nodes), 3, []);
            nodes_pos = reshape(softRobot.q(1:3*softRobot.n_nodes), 3, []);
            X = nodes_pos(1, :); Y = nodes_pos(2, :);
            Fx = F_res_vec(1, :); Fy = F_res_vec(2, :);
            F_history_motor_torque(timeStep) = abs(sum(X .* Fy - Y .* Fx)); 
        end
        
        ctime = ctime + sim_params.dt;
        softRobot.q0 = softRobot.q; 
    end
    
    % --- 提取并保存当前仿真评价指标 ---
    broken_count = sum(imc.is_broken);
    record_deicing_rate(i) = (broken_count / imc.num_ice) * 100;
    record_peak_torque(i)  = max(F_history_motor_torque);
    
    fprintf('本组测试结束: 除冰率 = %.1f%%, 峰值扭矩 = %.2f N.m\n', ...
        record_deicing_rate(i), record_peak_torque(i));
end

%% 3. 绘制最终论文用图 (双Y轴: 速度 vs 除冰率与扭矩)
figure('Name', 'Speed Optimization', 'Color', 'w', 'Position', [100, 100, 600, 450]);

yyaxis left
plot(test_spin_speeds, record_deicing_rate, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'Color', '#0072BD');
ylabel('Deicing Rate (%)', 'FontWeight', 'bold');
ylim([0, 110]); % 留点裕度

yyaxis right
plot(test_spin_speeds, record_peak_torque, '-s', 'LineWidth', 2, 'MarkerSize', 8, 'Color', '#D95319');
ylabel('Peak Motor Torque (N·m)', 'FontWeight', 'bold');
% ylim 根据你的实际最大扭矩自动调整即可

xlabel('Icicle Array Spin Speed \omega_{spin} (rad/s)', 'FontWeight', 'bold');
title('Effect of Spin Speed on Deicing Performance', 'FontSize', 12);
grid on;
ax = gca; ax.FontSize = 11;
legend({'Deicing Rate', 'Peak Motor Torque'}, 'Location', 'best');