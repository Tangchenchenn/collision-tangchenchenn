clc; clear all; close all;

projectRoot = fileparts(mfilename('fullpath')); 
cd(projectRoot); addpath(genpath(projectRoot));

%% 1. 定义要扫描的参数范围 (除冰距离)
% 假设轮毂半径0.05，绳长0.15，最大旋转半径0.20。
% 测试距离从 0.16(过近) 到 0.22(过远，可能打不到)
test_distances = [0.16, 0.17, 0.18, 0.19, 0.20, 0.21]; 
num_tests = length(test_distances);

% 用于存储真实数据
record_mean_peak_force = zeros(1, num_tests); % 记录有效击打冰柱的平均峰值力
record_max_peak_force  = zeros(1, num_tests); % 记录最大绝对峰值力

fprintf('开始批量仿真: 探究最合适除冰距离 (共 %d 组)...\n', num_tests);

%% 2. 批量仿真循环
for i = 1:num_tests
    current_dist = test_distances(i);
    fprintf('\n>>> 开始第 %d 组测试: array_center_dist = %.2f m <<<\n', i, current_dist);
    
    % --- 加载配置 ---
    robotDescriptionDeicing; 
    
    % [核心覆盖] 覆盖除冰距离参数
    env.contact_params.array_center_dist = current_dist;
    
    if ~isfield(sim_params, 'omega'), sim_params.omega = [0; 0; sim_params.omega_target]; end
    env.ext_force_list = unique([env.ext_force_list, "centrifugal", "coriolis"]);
    
    % 几何与环境初始化
    [nodes, edges, rod_edges, shell_edges, rod_shell_joint_edges, rod_shell_joint_total_edges, face_nodes, face_edges, face_shell_edges, ...
        elStretchRod, elStretchShell, elBendRod, elBendSign, elBendShell, sign_faces, face_unit_norms] = createGeometry(nodes, edges, face_nodes);
    twist_angles = zeros(size(rod_edges,1)+size(rod_shell_joint_total_edges,1),1);
    [environment, imc] = createEnvironmentAndIMCStructs(env, geom, material, sim_params);
    
    softRobot = MultiRod(geom, material, twist_angles, nodes, edges, rod_edges, shell_edges, ...
        rod_shell_joint_edges, rod_shell_joint_total_edges, face_nodes, sign_faces, face_edges, face_shell_edges, sim_params, environment);
    
    n_stretch = size(elStretchRod,1) + size(elStretchShell,1); n_bend_twist = size(elBendRod,1);
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
            
        ctime = ctime + sim_params.dt;
        softRobot.q0 = softRobot.q; 
    end
    
    % --- 提取数据 ---
    % imc.peak_force 记录了仿真中每根冰柱承受的峰值力
    valid_forces = imc.peak_force(imc.peak_force > 0); % 剔除根本没被打到的冰柱
    
    if isempty(valid_forces)
        record_mean_peak_force(i) = 0;
        record_max_peak_force(i) = 0;
    else
        record_mean_peak_force(i) = mean(valid_forces);
        record_max_peak_force(i)  = max(valid_forces);
    end
    
    fprintf('测试结束: 平均峰值力 = %.2f N, 绝对最大峰值力 = %.2f N\n', ...
        record_mean_peak_force(i), record_max_peak_force(i));
end

%% 3. 绘制最终论文用图 (单Y轴: 距离 vs 峰值接触力)
figure('Name', 'Distance Optimization', 'Color', 'w', 'Position', [750, 100, 600, 450]);

plot(test_distances, record_mean_peak_force, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'Color', '#77AC30', 'DisplayName', 'Mean Peak Force');
hold on;
plot(test_distances, record_max_peak_force, '--^', 'LineWidth', 2, 'MarkerSize', 8, 'Color', '#A2142F', 'DisplayName', 'Max Peak Force');

% [关键] 绘制冰柱屈服断裂所需的理论力阈值线，假设需要 50N (具体按你的 sigma_t 算出的值修改)
break_threshold = 50; 
yline(break_threshold, 'k:', 'LineWidth', 2, 'DisplayName', 'Fracture Threshold');

xlabel('Deicing Distance from Center (m)', 'FontWeight', 'bold');
ylabel('Impact Peak Force (N)', 'FontWeight', 'bold');
title('Effect of Deicing Distance on Impact Force', 'FontSize', 12);
grid on;
ax = gca; ax.FontSize = 11;
legend('Location', 'best');