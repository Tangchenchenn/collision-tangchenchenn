clc; clearvars; close all;

% 路径设置
projectRoot = fileparts(mfilename('fullpath')); 
cd(projectRoot);
addpath(genpath(projectRoot));

%% 1. 加载配置
robotDescriptionDeicing; 

% 强制参数修正
if ~isfield(sim_params, 'omega'), sim_params.omega = [0; 0; sim_params.omega_target]; end
env.ext_force_list = unique([env.ext_force_list, "centrifugal", "coriolis"]);

% 保存初始坐标
nodes_original = nodes;

% 创建几何与对象
[nodes, edges, rod_edges, shell_edges, rod_shell_joint_edges, rod_shell_joint_total_edges, face_nodes, face_edges, face_shell_edges, ...
    elStretchRod, elStretchShell, elBendRod, elBendSign, elBendShell, sign_faces, face_unit_norms] ...
    = createGeometry(nodes, edges, face_nodes);

twist_angles = zeros(size(rod_edges,1)+size(rod_shell_joint_total_edges,1),1);
[environment, imc] = createEnvironmentAndIMCStructs(env, geom, material, sim_params);

softRobot = MultiRod(geom, material, twist_angles, nodes, edges, rod_edges, shell_edges, ...
    rod_shell_joint_edges, rod_shell_joint_total_edges, face_nodes, sign_faces, ...
    face_edges, face_shell_edges, sim_params, environment);
softRobot.nodes_local = nodes_original; 

%% 创建弹簧
n_stretch = size(elStretchRod,1) + size(elStretchShell,1);
n_bend_twist = size(elBendRod,1);

if n_stretch==0
    stretch_springs=[]; 
else
    % 【修正】：恢复为正序循环，避免触发 MATLAB 对象数组隐式调用 0 参数构造函数的报错机制
    for s = 1:n_stretch
        if s <= size(elStretchRod,1)
            stretch_springs(s) = stretchSpring(softRobot.refLen(s), elStretchRod(s,:), softRobot);
        else
            stretch_springs(s) = stretchSpring(softRobot.refLen(s), elStretchShell(s-size(elStretchRod,1),:), softRobot, softRobot.ks(s)); 
        end
    end
end

if n_bend_twist==0
    bend_twist_springs=[]; 
else
    % 【修正】：恢复为正序循环
    for b = 1:n_bend_twist
        bend_twist_springs(b) = bendTwistSpring(elBendRod(b,:), elBendSign(b,:), [0 0], 0, softRobot);
    end
end

hinge_springs = []; triangle_springs = [];

%% 初始化系统
softRobot = computeSpaceParallel(softRobot);
theta = softRobot.q0(3*softRobot.n_nodes+1 : 3*softRobot.n_nodes+softRobot.n_edges_dof);
theta = theta(:); % 强制转换为列向量，防止后续运算发生维数不匹配

[softRobot.m1, softRobot.m2] = computeMaterialDirectors(softRobot.a1, softRobot.a2, theta);
bend_twist_springs = setkappa(softRobot, bend_twist_springs);
softRobot.undef_refTwist = computeRefTwist_bend_twist_spring(bend_twist_springs, softRobot.a1, softRobot.tangent, zeros(n_bend_twist,1));
softRobot.refTwist = computeRefTwist_bend_twist_spring(bend_twist_springs, softRobot.a1, softRobot.tangent, softRobot.undef_refTwist);

% 边界条件
softRobot.fixed_nodes = fixed_node_indices; 
softRobot.fixed_edges = [];

[softRobot.fixedDOF, softRobot.freeDOF] = FindFixedFreeDOF(softRobot.fixed_nodes, softRobot.fixed_edges, softRobot.n_DOF, softRobot.n_nodes);

%% ==========================================
%% 2. 物理仿真循环 (Run Simulation)
%% ==========================================
Nsteps = round(sim_params.totalTime/sim_params.dt);
ctime = 0; 
time_arr = linspace(0, sim_params.totalTime, Nsteps);
sim_params.log_data = true; 
sim_params.logStep = 10;   

dof_with_time = zeros(softRobot.n_DOF+1, Nsteps);
dof_with_time(1,:) = time_arr;

% 力和扭矩记录初始化
F_history = struct('stretch', zeros(Nsteps,1), 'bend', zeros(Nsteps,1), ...
                   'twist', zeros(Nsteps,1), 'cent', zeros(Nsteps,1), ...
                   'coriolis', zeros(Nsteps,1), 'contact', zeros(Nsteps,1), ...
                   'motor_torque', zeros(Nsteps,1), ...
                   'motor_omega', zeros(Nsteps,1), ...
                   'motor_power', zeros(Nsteps,1)); 

fprintf('开始物理仿真 (共 %d 步)...\n', Nsteps);

% --- [循环前初始化] ---
break_times = inf(1, imc.num_ice); % 记录每根冰柱的断裂时间

% 记录上一时刻的峰值力和碰撞位置
previous_peak_forces = zeros(1, imc.num_ice); 
collision_dists = zeros(1, imc.num_ice); 

% 定义冰柱出现的“稳定时间”，设定为马达加速时间加上一点缓冲（例如 0.5s 加速 + 0.1s 稳定）
% 定义冰柱出现的“稳定时间”
if isfield(env, 'contact_params') && isfield(env.contact_params, 'active_time')
    t_stable = env.contact_params.active_time; % 优先使用配置中的 0.8s
else
    t_stable = sim_params.ramp_time + 0.1;     % 备用方案
end

for timeStep = 1:Nsteps
    % 1. 马达转速控制
    if ctime < sim_params.ramp_time
        ratio = ctime / sim_params.ramp_time;
        current_omega_mag = ratio * sim_params.omega_target;
    else
        current_omega_mag = sim_params.omega_target;
    end
    
    sim_params.omega = [0; 0; current_omega_mag]; 
    imc.omega_mag = current_omega_mag;            
    imc.theta_accumulated = imc.theta_accumulated + current_omega_mag * sim_params.dt;
    
    if(sim_params.use_midedge), tau_0 = updatePreComp_without_sign(softRobot.q, softRobot); else, tau_0 = []; end
    
    % 2. 核心物理步进计算
    [softRobot, stretch_springs, bend_twist_springs, hinge_springs, force_now, imc] = ...
        timeStepper(softRobot, stretch_springs, bend_twist_springs, hinge_springs, ...
        triangle_springs, tau_0, environment, imc, sim_params, ctime);

    % --- [修正] 移除了原有的延时触发屏蔽逻辑，确保碰撞力实时反馈 ---

    % 3. 追踪峰值接触力与碰撞位置
    is_colliding_now = imc.peak_force > previous_peak_forces;
    actual_peak_force_this_step = sum(imc.peak_force(is_colliding_now));
    
    if any(is_colliding_now) && isfield(force_now, 'contact')
        contact_forces_vec = reshape(force_now.contact(1:3*softRobot.n_nodes), 3, []);
        node_contact_mags = vecnorm(contact_forces_vec);
        [~, max_node] = max(node_contact_mags);
        
        nodes_pos = reshape(softRobot.q(1:3*softRobot.n_nodes), 3, []);
        hit_z = abs(nodes_pos(3, max_node)); 
        collision_dists(is_colliding_now) = hit_z;
    end
    
    previous_peak_forces = imc.peak_force; 

    % 4. 断裂事件捕获
    if any(imc.is_broken)
        newly_broken = imc.is_broken & isinf(break_times);
        if any(newly_broken)
            broken_indices = find(newly_broken);
            for idx = broken_indices
                break_times(idx) = ctime;
                fprintf('>>> 冰柱 #%02d 断裂! 时刻: %.4fs | 峰值力: %6.1f N | 距根部距离: %.3f m\n', ...
                    idx, ctime, imc.peak_force(idx), collision_dists(idx));
            end
        end
    end

    % 5. 接触力与扭矩记录 [修正版]
    if isfield(force_now, 'contact')
        % 包含接触力的阻力总和
        F_resistive = force_now.drag + force_now.coriolis + force_now.contact; 
        F_res_vec = reshape(F_resistive(1:3*softRobot.n_nodes), 3, []);
        
        nodes_pos = reshape(softRobot.q(1:3*softRobot.n_nodes), 3, []);
        X = nodes_pos(1, :);
        Y = nodes_pos(2, :);
        Fx = F_res_vec(1, :);
        Fy = F_res_vec(2, :);
        
        % 计算总驱动扭矩
        node_torques_z = X .* Fy - Y .* Fx;
        F_history.motor_torque(timeStep) = abs(sum(node_torques_z)); 
        
        % 记录接触力图表数据
        contact_forces_vec = reshape(force_now.contact(1:3*softRobot.n_nodes), 3, []);
        total_contact_force_mag = sum(vecnorm(contact_forces_vec));
        F_history.contact(timeStep) = max(total_contact_force_mag, actual_peak_force_this_step);
    end
    
    % 6. 记录马达转速与计算功率
    F_history.motor_omega(timeStep) = current_omega_mag;
    F_history.motor_power(timeStep) = F_history.motor_torque(timeStep) * current_omega_mag;
    
    ctime = ctime + sim_params.dt;
    softRobot.q0 = softRobot.q; 
    
    if sim_params.log_data && mod(timeStep, sim_params.logStep) == 0
        dof_with_time(2:end,timeStep) = softRobot.q;
    end
end

% for timeStep = 1:Nsteps
%     % 1. 马达转速控制
%     if ctime < sim_params.ramp_time
%         ratio = ctime / sim_params.ramp_time;
%         current_omega_mag = ratio * sim_params.omega_target;
%     else
%         current_omega_mag = sim_params.omega_target;
%     end
% 
%     sim_params.omega = [0; 0; current_omega_mag]; 
%     imc.omega_mag = current_omega_mag;            
%     imc.theta_accumulated = imc.theta_accumulated + current_omega_mag * sim_params.dt;
% 
%     if(sim_params.use_midedge), tau_0 = updatePreComp_without_sign(softRobot.q, softRobot); else, tau_0 = []; end
% 
%     % 2. 核心物理步进计算
%     [softRobot, stretch_springs, bend_twist_springs, hinge_springs, force_now, imc] = ...
%         timeStepper(softRobot, stretch_springs, bend_twist_springs, hinge_springs, ...
%         triangle_springs, tau_0, environment, imc, sim_params, ctime);
% 
%     % ========================================================
%     % 延时触发机制：在绳索稳定前，强制屏蔽接触力
%     % % ========================================================
%     % if ctime < t_stable
%     %     imc.peak_force(:) = 0;        % 强制重置底层记录的接触峰值
%     %     imc.is_broken(:) = false;     % 强制防止冰柱被早期刮擦判定为断裂
%     %     if isfield(force_now, 'contact')
%     %         force_now.contact(:) = 0; % 强制接触力数组归零
%     %     end
%     % end
% 
%     % 3. 追踪峰值接触力与碰撞位置
%     is_colliding_now = imc.peak_force > previous_peak_forces;
%     actual_peak_force_this_step = sum(imc.peak_force(is_colliding_now));
% 
%     if any(is_colliding_now) && isfield(force_now, 'contact')
%         contact_forces_vec = reshape(force_now.contact(1:3*softRobot.n_nodes), 3, []);
%         node_contact_mags = vecnorm(contact_forces_vec);
%         [~, max_node] = max(node_contact_mags);
% 
%         nodes_pos = reshape(softRobot.q(1:3*softRobot.n_nodes), 3, []);
%         hit_z = abs(nodes_pos(3, max_node)); 
%         collision_dists(is_colliding_now) = hit_z;
%     end
% 
%     previous_peak_forces = imc.peak_force; 
% 
%     % 4. 断裂事件捕获
%     if any(imc.is_broken)
%         newly_broken = imc.is_broken & isinf(break_times);
%         if any(newly_broken)
%             broken_indices = find(newly_broken);
%             for idx = broken_indices
%                 break_times(idx) = ctime;
%                 fprintf('>>> 冰柱 #%02d 断裂! 时刻: %.4fs | 峰值力: %6.1f N | 距根部距离: %.3f m\n', ...
%                     idx, ctime, imc.peak_force(idx), collision_dists(idx));
%             end
%         end
%     end
% 
%     % 5. 接触力与扭矩记录
%     if isfield(force_now, 'contact')
%         contact_forces_vec = reshape(force_now.contact(1:3*softRobot.n_nodes), 3, []);
%         total_contact_force = sum(vecnorm(contact_forces_vec));
% 
%         F_history.contact(timeStep) = max(total_contact_force, actual_peak_force_this_step);
% 
%        % === [最终修正版] 方案 A: 提取根部单元内部扭转力矩 ===
%         n_nodes = softRobot.n_nodes;
%         twist_offset = 3 * n_nodes; % q 向量中 theta 角的起始索引
% 
%         % 确定每根绳索的边缘数量（总数 32，每根 16）
%         num_edges_per_rod = softRobot.n_edges_rod_only / 2;
%         % 第一根绳根部索引为 1，第二根绳根部索引为 17
%         root_edge_indices = [1, num_edges_per_rod + 1];
% 
%         total_internal_torque = 0;
%         for e_idx = root_edge_indices
%             % 确保索引在 q 向量和 refTwist 范围内
%             if (twist_offset + e_idx) <= length(softRobot.q) && e_idx <= length(softRobot.refTwist)
%                 % 1. 获取当前根部边缘的离散扭转角 theta
%                 theta_current = softRobot.q(twist_offset + e_idx);
% 
%                 % 2. 计算扭转偏差 (减去马达端固定角 0 和参考扭率)
%                 delta_theta = theta_current - 0 - softRobot.refTwist(e_idx);
% 
%                 % 3. M = (GJ / L) * delta_theta
%                 % 注意：这里去掉了 GJ 的索引 (e_idx)，因为 GJ 是标量
%                 % 使用 refLen(e_idx) 代表该段杆件的长度
%                 torque_i = (softRobot.GJ / softRobot.refLen(e_idx)) * delta_theta;
% 
%                 total_internal_torque = total_internal_torque + abs(torque_i);
%             end
%         end
% 
%         % 记录马达感受到的真实传导扭矩
%         F_history.motor_torque(timeStep) = total_internal_torque;      
%     end
    
%     % 6. 记录马达转速与计算功率
%     F_history.motor_omega(timeStep) = current_omega_mag;
%     F_history.motor_power(timeStep) = F_history.motor_torque(timeStep) * current_omega_mag;
% 
%     ctime = ctime + sim_params.dt;
%     softRobot.q0 = softRobot.q; 
% 
%     if sim_params.log_data && mod(timeStep, sim_params.logStep) == 0
%         dof_with_time(2:end,timeStep) = softRobot.q;
%     end
% end

% --- [循环结束后打印最终统计结果] ---
fprintf('\n====== 仿真分析结果 (自转速度: %d rad/s) ======\n', imc.omega_spin);
broken_count = sum(imc.is_broken);
fprintf('共计击碎冰柱: %d / %d\n', broken_count, imc.num_ice);
for j = 1:imc.num_ice
    if imc.is_broken(j)
        fprintf('  - 冰柱 #%02d: 断裂时刻 %.4fs | 峰值力: %6.2f N | 碰撞点距根部: %.3f m\n', ...
            j, break_times(j), imc.peak_force(j), collision_dists(j));
    else
        if imc.peak_force(j) > 0
            fprintf('  - 冰柱 #%02d: 完好(受冲击)  | 碰撞峰值力: %6.2f N | 碰撞点距根部: %.3f m\n', ...
                j, imc.peak_force(j), collision_dists(j));
        else
            fprintf('  - 冰柱 #%02d: 完好(未碰撞)  | 受力: 0.00 N\n', j);
        end
    end
end
fprintf('除冰率: %.1f%%\n', (broken_count/imc.num_ice)*100);

% 马达性能数据统计打印
fprintf('\n====== 气动马达驱动性能校核 ======\n');
energy_consumption = trapz(time_arr, F_history.motor_power);
peak_torque = max(F_history.motor_torque);
peak_power = max(F_history.motor_power);
fprintf('总除冰作业能耗估计: %.2f J\n', energy_consumption);
fprintf('马达输出峰值扭矩: %.2f N·m\n', peak_torque);
fprintf('马达瞬态峰值功率: %.2f W\n', peak_power);
fprintf('==========================\n');

%% 3. 动画生成 (Post-Processing) 
%% ==========================================
fprintf('仿真完成，开始生成动画视频...\n');
videoFileName = 'Deicing_Spin_Fixed.mp4';
v = VideoWriter(videoFileName, 'MPEG-4'); 
v.FrameRate = 30; 
open(v);

h_fig = figure('Renderer', 'opengl', 'Color', 'w'); 
set(h_fig, 'Position', [100, 100, 1024, 768]); 

plot_limit = 0.3; 
x_lims = [-plot_limit, plot_limit];
y_lims = [-plot_limit, plot_limit];
z_lims = [-0.1, 0.4]; 

targetH = [];  
targetW = [];  

for k = sim_params.logStep : sim_params.logStep : Nsteps
    t_frame = dof_with_time(1, k);
    q_frame = dof_with_time(2:end, k);

    if norm(q_frame) == 0, continue; end

    softRobot.q = q_frame; 
    figure(h_fig); 

    if t_frame < sim_params.ramp_time
        imc.theta_accumulated = 0.5 * (sim_params.omega_target / sim_params.ramp_time) * (t_frame^2);
    else
        theta_ramp = 0.5 * sim_params.omega_target * sim_params.ramp_time;
        imc.theta_accumulated = theta_ramp + sim_params.omega_target * (t_frame - sim_params.ramp_time);
    end

    imc.is_broken = (t_frame >= break_times);
    try
        plot_MultiRod(softRobot, t_frame, sim_params, environment, imc);
    catch ME
        warning('绘图出错: %s', ME.message);
        clf; plot_MultiRod(softRobot, t_frame, sim_params, environment, imc);
    end

    xlim(x_lims); ylim(y_lims); zlim(z_lims);     
    axis equal; view(3); grid on;

    ax = gca;
    ax.XLimMode = 'manual'; ax.YLimMode = 'manual'; ax.ZLimMode = 'manual';
    ax.DataAspectRatioMode = 'manual'; ax.PlotBoxAspectRatioMode = 'manual';
    ax.CameraPositionMode = 'manual'; ax.CameraTargetMode = 'manual'; ax.CameraViewAngleMode = 'manual';

    current_force = 0;
    if k <= length(F_history.contact)
        current_force = F_history.contact(k);
    end
    title(sprintf('Time: %.3fs | Force: %.1f N', t_frame, current_force));

    frame = getframe(h_fig); 
    if isempty(targetH)
        targetH = size(frame.cdata, 1);
        targetW = size(frame.cdata, 2);
    end
    if ~isequal(size(frame.cdata, 1), targetH) || ~isequal(size(frame.cdata, 2), targetW)
        frame.cdata = imresize(frame.cdata, [targetH, targetW]);
    end
    writeVideo(v, frame);
end

close(v);
fprintf('动画已保存: %s\n', videoFileName);

%% ==========================================
%% 4. 绘制马达性能综合评估曲线
%% ==========================================
figure('Name', '气动马达除冰性能评估', 'Position', [150, 50, 900, 1000]); 

cn_font = 'Microsoft YaHei'; 

% 子图1：冰层接触力与断裂点标记
subplot(4,1,1);
plot(time_arr, F_history.contact, 'r-', 'LineWidth', 1.5); 
hold on;
valid_breaks = isfinite(break_times);
scatter(break_times(valid_breaks), imc.peak_force(valid_breaks), 60, 'k', 'p', 'filled');
title('瞬态接触力与冰柱断裂', 'FontName', cn_font); 
xlabel('时间 [s]', 'FontName', cn_font); ylabel('力 [N]', 'FontName', cn_font);
legend('接触力', '冰柱断裂', 'Location', 'best', 'FontName', cn_font);
grid on;
set(gca, 'FontName', cn_font); 

% 子图2：转速曲线 (马达启动响应)
subplot(4,1,2);
plot(time_arr, F_history.motor_omega, 'b-', 'LineWidth', 1.5);
hold on;
yline(sim_params.omega_target, 'k--', '目标转速', 'FontName', cn_font);
title('马达转速响应 (\omega)', 'FontName', cn_font); 
xlabel('时间 [s]', 'FontName', cn_font); ylabel('转速 [rad/s]', 'FontName', cn_font);
grid on;
set(gca, 'FontName', cn_font);

% 子图3：驱动扭矩曲线 (冲击载荷校核)
subplot(4,1,3);
plot(time_arr, F_history.motor_torque, 'm-', 'LineWidth', 1.5); 
hold on;
[max_torque, idx_t] = max(F_history.motor_torque);
plot(time_arr(idx_t), max_torque, 'mo', 'MarkerSize', 8, 'MarkerFaceColor', 'm');
title('马达驱动扭矩', 'FontName', cn_font); 
xlabel('时间 [s]', 'FontName', cn_font); ylabel('扭矩 [N·m]', 'FontName', cn_font);
legend('输出扭矩', sprintf('峰值: %.2f N·m', max_torque), 'Location', 'best', 'FontName', cn_font);
grid on;
set(gca, 'FontName', cn_font);

% 子图4：动态功率消耗 (供气选型依据)
subplot(4,1,4);
plot(time_arr, F_history.motor_power, 'g-', 'LineWidth', 1.5);
hold on;
[max_power, idx_p] = max(F_history.motor_power);
plot(time_arr(idx_p), max_power, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
title(sprintf('马达输出功率 (总消耗能量: %.2f J)', energy_consumption), 'FontName', cn_font); 
xlabel('时间 [s]', 'FontName', cn_font); ylabel('功率 [W]', 'FontName', cn_font);
legend('输出功率', sprintf('峰值: %.2f W', max_power), 'Location', 'best', 'FontName', cn_font);
grid on;
set(gca, 'FontName', cn_font);

sgtitle('气动马达性能与除冰动力学评估', 'FontWeight', 'bold', 'FontSize', 14, 'FontName', cn_font);