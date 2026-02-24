clc; clear all; close all;

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

%% 初始化系统
softRobot = computeSpaceParallel(softRobot);
theta = softRobot.q0(3*softRobot.n_nodes+1 : 3*softRobot.n_nodes+softRobot.n_edges_dof);
[softRobot.m1, softRobot.m2] = computeMaterialDirectors(softRobot.a1, softRobot.a2, theta);
bend_twist_springs = setkappa(softRobot, bend_twist_springs);
softRobot.undef_refTwist = computeRefTwist_bend_twist_spring(bend_twist_springs, softRobot.a1, softRobot.tangent, zeros(n_bend_twist,1));
softRobot.refTwist = computeRefTwist_bend_twist_spring(bend_twist_springs, softRobot.a1, softRobot.tangent, softRobot.undef_refTwist);

% 边界条件
% n_nodes_per_rod = 13; 
% fixed_node_indices = [1, n_nodes_per_rod + 1]; 
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
                   'motor_torque', zeros(Nsteps,1)); % <--- 新增这一个字段

fprintf('开始物理仿真 (共 %d 步)...\n', Nsteps);

% --- [循环前初始化] ---
breaking_event_data = struct('time', NaN, 'force', 0, 'z_pos', NaN);
has_captured_break = false;
for timeStep = 1:Nsteps
    % % --- 1. 马达驱动修正 ---
    % % 【修正】在随动坐标系下，根节点不需要旋转！它相对于坐标系是静止的。
    % % 我们只更新时间，不更新 softRobot.q 的固定点坐标。
    % 
    % % --- 2. 调用时间步进器 ---
    % if(sim_params.use_midedge), tau_0 = updatePreComp_without_sign(softRobot.q, softRobot); else, tau_0 = []; end
    % 
    % % 传入 ctime 用于计算冰柱位置 (IceContact 需要)
    % % [softRobot, stretch_springs, bend_twist_springs, hinge_springs, force_now] = ...
    % %     timeStepper(softRobot, stretch_springs, bend_twist_springs, hinge_springs, ...
    % %     triangle_springs, tau_0, environment, imc, sim_params, ctime);
    % % [softRobot, stretch_springs, bend_twist_springs, hinge_springs, force_now, imc] = ...
    % %     timeStepper(softRobot, stretch_springs, bend_twist_springs, hinge_springs, ...
    % %     triangle_springs, tau_0, environment, imc, sim_params, ctime);
    % [softRobot, stretch_springs, bend_twist_springs, hinge_springs, force_now, imc] = ...
    %     timeStepper(softRobot, stretch_springs, bend_twist_springs, hinge_springs, ...
    %     triangle_springs, tau_0, environment, imc, sim_params, ctime);
    % 1. 计算当前瞬时角速度 (线性爬升)
    if ctime < sim_params.ramp_time
        ratio = ctime / sim_params.ramp_time;
        current_omega_mag = ratio * sim_params.omega_target;
    else
        current_omega_mag = sim_params.omega_target;
    end
    
    % 2. 统一更新所有结构体中的转速参数
    sim_params.omega = [0; 0; current_omega_mag]; % 供离心力、空气阻力使用
    imc.omega_mag = current_omega_mag;            % 供 IceContact 速度项使用
    
    % 3. [关键] 精确积分计算当前冰柱的角度位置
    % 角度 = 上一步角度 + 当前转速 * 时间步长
    imc.theta_accumulated = imc.theta_accumulated + current_omega_mag * sim_params.dt;

    % 4. 调用时间步进器
    if(sim_params.use_midedge), tau_0 = updatePreComp_without_sign(softRobot.q, softRobot); else, tau_0 = []; end
    
    [softRobot, stretch_springs, bend_twist_springs, hinge_springs, force_now, imc] = ...
        timeStepper(softRobot, stretch_springs, bend_twist_springs, hinge_springs, ...
        triangle_springs, tau_0, environment, imc, sim_params, ctime);
    if imc.is_broken && ~has_captured_break
        breaking_event_data.time = ctime;
        breaking_event_data.force = imc.peak_force;
        % 记录发生断裂的节点位置（可选）
        breaking_event_data.z_pos = softRobot.q(3*softRobot.n_nodes); 
        
        has_captured_break = true; % 确保只记录第一次断裂
    end
    % --- 3. 数据记录修正 ---
    % 【修正】记录所有节点的接触力总和，而不仅仅是末端节点
    if isfield(force_now, 'contact')
        % 将力向量 reshape 为 [3 x n_nodes]，计算每列的模，再求和
        contact_forces_vec = reshape(force_now.contact(1:3*softRobot.n_nodes), 3, []);
        total_contact_force = sum(vecnorm(contact_forces_vec));
        F_history.contact(timeStep) = total_contact_force;

        % === [新增] 计算马达驱动扭矩 (Motor Torque about Z-axis) ===
        % 阻碍旋转的力包括：空气阻力 (drag)、科氏力 (coriolis)、接触摩擦力 (contact)
        F_resistive = force_now.drag + force_now.coriolis + force_now.contact; 
        F_res_vec = reshape(F_resistive(1:3*softRobot.n_nodes), 3, []);
        
        % 提取当前所有节点的坐标 (x, y)
        nodes_pos = reshape(softRobot.q(1:3*softRobot.n_nodes), 3, []);
        X = nodes_pos(1, :);
        Y = nodes_pos(2, :);
        
        % 提取阻力在 X 和 Y 方向的分量
        Fx = F_res_vec(1, :);
        Fy = F_res_vec(2, :);
        
        % 扭矩 = r x F，关于 Z 轴的扭矩 Mz = X * Fy - Y * Fx
        % 马达需要提供的扭矩与阻力扭矩大小相等、方向相反
        node_torques_z = X .* Fy - Y .* Fx;
        F_history.motor_torque(timeStep) = abs(sum(node_torques_z)); 
        % ==========================================================
    end
    
    ctime = ctime + sim_params.dt;
    softRobot.q0 = softRobot.q; 
    
    if sim_params.log_data && mod(timeStep, sim_params.logStep) == 0
        dof_with_time(2:end,timeStep) = softRobot.q;
        fprintf('进度: %.1f%% | 接触力总和: %.2e N\n', (timeStep/Nsteps)*100, F_history.contact(timeStep));
    end
end
% --- [循环结束后打印结果] ---
if has_captured_break
    fprintf('\n====== 仿真分析结果 ======\n');
    fprintf('断裂时刻: %.4f s\n', breaking_event_data.time);
    fprintf('断裂峰值力: %.2f N\n', breaking_event_data.force);
    fprintf('==========================\n');
end
%% ==========================================
% %% 3. 动画生成 (Post-Processing)
% %% ==========================================
% fprintf('仿真完成，开始生成动画视频...\n');
% videoFileName = 'Deicing_Spin_Fixed.mp4';
% v = VideoWriter(videoFileName, 'MPEG-4'); v.FrameRate = 30; open(v);
% 
% h_fig = figure('Renderer', 'opengl', 'Color', 'w'); 
% set(h_fig, 'Position', [100, 100, 800, 600]); % 设置窗口大小
% 
% for k = sim_params.logStep : sim_params.logStep : Nsteps
%     t_frame = dof_with_time(1, k);
%     q_frame = dof_with_time(2:end, k);
% 
%     % 跳过全零帧（如果有）
%     if norm(q_frame) == 0, continue; end
% 
%     softRobot.q = q_frame; 
% 
%     % ----------------------------------------------------
%     % 调用绘图函数
%     % plot_MultiRod 内部现在会自动处理：
%     % 1. 绘制不动的绳子 (q 本身就是随动坐标)
%     % 2. 绘制绕圈的冰柱 (根据 t_frame 计算角度)
%     % ----------------------------------------------------
%     figure(h_fig);
%     try
%         plot_MultiRod(softRobot, t_frame, sim_params, environment, imc);
%     catch ME
%         warning('绘图出错: %s', ME.message);
%         clf; plot_MultiRod(softRobot, t_frame, sim_params, environment, imc);
%     end
% 
%     % 添加额外的标题信息 (力的大小)
%     current_force = 0;
%     if k <= length(F_history.contact)
%         current_force = F_history.contact(k);
%     end
%     title(sprintf('Time: %.3fs | Contact Force: %.1f N', t_frame, current_force));
% 
%     % 写入视频
%     frame = getframe(gcf);
%     writeVideo(v, frame);
% end
% 
% close(v);
% fprintf('动画已保存: %s\n', videoFileName);
% 
% % 绘制总接触力曲线
% figure; 
% plot(time_arr, F_history.contact, 'r-', 'LineWidth', 1.5); 
% title('Total Contact Force over Time'); 
% xlabel('Time [s]'); ylabel('Force [N]');
% grid on;
%% ==========================================
%% 3. 动画生成 (Post-Processing) - 修正版
%% ==========================================
fprintf('仿真完成，开始生成动画视频...\n');
videoFileName = 'Deicing_Spin_Fixed.mp4';
v = VideoWriter(videoFileName, 'MPEG-4'); 
v.FrameRate = 30; 
open(v);

h_fig = figure('Renderer', 'opengl', 'Color', 'w'); 
set(h_fig, 'Position', [100, 100, 1024, 768]); % 设置固定的窗口大小

% 定义固定的坐标轴范围 (根据你的场景大小调整)
% 考虑到冰柱距离中心 0.15m，绳子甩起来可能到 0.2m+
plot_limit = 0.3; 
x_lims = [-plot_limit, plot_limit];
y_lims = [-plot_limit, plot_limit];
z_lims = [-0.1, 0.4]; % Z轴通常需要根据绳子高度调整

for k = sim_params.logStep : sim_params.logStep : Nsteps
    t_frame = dof_with_time(1, k);
    q_frame = dof_with_time(2:end, k);
    
    % 跳过全零帧
    if norm(q_frame) == 0, continue; end
    
    softRobot.q = q_frame; 
    
    figure(h_fig); % 激活当前窗口
    
    % 1. 调用绘图
    % 注意：plot_MultiRod 内部有 clf，会清除上一帧，这很好
    try
        plot_MultiRod(softRobot, t_frame, sim_params, environment, imc);
    catch ME
        warning('绘图出错: %s', ME.message);
        clf; plot_MultiRod(softRobot, t_frame, sim_params, environment, imc);
    end
    
    % 2. 【关键修改】强制锁定坐标轴和视角
    % 必须在 plot_MultiRod 之后执行，覆盖它的默认设置
    axis equal;       % 保证比例一致
    xlim(x_lims);     % 强制 X 范围
    ylim(y_lims);     % 强制 Y 范围
    zlim(z_lims);     % 强制 Z 范围
    grid on;
    view(3);          % 锁定三维视角
    
    % 3. 更新标题
    current_force = 0;
    if k <= length(F_history.contact)
        current_force = F_history.contact(k);
    end
    title(sprintf('Time: %.3fs | Force: %.1f N', t_frame, current_force));
    
    % 4. 写入视频
    frame = getframe(gcf);
    writeVideo(v, frame);
end

close(v);
fprintf('动画已保存: %s\n', videoFileName);

% 绘制力曲线 (建议单独画一个清晰的图)
% 绘制总接触力曲线
figure; 
subplot(2,1,1);
plot(time_arr, F_history.contact, 'r-', 'LineWidth', 1.5); 
title('Total Contact Force (Sampled)'); 
xlabel('Time [s]'); ylabel('Force [N]');
grid on;

% 绘制马达扭矩曲线
subplot(2,1,2);
plot(time_arr, F_history.motor_torque, 'b-', 'LineWidth', 1.5); 
title('Motor Driving Torque (Z-axis)'); 
xlabel('Time [s]'); ylabel('Torque [N·m]');
grid on;