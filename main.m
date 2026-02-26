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
break_times = inf(1, env.contact_params.num_ice); % 记录每根冰柱的断裂时间
for timeStep = 1:Nsteps
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
    % ================= [多冰柱阵列：断裂事件捕获] =================
    if any(imc.is_broken)
        % 找出这一帧新断裂的冰柱
        newly_broken = imc.is_broken & isinf(break_times);
        if any(newly_broken)
            broken_indices = find(newly_broken);
            for idx = broken_indices
                break_times(idx) = ctime;
                fprintf('>>> 冰柱 #%d 断裂! 时刻: %.4fs | 峰值力: %.1fN\n', ...
                    idx, ctime, imc.peak_force(idx));
            end
        end
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
% --- [循环结束后打印结果] ---
fprintf('\n====== 仿真分析结果 (自转速度: %d rad/s) ======\n', imc.omega_spin);
broken_count = sum(imc.is_broken);
fprintf('共计击碎冰柱: %d / %d\n', broken_count, imc.num_ice);
for j = 1:imc.num_ice
    if imc.is_broken(j)
        fprintf('  - 冰柱 #%d: 断裂时刻 %.4fs | 峰值力: %.2f N\n', j, break_times(j), imc.peak_force(j));
    else
        fprintf('  - 冰柱 #%d: 完好\n', j);
    end
end
fprintf('除冰率: %.1f%%\n', (broken_count/imc.num_ice)*100);
fprintf('==========================\n');
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
    % 放在调用 plot_MultiRod 之前
    if t_frame < sim_params.ramp_time
        imc.theta_accumulated = 0.5 * (sim_params.omega_target / sim_params.ramp_time) * (t_frame^2);
    else
        theta_ramp = 0.5 * sim_params.omega_target * sim_params.ramp_time;
        imc.theta_accumulated = theta_ramp + sim_params.omega_target * (t_frame - sim_params.ramp_time);
    end
    % === 【核心修复】：动态重构当前帧的冰柱断裂状态 ===
    imc.is_broken = (t_frame >= break_times);
    try
        plot_MultiRod(softRobot, t_frame, sim_params, environment, imc);
    catch ME
        warning('绘图出错: %s', ME.message);
        clf; plot_MultiRod(softRobot, t_frame, sim_params, environment, imc);
    end
    
% 2. 【关键修改】彻底锁死坐标轴和相机视角 (杜绝一切画面缩放和抖动)
    % 必须在 plot_MultiRod 之后执行，覆盖它的默认设置
    xlim(x_lims);     
    ylim(y_lims);     
    zlim(z_lims);     
    axis equal;       
    view(3);          
    grid on;

    % 获取当前坐标轴句柄
    ax = gca;
    % 锁死数据范围
    ax.XLimMode = 'manual';
    ax.YLimMode = 'manual';
    ax.ZLimMode = 'manual';
    % 锁死长宽比例
    ax.DataAspectRatioMode = 'manual';
    ax.PlotBoxAspectRatioMode = 'manual';
    % 锁死相机的距离、目标点和视野角度 (解决画面忽大忽小的根本方法)
    ax.CameraPositionMode = 'manual';
    ax.CameraTargetMode = 'manual';
    ax.CameraViewAngleMode = 'manual';

    % 3. 更新标题
    current_force = 0;
    if k <= length(F_history.contact)
        current_force = F_history.contact(k);
    end
    title(sprintf('Time: %.3fs | Force: %.1f N', t_frame, current_force));
    
        % 4. 写入视频
        frame = getframe(gcf);
        
        % --- [新增：强制统一分辨率以避开报错] ---
        % 将每一帧强制缩放到 v.Height x v.Width (即第一帧或默认的大小)
        if size(frame.cdata, 1) ~= v.Height || size(frame.cdata, 2) ~= v.Width
            frame.cdata = imresize(frame.cdata, [v.Height, v.Width]);
        end
        % ---------------------------------------
        
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