% 添加项目根目录及其所有子文件夹到路径
currentFile = mfilename('fullpath');
[currentPath, ~, ~] = fileparts(currentFile);
projectRoot = fullfile(currentPath, '..', '..'); % 根据你的文件夹深度调整
addpath(genpath(projectRoot));

sim_params.static_sim = false; % 必须设为动态仿真以模拟旋转动力学
sim_params.TwoDsim = false;    % 建议设为false，以捕捉三维摆动和扭转耦合
sim_params.use_lineSearch = true; % 高速旋转或碰撞时，建议开启线搜索以增强稳定性
% --- robotDescriptionDeicing.m 中添加 ---
sim_params.use_midedge = false; % 用于控制壳体弯曲模型，除冰绳仿真设为 false 即可
% 在 robotDescriptionDeicing.m 中添加
sim_params.showFrames = false; % 设为 false 表示不显示 Bishop 框架，仅观察绳索构型

% 气动马达驱动参数 (对应你之前的运动方程)
sim_params.ramp_time = 0.5;% 设定 0.5s 从 0 加速到目标转速
RPM_target = 1000; % 定义马达转速为 2000 r/min
sim_params.omega_target = (RPM_target * 2 * pi) / 60;
sim_params.alpha = sim_params.omega_target / 3;
sim_params.hub_radius = 0.05;      % 轮毂半径

% % 将标量转速转化为向量（假设绕 Z 轴旋转）
% sim_params.omega = [0; 0; sim_params.omega_target];
% 初始化当前角速度和角度
sim_params.omega = [0; 0; 0];       % 初始角速度向量为 0
imc.theta_accumulated = 0;          % [新增] 用于精确记录冰柱转过的角度
% 时间步长 (高速旋转下建议调小步长)
sim_params.dt = 1e-4; 
sim_params.totalTime = 5.0; % 仿真时长

% --- robotDescriptionDeicing.m 中添加 ---
sim_params.tol = 1e-3;           % 力的平衡容差 (Newton 迭代停止标准)
sim_params.ftol = 1e-3;         % 相对力容差
sim_params.dtol = 0.01;          % 节点位移速度容差
sim_params.maximum_iter = 20;    % Newton 迭代的最大次数，防止死循环

% --- robotDescriptionDeicing.m 中添加 ---
sim_params.log_data = true;   % 是否记录所有自由度的数据
sim_params.logStep = 10;       % 记录频率，1 表示每步都记录
sim_params.plotStep =100;     % 绘图频率，建议设大一点（如 10）以加快运行速度

%% 引用刚才生成的几何文件
geom.inputFileName = 'deicing_cable_symmetric.txt'; 
[nodes, edges, face_nodes] = inputProcessorNew(geom.inputFileName);

%% 除冰绳物理参数
geom.rod_r0 = 0.002; % 绳索截面半径 (m)

geom.shell_h = 0;    % 即使没有壳体，也需初始化为0
% 计算截面属性
r = geom.rod_r0;
geom.Axs = pi * r^2;
geom.Ixs1 = (pi * r^4) / 4;
geom.Ixs2 = (pi * r^4) / 4;
geom.Jxs = (pi * r^4) / 2;

% 材料参数 (需根据实际材料如钢索或复合绳调整)
material.youngs_shell = 0;  % 壳体杨氏模量，设为0即可
material.poisson_shell = 0; % 壳体泊松比，设为0即可
material.density = 1100;      % 密度 (kg/m^3)
material.youngs_rod = 1e9;   % 杨氏模量 E (Pa)
material.poisson_rod = 0.35;   % 泊松比 (计算剪切模量 G 用)
material.contact_stiffness = 50000; % 接触刚度 (建议值 1000 - 10000)
material.mu = 0.3;                 % 摩擦系数 (防止未来开启 selfFriction 报错)
material.velTol = 1e-4;            % 速度容差 (防止未来开启 selfFriction 报错)
% --- robotDescriptionDeicing.m 中添加 ---
env.eta = 0.1; % 粘性阻尼系数 (推荐值 0.01-0.1 之间，用于数值稳定性)

%% 外部力设置
% 必须包含 "aerodynamic" 以模拟高速下的空气阻力，"gravity" 模拟重力
% 确保环境参数包含这些力
env.ext_force_list = ["gravity", "viscous", "aerodynamic", "centrifugal", "coriolis", "selfContact"];
% viscous（粘性力）有助于数值稳定，aerodynamic（空气动力学）模拟除冰绳高速旋转时的风阻
env.g = [0; 0; -9.81]; % 重力加速度
env.air_density = 1.225; % 空气密度 (kg/m^3)
env.Cd = 1.0; % 阻力系数
env.rho = 1.225; % 空气密度 (kg/m^3)
% ==============================
% 冰柱几何与运动参数（圆周阵列配置）
env.contact_params.ice_radius = 0.015;       % 单根冰柱半径 (16mm)
env.contact_params.num_ice = 10;             % [修改] 冰柱数量 m=10
env.contact_params.array_radius = 0.1;       % [新增] 圆周的半径 R=0.1
% [修改] 圆周的圆心距除冰绳旋转中心距离 = 0.1 + 0.15 = 0.25
env.contact_params.array_center_dist = env.contact_params.array_radius + 0.15; 

env.contact_params.omega_mag = sim_params.omega_target; % 公转角速度
env.contact_params.omega_spin = 1;         % [新增] 冰柱自转角速度(rad/s) -> 后续改这个值作对比分析

env.contact_params.compute_friction = true; 
env.contact_params.active_time = 0.8;
% ==============================
env.contact_params.sigma_t = 1.5e6;  
env.contact_params.z_root = 0.07;     
env.contact_params.rod_radius = geom.rod_r0;

% [修改] 将断裂状态和峰值力记录改为数组，每根冰柱独立记录
env.contact_params.is_broken = false(1, env.contact_params.num_ice); 
env.contact_params.peak_force = zeros(1, env.contact_params.num_ice); 
% =============================================
%% 边界条件
fixed_node_indices = [1, 14]; % 同时固定两根绳子的根部节点
input_log_node = size(nodes, 1);


%% Plot dimensions (在 robotDescriptionDeicing.m 中修改)
% 5. 坐标轴范围补全 (防止报错：无法识别字段 plot_x)
% 如果已经在 robotDescription 里改过，这里可以不改，或者在此处强制指定：
sim_params.plot_x = [-0.3, 0.3]; 
sim_params.plot_y = [-0.3, 0.3];
sim_params.plot_z = [-0.1, 0.5]; % 考虑到 3D 甩起的高度
