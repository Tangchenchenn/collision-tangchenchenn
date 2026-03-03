% robot2DescriptionDeicing.m
% 用途：纯绳子高速旋转测试（无冰柱、无接触）
% 与 length_vs_r.m 配合使用，测试不同轮毂半径下的稳定扫掠半径

%% 1. 项目路径
currentFile = mfilename('fullpath');
[currentPath, ~, ~] = fileparts(currentFile);
projectRoot = fullfile(currentPath, '..', '..');
addpath(genpath(projectRoot));

%% 2. 仿真控制参数
sim_params.static_sim     = false;          % 必须动态
sim_params.TwoDsim        = false;          % 允许3D摆动+扭转
sim_params.use_lineSearch = true;           % 高速旋转强烈建议开启
sim_params.use_midedge    = false;          % 绳子不需要中边弯曲模型
sim_params.showFrames     = false;          % 不显示参考框架，加快速度

% 时间与精度
sim_params.dt         = 1e-4;               % 高速旋转建议小步长
sim_params.totalTime  = 5.0;                % 5秒通常足够观察稳定
sim_params.tol        = 5e-4;               % 稍放宽一点提高收敛率
sim_params.ftol       = 5e-3;
sim_params.dtol       = 0.005;
sim_params.maximum_iter = 25;

% 日志与绘图
sim_params.log_data   = true;
sim_params.logStep    = 20;                 % 每20步记录一次，节省内存
sim_params.plotStep   = 100;                % 绘图频率

% 转速爬坡
sim_params.ramp_time  = 0.3;                % 0.3秒线性加速到目标转速

RPM_target            = 1000;
sim_params.omega_target = RPM_target * 2 * pi / 60;
sim_params.omega      = [0; 0; 0];          % 初始为0，由主程序逐步爬坡

%% 3. 几何 ── 程序生成简单径向单根绳子
% 轮毂半径由外部传入（length_vs_r.m 中设置 geom.hub_radius）
if ~isfield(geom, 'hub_radius')
    geom.hub_radius = 0.025;               % 默认值，仅用于单独测试
end

R_hub   = geom.hub_radius;
L_rope  = R_hub * 7;                       % 核心比例 7:1

% 固定节点数，减少不同案例间的离散化差异
n_nodes = 81;                              % 固定为81个节点

% 初始位形：沿x轴径向直线（离心力会甩开）
s = linspace(0, L_rope, n_nodes)';
nodes = [R_hub + s, zeros(n_nodes,1), zeros(n_nodes,1)];  % [x, y, z]

% 边连接（单根连续杆）
edges = [(1:n_nodes-1)', (2:n_nodes)'];

% 无壳体、无面
face_nodes = [];
face_edges = [];

geom.n_rod         = 1;
geom.rod_start     = 1;
geom.n_shell       = 0;
geom.inputFileName = '';                   % 不再使用输入文件

%% 4. 几何与材料参数（纯绳子版，强制补全壳体相关字段）

% 杆（绳子）参数
geom.rod_r0     = 0.002;                    % 绳半径 2mm
geom.Axs        = pi * geom.rod_r0^2;
geom.Ixs1       = pi * geom.rod_r0^4 / 4;
geom.Ixs2       = pi * geom.rod_r0^4 / 4;
geom.Jxs        = pi * geom.rod_r0^4 / 2;

% 壳体参数 - 全部清零（防止 MultiRod 构造函数报错）
geom.shell_h    = 0;                        % 厚度，必须补上！
geom.shell_r    = 0;
geom.Axs_shell  = 0;
geom.Ixs_shell1 = 0;
geom.Ixs_shell2 = 0;
geom.Jxs_shell  = 0;

% 材料
material.density          = 1100;           % kg/m³
material.youngs_rod       = 8e8;            % Pa （可调：更软→更易甩开）
material.poisson_rod      = 0.35;
material.youngs_shell     = 0;
material.poisson_shell    = 0;
material.contact_stiffness = 0;             % 无接触
material.mu               = 0;
material.velTol           = 1e-4;

% 数值阻尼（强烈建议保留少量）
env.eta = 0.05;                             % 粘性阻尼系数，0.01~0.1 常用范围

%% 5. 外部力（只保留旋转相关）
env.ext_force_list = ["gravity", "viscous", "aerodynamic", "centrifugal", "coriolis"];

env.g = [0; 0; -9.81];
env.air_density = 1.225;
env.Cd          = 1.2;                      % 圆柱阻力系数典型值
env.rho         = 1.225;

% 关闭所有接触
env.selfContact   = false;
env.floorContact  = false;

% 清空冰柱相关参数（防止报错）
env.contact_params.num_ice         = 0;
env.contact_params.is_broken       = [];
env.contact_params.peak_force      = [];
env.contact_params.ice_radius      = 0;
env.contact_params.array_radius    = 0;

%% 6. 边界条件
fixed_node_indices = 1;                     % 只固定根部节点

input_log_node = size(nodes, 1);            % 末端节点，用于记录

%% 7. 绘图范围（根据当前尺度自动调整）
plot_limit = L_rope * 1.1 + 0.05;
sim_params.plot_x = [-plot_limit, plot_limit];
sim_params.plot_y = [-plot_limit, plot_limit];
sim_params.plot_z = [-0.15, plot_limit * 0.7];

%% 8. 提示
fprintf('robot2DescriptionDeicing 加载完成\n');
fprintf('  轮毂半径   : %.4f m\n', R_hub);
fprintf('  绳长       : %.4f m (比例 %.1f:1)\n', L_rope, L_rope/R_hub);
fprintf('  节点数     : %d (固定)\n', n_nodes);
fprintf('  目标转速   : %d rpm (%.1f rad/s)\n\n', RPM_target, sim_params.omega_target);