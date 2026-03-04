% robot2DescriptionDeicing.m
% 用途：纯绳子高速旋转测试（无冰柱、无接触），但已改为双绳对称构型
% 与 length_vs_r.m 配合使用，测试不同轮毂半径下的稳定扫掠半径

%% 1. 项目路径
currentFile = mfilename('fullpath');
[currentPath, ~, ~] = fileparts(currentFile);
projectRoot = fullfile(currentPath, '..', '..');
addpath(genpath(projectRoot));

%% 2. 仿真控制参数 (已对齐 robotDescriptionDeicing.m)
sim_params.static_sim     = false;          % 必须动态
sim_params.TwoDsim        = false;          % 允许3D摆动+扭转
sim_params.use_lineSearch = true;           % 高速旋转强烈建议开启
sim_params.use_midedge    = false;          % 绳子不需要中边弯曲模型
sim_params.showFrames     = false;          % 不显示参考框架，加快速度

% 时间与精度
sim_params.dt         = 1e-4;               % 高速旋转建议小步长
sim_params.totalTime  = 5.0;                % 仿真时长
sim_params.tol        = 1e-3;               % 【修改】对齐主配置的力容差
sim_params.ftol       = 1e-3;               % 【修改】
sim_params.dtol       = 0.01;               % 【修改】位移速度容差
sim_params.maximum_iter = 20;               % 【修改】最大迭代次数

% 日志与绘图
sim_params.log_data   = true;
sim_params.logStep    = 10;                 % 【修改】对齐主配置
sim_params.plotStep   = 100;                

% 转速爬坡
sim_params.ramp_time  = 0.5;                % 【修改】对齐主配置的0.5秒加速

RPM_target            = 1000;
sim_params.omega_target = RPM_target * 2 * pi / 60;
sim_params.alpha      = sim_params.omega_target / 3;
sim_params.omega      = [0; 0; 0];          

%% 3. 几何 ── 程序生成对称的两根径向绳子
% 轮毂半径由外部传入（length_vs_r.m 中设置 geom.hub_radius）
if ~isfield(geom, 'hub_radius')
    geom.hub_radius = 0.05;                % 【修改】默认值对齐主配置的 0.05m
end

R_hub   = geom.hub_radius;
L_rope  = 0.246;                           % 【修改】统一定义绳长为 0.246m (参考真实输入)

% 固定节点数
n_nodes = 17;                              % 单根绳子的节点数 (可调)

% 【修改】生成两根对称的绳子 (类似 txt 文件中的拓扑)
s = linspace(0, L_rope, n_nodes)';
% 第一根绳子：沿 +x 轴径向向外
nodes1 = [R_hub + s, zeros(n_nodes,1), zeros(n_nodes,1)];  % [x, y, z]
% 第二根绳子：沿 -x 轴径向向外 (对称)
nodes2 = [-R_hub - s, zeros(n_nodes,1), zeros(n_nodes,1)]; 

nodes = [nodes1; nodes2];

% 边连接（建立两条独立的链段）
edges1 = [(1:n_nodes-1)', (2:n_nodes)'];
edges2 = [(n_nodes+1:2*n_nodes-1)', (n_nodes+2:2*n_nodes)'];
edges = [edges1; edges2];

% 无壳体、无面
face_nodes = [];
face_edges = [];

geom.n_rod         = 2;                    % 【修改】两根绳索
geom.n_shell       = 0;
geom.inputFileName = '';                   

%% 4. 几何与材料参数 (已完全对齐 robotDescriptionDeicing.m)

% 杆（绳子）参数
geom.rod_r0     = 0.002;                    % 绳半径 2mm
geom.Axs        = pi * geom.rod_r0^2;
geom.Ixs1       = pi * geom.rod_r0^4 / 4;
geom.Ixs2       = pi * geom.rod_r0^4 / 4;
geom.Jxs        = pi * geom.rod_r0^4 / 2;

% 壳体参数 - 全部清零
geom.shell_h    = 0;                        
geom.shell_r    = 0;
geom.Axs_shell  = 0;
geom.Ixs_shell1 = 0;
geom.Ixs_shell2 = 0;
geom.Jxs_shell  = 0;

% 材料
material.density          = 1100;           % kg/m³
material.youngs_rod       = 1e9;            % 【修改】对齐主配置的 1e9 Pa
material.poisson_rod      = 0.35;
material.youngs_shell     = 0;
material.poisson_shell    = 0;
material.contact_stiffness = 50000;         % 【修改】对齐主配置 (尽管未开启接触计算)
material.mu               = 0.3;            % 【修改】对齐主配置
material.velTol           = 1e-4;

% 数值阻尼
env.eta = 0.1;                              % 【修改】粘性阻尼对齐主配置的 0.1

%% 5. 外部力（只保留旋转相关，无接触）
env.ext_force_list = ["gravity", "viscous", "aerodynamic", "centrifugal", "coriolis"];

env.g = [0; 0; -9.81];
env.air_density = 1.225;
env.Cd          = 1.0;                      % 【修改】阻力系数对齐 1.0
env.rho         = 1.225;

% 关闭所有接触 (纯绳子旋转测试不需要)
env.selfContact   = false;
env.floorContact  = false;

% 清空冰柱相关参数（防止底层代码读取报错）
env.contact_params.num_ice         = 0;
env.contact_params.is_broken       = [];
env.contact_params.peak_force      = [];
env.contact_params.ice_radius      = 0;
env.contact_params.array_radius    = 0;

%% 6. 边界条件
% 【修改】同时固定两根绳子的根部节点 (第1个节点 和 第n_nodes+1个节点)
fixed_node_indices = [1, n_nodes + 1];                     

% 记录节点选为第一根绳子的末端
input_log_node = n_nodes;            

%% 7. 绘图范围
plot_limit = L_rope * 1.1 + 0.05;
sim_params.plot_x = [-plot_limit, plot_limit];
sim_params.plot_y = [-plot_limit, plot_limit];
sim_params.plot_z = [-0.1, plot_limit * 0.7];

%% 8. 提示
fprintf('robot2DescriptionDeicing 加载完成 (双绳对称版)\n');
fprintf('  轮毂半径   : %.4f m\n', R_hub);
fprintf('  单绳长度   : %.4f m\n', L_rope);
fprintf('  节点数/绳  : %d \n', n_nodes);
fprintf('  目标转速   : %d rpm (%.1f rad/s)\n\n', RPM_target, sim_params.omega_target);