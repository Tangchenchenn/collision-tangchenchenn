% create_deicing_cable_input_symmetric.m
n_nodes_per_rod = 13; % 每根绳子的节点数
R_hub = 0.016;        % 轮毂半径 [cite: 295]
L_total = 0.246;      % 绳索总长 [cite: 295]
phi_0A = 0;           % 第一根绳子的初始安装角
l_bar = L_total / (n_nodes_per_rod - 1); 

%% 1. 生成第一根绳子的节点 (Rod 1)
nodes1 = zeros(n_nodes_per_rod, 3);
nodes1(1, :) = [R_hub * cos(phi_0A), R_hub * sin(phi_0A), 0]; % 起点 x0 [cite: 41]

for k = 2:n_nodes_per_rod
    sum_val = [0, 0, 0];
    for i = 0:(k-2)
        % 按照文档公式(3.40)计算非线性构型分量 [cite: 40]
        term_x = cos(pi/4 - (i*pi)/(4*(n_nodes_per_rod-1))) * cos((i*pi)/(n_nodes_per_rod-1));
        term_y = cos(pi/4 - (i*pi)/(4*(n_nodes_per_rod-1))) * sin((i*pi)/(n_nodes_per_rod-1));
        sum_val = sum_val + l_bar * [term_x, term_y, 0];
    end
    nodes1(k, :) = nodes1(1, :) + sum_val;
end

%% 2. 生成第二根绳子的节点 (Rod 2 - 对称镜像)
% 因为对称安装在轮毂两边，坐标直接取负即可实现 180 度对称
nodes2 = -nodes1; 

% 合并所有节点
total_nodes = [nodes1; nodes2];

%% 3. 生成边连接关系
edges1 = zeros(n_nodes_per_rod-1, 2);
for i = 1:n_nodes_per_rod-1
    edges1(i, :) = [i, i+1];
end

edges2 = zeros(n_nodes_per_rod-1, 2);
for i = 1:n_nodes_per_rod-1
    % 第二根绳子的节点索引从 n_nodes_per_rod + 1 开始
    edges2(i, :) = [i + n_nodes_per_rod, i+1 + n_nodes_per_rod];
end

% 合并所有边
total_edges = [edges1; edges2];

%% 4. 保存为文本文件
filename = 'deicing_cable_symmetric.txt';
fid = fopen(filename, 'w');
if fid ~= -1
    fprintf(fid, '*Nodes\n');
    fclose(fid);
end
writematrix(total_nodes, filename, 'WriteMode', 'append')
fid = fopen(filename, 'at');
if fid ~= -1
    fprintf(fid, '*Edges\n');
    fclose(fid);
end
writematrix(total_edges, filename, 'WriteMode', 'append')

% %% 5. 预览与检查
% figure();
% plot(total_nodes(1:n_nodes_per_rod, 1), total_nodes(1:n_nodes_per_rod, 2), 'b-o', 'DisplayName', 'Rod 1'); hold on;
% plot(total_nodes(n_nodes_per_rod+1:end, 1), total_nodes(n_nodes_per_rod+1:end, 2), 'r-o', 'DisplayName', 'Rod 2');
% axis equal; grid on;
% legend; title('双侧对称除冰绳初始构型预览');