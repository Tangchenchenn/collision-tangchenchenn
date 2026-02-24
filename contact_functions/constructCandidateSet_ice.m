function [candidate_edges, closest_dist] = constructCandidateSet_ice(q, rod_edges, candidate_limit, scale)
% constructCandidateSet_ice: 构建除冰绳与固定冰柱的接触候选集
% q: 当前所有自由度的坐标向量
% rod_edges: 除冰绳的边连接信息 (n_edges x 2)
% candidate_limit: 候选判定阈值 (已包含 scale)
% scale: 长度缩放因子

% --- 几何参数 (基于 robotDescriptionDeicing.m) ---
ice_center_dist = 0.15 * scale;  % 冰柱中心距离旋转中心 0.15m
ice_radius = 0.004 * scale;      % 冰柱半径 4mm
rod_radius = 0.002 * scale;      % 绳索半径 2mm
contact_threshold = ice_radius + rod_radius; % 物理接触界限 (0.006m)

n_edges = size(rod_edges, 1);
minDs = zeros(n_edges, 1);

% 假设冰柱轴线平行于 Z 轴，中心位于随动系的 X 轴正方向 [0.15, 0, z]
ice_pos_xy = [ice_center_dist, 0];

% 遍历除冰绳的每一条边
for i = 1:n_edges
    % 获取边的两个节点索引及坐标
    node1_idx = rod_edges(i, 1);
    node2_idx = rod_edges(i, 2);
    
    % 提取节点在 XY 平面的坐标 (假设 mapNodetoDOF 逻辑)
    p1 = q(3*node1_idx-2 : 3*node1_idx-1)'; % [x1, y1]
    p2 = q(3*node2_idx-2 : 3*node2_idx-1)'; % [x2, y2]
    
    % 计算线段 (p1, p2) 到点 ice_pos_xy 的最短距离
    % 这是标准的点到线段距离算法
    v = p2 - p1;
    w = ice_pos_xy - p1;
    c1 = dot(w, v);
    if c1 <= 0
        dist_sq = sum((ice_pos_xy - p1).^2);
    else
        c2 = dot(v, v);
        if c2 <= c1
            dist_sq = sum((ice_pos_xy - p2).^2);
        else
            b = c1 / c2;
            pb = p1 + b * v;
            dist_sq = sum((ice_pos_xy - pb).^2);
        end
    end
    
    % 减去半径，得到表面到表面的距离
    minDs(i) = sqrt(dist_sq) - contact_threshold;
end

% 计算全局最小距离（用于返回信息）
closest_dist = min(minDs) / scale;

% 筛选在 candidate_limit 范围内的边
% 注意：此处 minDs 是物体表面的实际距离，需与 candidate_limit 匹配
col_indices = find(minDs < candidate_limit);
candidate_edges = rod_edges(col_indices, :);

end