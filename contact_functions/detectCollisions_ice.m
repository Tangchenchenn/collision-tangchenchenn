function [in_contact_edges, closest_distance, num_coll] = detectCollisions_ice(q, candidate_edges, imc, delta, contact_len, scale)
% 从 imc 结构体动态读取参数，不再写死
ice_radius = imc.ice_radius * scale;      % 动态获取半径
rod_radius = imc.rod_radius * scale;      % 动态获取绳索半径
ice_center_dist = imc.ice_center_dist * scale; % 动态获取中心距离

contact_threshold = ice_radius + rod_radius; 
ice_pos_xy = [ice_center_dist, 0];

n_edges = size(candidate_edges, 1);
minDs = zeros(n_edges, 1);

% 1. 计算所有候选边到冰柱表面的最小距离
for i = 1:n_edges
    node1_idx = candidate_edges(i, 1);
    node2_idx = candidate_edges(i, 2);
    
    % 提取节点在 XY 平面的坐标 (假设 mapNodetoDOF 逻辑)
    p1 = q(3*node1_idx-2 : 3*node1_idx-1)'; 
    p2 = q(3*node2_idx-2 : 3*node2_idx-1)';
    
    % 计算线段 (p1, p2) 到点 ice_pos_xy 的最短距离 (点到线段距离算法)
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
    
    % 表面到表面的距离 (物体表面间距)
    minDs(i) = sqrt(dist_sq) - contact_threshold;
end

% 2. 确定在碰撞激活区内的索引 [参照原文判定逻辑: < scale*(delta + contact_len)]
% 这里 minDs 已经是相对于表面的距离，故直接与激活阈值比较
col_indices = find(minDs < scale*(delta + contact_len));

% 3. 封装返回结果
in_contact_edges = candidate_edges(col_indices, :); 
if ~isempty(minDs)
    closest_distance = min(minDs) / scale;
else
    closest_distance = inf;
end
num_coll = size(in_contact_edges, 1);

end