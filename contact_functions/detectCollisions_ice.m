function [in_contact_edges, closest_distance, num_coll] = detectCollisions_ice(q, candidate_edges, imc, delta, contact_len, scale)
    ice_radius = imc.ice_radius * scale;      
    rod_radius = imc.rod_radius * scale;      
    contact_threshold = ice_radius + rod_radius; 
    
    num_ice = imc.num_ice;
    R_array = imc.array_radius * scale;
    L_center = imc.array_center_dist * scale;

    if isfield(imc, 'theta_accumulated')
        theta_orb = -imc.theta_accumulated; 
    else
        theta_orb = 0; 
    end
    rot_mat = [cos(theta_orb), -sin(theta_orb); sin(theta_orb), cos(theta_orb)];

    n_edges = size(candidate_edges, 1);
    minDs_all = inf(n_edges, 1);

    % 遍历所有的候选边与所有的冰柱
    for j = 1:num_ice
        if isfield(imc, 'is_broken') && imc.is_broken(j)
            continue; % 忽略已碎冰柱
        end
        
        phi_j = 2 * pi * (j - 1) / num_ice;
        P0_j_xy = [L_center + R_array * cos(phi_j); R_array * sin(phi_j)];
        ice_pos_xy = (rot_mat * P0_j_xy)'; % 转为行向量供下文计算
        
        for i = 1:n_edges
            node1_idx = candidate_edges(i, 1);
            node2_idx = candidate_edges(i, 2);
            p1 = q(3*node1_idx-2 : 3*node1_idx-1)'; 
            p2 = q(3*node2_idx-2 : 3*node2_idx-1)';
            
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
            
            dist_surface = sqrt(dist_sq) - contact_threshold;
            minDs_all(i) = min(minDs_all(i), dist_surface);
        end
    end

    col_indices = find(minDs_all < scale*(delta + contact_len));
    in_contact_edges = candidate_edges(col_indices, :); 
    
    if ~isempty(minDs_all) && any(minDs_all ~= inf)
        closest_distance = min(minDs_all) / scale;
    else
        closest_distance = inf;
    end
    num_coll = size(in_contact_edges, 1);
end