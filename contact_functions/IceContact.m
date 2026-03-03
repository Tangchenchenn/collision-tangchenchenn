function [Fc, Jc, Ffr, Jfr, imc] = IceContact(imc, q, q0, rod_edges, iter, dt, current_time)
    n_dof = size(q, 1);
    if isfield(imc, 'active_time'), t_start = imc.active_time; else, t_start = 0; end
    
    if current_time < t_start
        Fc = zeros(n_dof, 1); Jc = zeros(n_dof, n_dof);
        Ffr = zeros(n_dof, 1); Jfr = zeros(n_dof, n_dof);
        return; 
    end

    % 若所有冰柱均断裂，直接返回
    if isfield(imc, 'is_broken') && all(imc.is_broken)
        Fc = zeros(n_dof, 1); Jc = zeros(n_dof, n_dof);
        Ffr = zeros(n_dof, 1); Jfr = zeros(n_dof, n_dof);
        return; 
    end

    n_edges = size(rod_edges, 1);
    k_c = imc.k_c; mu_k = imc.mu_k; delta = imc.delta;
    
    % 读取几何与运动参数
    R_ice = imc.ice_radius; 
    R_rod = imc.rod_radius; 
    
    % 判断是否使用了自定义的冰柱坐标矩阵
    use_custom_positions = isfield(imc, 'ice_positions');
    if use_custom_positions
        num_ice = size(imc.ice_positions, 1);
        imc.num_ice = num_ice; 
    else
        num_ice = imc.num_ice;
        R_array = imc.array_radius;
        L_center = imc.array_center_dist;
    end
    
    omega_orbit = imc.omega_mag; % 公转速度
    if isfield(imc, 'omega_spin'), omega_spin = imc.omega_spin; else, omega_spin = 0; end
    
    sigma_t = imc.sigma_t; z_root = imc.z_root; 
    d_ice = 2 * R_ice;
    I_ice = (pi * d_ice^4) / 64; 
    c_ice = d_ice / 2;           

    Fc = zeros(n_dof, 1); Jc = zeros(n_dof, n_dof);
    Ffr = zeros(n_dof, 1); Jfr = zeros(n_dof, n_dof);
    contact_dist_limit = R_ice + R_rod + delta;
    
    % 获取当前的公转累加角度
    if isfield(imc, 'theta_accumulated')
        theta_orb = -imc.theta_accumulated; 
    else
        theta_orb = -omega_orbit * current_time; 
    end

    % 公转旋转矩阵 
    rot_mat = [cos(theta_orb), -sin(theta_orb); sin(theta_orb), cos(theta_orb)];

    %% 遍历所有冰柱
    for j = 1:num_ice
        if isfield(imc, 'is_broken') && imc.is_broken(j), continue; end 
        
        % 计算第 j 根冰柱的当前绝对位置
        if use_custom_positions
            P0_j_xy = imc.ice_positions(j, :)';
        else
            theta_spin = omega_spin * current_time;
            phi_j = 2 * pi * (j - 1) / num_ice + theta_spin;
            P0_j_xy = [L_center + R_array * cos(phi_j); R_array * sin(phi_j)];
        end
        
        P_ice_xy = rot_mat * P0_j_xy; 
        
        % 计算该冰柱中心点由于公转产生的线速度 
        V_orbit_xy = [-omega_orbit * P_ice_xy(2); omega_orbit * P_ice_xy(1)];

        % 【核心修改】：基于边（edge）检测碰撞并分配接触力
        for i = 1:n_edges
            node1_idx = rod_edges(i, 1);
            node2_idx = rod_edges(i, 2);
            
            idx1 = (3*node1_idx-2):(3*node1_idx);
            idx2 = (3*node2_idx-2):(3*node2_idx);
            
            p1 = q(idx1); 
            p2 = q(idx2);
            p1_xy = p1(1:2);
            p2_xy = p2(1:2);
            
            % 计算线段上距离冰柱中心最近的点
            v_edge = p2_xy - p1_xy;
            w = P_ice_xy - p1_xy;
            c1 = dot(w, v_edge);
            c2 = dot(v_edge, v_edge);
            
            if c1 <= 0
                b = 0; % 最近点在节点 1
                pb_xy = p1_xy;
            elseif c2 <= c1
                b = 1; % 最近点在节点 2
                pb_xy = p2_xy;
            else
                b = c1 / c2; % 最近点在线段内部
                pb_xy = p1_xy + b * v_edge;
            end
            
            diff = pb_xy - P_ice_xy; % 最近接触点相对于冰柱的矢量
            dist = norm(diff);
            
            if dist < contact_dist_limit
                % === 【修改代码：提前计算接触点 Z 坐标并拦截】 ===
                z_b = p1(3) + b * (p2(3) - p1(3)); % 计算接触点在 Z 轴的高度
                
                % 检查：如果接触点的高度不在冰柱的实体范围内，则视为“挥空”，跳过不产生力
                if z_b < imc.ice_z_min || z_b > imc.ice_z_max
                    continue; 
                end
                % ==================================================

                if dist < 1e-9, normal = [1;0;0]; else, normal_xy = diff / dist; normal = [normal_xy; 0]; end
                pen_depth = contact_dist_limit - dist;
                f_mag = k_c * pen_depth; 
                
                % 断裂判据：使用接触点的插值 Z 坐标计算力臂
                lever_arm = abs(z_root - z_b); 
                safe_lever = max(lever_arm, 0.001);
                F_break = (sigma_t * I_ice) / (c_ice * safe_lever);
                
                if f_mag >= F_break
                    imc.is_broken(j) = true; 
                    imc.peak_force(j) = f_mag; 
                    break; % 冰柱已断，无需继续计算该冰柱与其他边的碰撞
                end
                
                % 【关键点】：将力通过插值系数分配给两端节点
                F_contact = f_mag * normal;
                Fc(idx1) = Fc(idx1) + (1 - b) * F_contact;
                Fc(idx2) = Fc(idx2) + b * F_contact;
                
                % 【关键点】：基于形函数分配雅可比矩阵
                J_base = -k_c * (normal * normal');
                Jc(idx1, idx1) = Jc(idx1, idx1) + (1 - b)^2 * J_base;
                Jc(idx1, idx2) = Jc(idx1, idx2) + (1 - b) * b * J_base;
                Jc(idx2, idx1) = Jc(idx2, idx1) + b * (1 - b) * J_base;
                Jc(idx2, idx2) = Jc(idx2, idx2) + b^2 * J_base;
                
                if imc.compute_friction
                    % 插值计算接触点当前速度
                    v1 = (q(idx1) - q0(idx1)) / dt;
                    v2 = (q(idx2) - q0(idx2)) / dt;
                    v_b = v1 + b * (v2 - v1);

                    if use_custom_positions
                        V_spin_xy = [0; 0]; 
                    else
                        V_spin_xy = [-omega_spin * diff(2); omega_spin * diff(1)];
                    end
                    
                    V_ice_total = [V_orbit_xy + V_spin_xy; 0];
                    
                    % 计算摩擦力方向并同样按比例分配到节点
                    v_rel = v_b - V_ice_total;
                    v_tan = v_rel - dot(v_rel, normal) * normal;
                    if norm(v_tan) > 1e-6
                        tangent_dir = v_tan / norm(v_tan);
                        f_fr = -mu_k * f_mag * tangent_dir;
                        
                        Ffr(idx1) = Ffr(idx1) + (1 - b) * f_fr;
                        Ffr(idx2) = Ffr(idx2) + b * f_fr;
                    end
                end
            end
        end
    end
    imc.C = []; 
end