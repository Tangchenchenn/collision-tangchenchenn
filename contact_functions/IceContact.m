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

    n_nodes = n_dof / 3;
    k_c = imc.k_c; mu_k = imc.mu_k; delta = imc.delta;
    
    % 读取几何与运动参数
    R_ice = imc.ice_radius; 
    R_rod = imc.rod_radius; 
    
    % 【核心修改 1】：判断是否使用了自定义的冰柱坐标矩阵
    use_custom_positions = isfield(imc, 'ice_positions');
    if use_custom_positions
        num_ice = size(imc.ice_positions, 1);
        imc.num_ice = num_ice; % 强制同步数量
    else
        num_ice = imc.num_ice;
        R_array = imc.array_radius;
        L_center = imc.array_center_dist;
    end
    
    omega_orbit = imc.omega_mag; % 公转速度 (主轴转速)
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

    % 公转旋转矩阵 (用于将初始坐标旋转到当前时刻)
    rot_mat = [cos(theta_orb), -sin(theta_orb); sin(theta_orb), cos(theta_orb)];

    %% 遍历所有冰柱
    for j = 1:num_ice
        if isfield(imc, 'is_broken') && imc.is_broken(j), continue; end % 已断裂则跳过
        
        % 【核心修改 2】：计算第 j 根冰柱的当前绝对位置
        if use_custom_positions
            % 直接读取自定义的绝对坐标
            P0_j_xy = imc.ice_positions(j, :)';
        else
            % 原有的环形阵列逻辑
            theta_spin = omega_spin * current_time;
            phi_j = 2 * pi * (j - 1) / num_ice + theta_spin;
            P0_j_xy = [L_center + R_array * cos(phi_j); R_array * sin(phi_j)];
        end
        
        P_ice_xy = rot_mat * P0_j_xy; % 当前时刻的公转位置
        
        % 计算该冰柱中心点由于公转产生的线速度 (v = w × r)
        V_orbit_xy = [-omega_orbit * P_ice_xy(2); omega_orbit * P_ice_xy(1)];

        % 检测与所有节点的碰撞
        for i = 1:n_nodes
            idx = (3*i-2):(3*i);
            p_node = q(idx);        
            v_node = (q(idx) - q0(idx)) / dt; 
            
            p_node_xy = p_node(1:2);
            diff = p_node_xy - P_ice_xy; % 节点相对冰柱中心的矢量
            dist = norm(diff);
            
            if dist < contact_dist_limit
                if dist < 1e-9, normal = [1;0;0]; else, normal_xy = diff / dist; normal = [normal_xy; 0]; end
                pen_depth = contact_dist_limit - dist;
                f_mag = k_c * pen_depth; 
                
                % 断裂判据
                lever_arm = abs(z_root - p_node(3)); 
                safe_lever = max(lever_arm, 0.001);
                F_break = (sigma_t * I_ice) / (c_ice * safe_lever);
                
                if f_mag >= F_break
                    imc.is_broken(j) = true; 
                    imc.peak_force(j) = f_mag; 
                end
                
                Fc(idx) = Fc(idx) + f_mag * normal;
                J_node = -k_c * (normal * normal');
                Jc(idx, idx) = Jc(idx, idx) + J_node;
                
                if imc.compute_friction
                    % 运动学叠加：自转线速度 + 公转线速度
                    if use_custom_positions
                        V_spin_xy = [0; 0]; % 离散分布不再计算阵列自转
                    else
                        V_spin_xy = [-omega_spin * diff(2); omega_spin * diff(1)];
                    end
                    
                    % 表面接触点的绝对速度
                    V_ice_total = [V_orbit_xy + V_spin_xy; 0];
                    
                    % 相对速度与切向摩擦
                    v_rel = v_node - V_ice_total;
                    v_tan = v_rel - dot(v_rel, normal) * normal;
                    if norm(v_tan) > 1e-6
                        tangent_dir = v_tan / norm(v_tan);
                        f_fr = -mu_k * f_mag * tangent_dir;
                        Ffr(idx) = Ffr(idx) + f_fr;
                    end
                end
            end
        end
    end
    imc.C = []; 
end