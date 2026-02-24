function [Fc, Jc, Ffr, Jfr, imc] = IceContact(imc, q, q0, rod_edges, iter, dt, current_time)
    n_dof = size(q, 1);
    if isfield(imc, 'active_time')
        t_start = imc.active_time;
    else
        t_start = 0; % 默认一开始就存在
    end
    
    if current_time < t_start
        Fc = zeros(n_dof, 1); Jc = zeros(n_dof, n_dof);
        Ffr = zeros(n_dof, 1); Jfr = zeros(n_dof, n_dof);
        return; % 提前结束函数
    end
    %% 0. [新增] 检查全局破碎状态
    % 如果上一帧已经碎了，这帧直接返回0力，不再计算
    if isfield(imc, 'is_broken') && imc.is_broken
        n_dof = size(q, 1);
        Fc = zeros(n_dof, 1); Jc = zeros(n_dof, n_dof);
        Ffr = zeros(n_dof, 1); Jfr = zeros(n_dof, n_dof);
        return; 
    end

    %% 1. 参数提取
    n_dof = size(q, 1);
    n_nodes = n_dof / 3;
    
    k_c = imc.k_c;          
    mu_k = imc.mu_k;       
    delta = imc.delta;
    
    if isfield(imc, 'ice_radius'), R_ice = imc.ice_radius; else, R_ice = 0.004; end
    if isfield(imc, 'rod_radius'), R_rod = imc.rod_radius; else, R_rod = 0.002; end
    if isfield(imc, 'ice_center_dist'), L_dist = imc.ice_center_dist; else, L_dist = 0.15; end
    if isfield(imc, 'omega_mag'), omega_scalar = imc.omega_mag; else, omega_scalar = 0; end

    % [新增] 读取力学判据参数
    if isfield(imc, 'sigma_t'), sigma_t = imc.sigma_t; else, sigma_t = 1.5e6; end
    if isfield(imc, 'z_root'), z_root = imc.z_root; else, z_root = 0.1; end

    % [新增] 预计算截面属性 (假设圆形截面) [cite: 7, 8]
    d_ice = 2 * R_ice;
    I_ice = (pi * d_ice^4) / 64; % 惯性矩
    c_ice = d_ice / 2;           % 距离中性轴最远距离
    % 
    % %% 2. 冰柱运动计算 (保持不变)
    % theta_ice = -omega_scalar * current_time;
    % P_ice_xy = [L_dist * cos(theta_ice); L_dist * sin(theta_ice)];
    % V_ice_xy = [-omega_scalar * P_ice_xy(2); omega_scalar * P_ice_xy(1)];
    % V_ice = [V_ice_xy; 0]; 
    %% 2. 冰柱运动计算 (修正位置逻辑)
    % 使用从 main.m 传进来的累加角度，确保加速过程位置精确
    if isfield(imc, 'theta_accumulated')
        theta_ice = -imc.theta_accumulated; 
    else
        theta_ice = -imc.omega_mag * current_time; % 兜底逻辑
    end

    % 基于当前角度计算冰柱中心坐标
    P_ice_xy = [L_dist * cos(theta_ice); L_dist * sin(theta_ice)]; %
    
    % 基于当前瞬时转速计算冰柱切向速度 (用于摩擦力)
    V_ice_xy = [-imc.omega_mag * P_ice_xy(2); imc.omega_mag * P_ice_xy(1)]; %
    V_ice = [V_ice_xy; 0];
    %% 3. 初始化输出 (保持不变)
    Fc = zeros(n_dof, 1); Jc = zeros(n_dof, n_dof);
    Ffr = zeros(n_dof, 1); Jfr = zeros(n_dof, n_dof);
    
    contact_dist_limit = R_ice + R_rod + delta;
    
    %% 4. 遍历检测
    for i = 1:n_nodes
        idx = (3*i-2):(3*i);
        p_node = q(idx);        
        v_node = (q(idx) - q0(idx)) / dt; 
        
        p_node_xy = p_node(1:2);
        diff = p_node_xy - P_ice_xy; 
        dist = norm(diff);
        
        if dist < contact_dist_limit
            if dist < 1e-9, normal = [1;0;0]; else, normal_xy = diff / dist; normal = [normal_xy; 0]; end
            pen_depth = contact_dist_limit - dist;
            f_mag = k_c * pen_depth; % 当前计算出的弹性接触力
            
            % ================= [核心修改] 破碎判据 =================
            % 1. 计算力臂 (接触点到根部的高度差) [cite: 5]
            % p_node(3) 是节点的 Z 坐标
            lever_arm = abs(z_root - p_node(3)); 
            
            % 2. 计算临界断裂力 F_break [cite: 15]
            % F_break = (sigma_t * I) / (c * lever_arm)
            % 避免 lever_arm 为 0 (假设最小臂长 1mm)
            safe_lever = max(lever_arm, 0.001);
            F_break = (sigma_t * I_ice) / (c_ice * safe_lever);
            

            % ... 前面计算 f_mag 的代码 ...
            
            % 3. 判断是否达到断裂条件
            if f_mag >= F_break
                % [核心修改] 记录峰值力到 imc 结构体中
                imc.is_broken = true; 
                imc.peak_force = f_mag; % 将瞬时力保存起来
                
                fprintf('>>> 冰柱断裂! 时间:%.4fs | 峰值力: %.1fN | 位置: %.3fm <<<\n', ...
                        current_time, f_mag, p_node(3));
                
                % 为了让主程序记录到这个力，断裂这一步不应该直接清零
                % 我们让它先计算完这一步的 Fc，从下一步开始才会因为 is_broken 而返回 0
            end
            
            % 正常施加力（即使断裂这一帧也要施加，否则主程序收不到值）
            Fc(idx) = Fc(idx) + f_mag * normal;
            J_node = -k_c * (normal * normal');
            Jc(idx, idx) = Jc(idx, idx) + J_node;
            % ... 后面逻辑保持不变 ...


    %         % 3. 判断是否达到断裂条件
    %         if f_mag >= F_break
    %             % [断裂发生]
    %             imc.is_broken = true; % 永久标记破碎状态
    % 
    %             fprintf('>>> 冰柱断裂! Time:%.4fs | Force: %.1fN > Limit: %.1fN | Z-Loc: %.3fm <<<\n', ...
    %                     current_time, f_mag, F_break, p_node(3));
    % 
    %             % % 立即将所有力清零并返回，因为冰柱消失了
    %             Fc(:) = 0; Jc(:) = 0; Ffr(:) = 0; Jfr(:) = 0;
    %             return; 
    %             % [调试] 先不要立即 return，先打印出来看看数值是多少
    % % fprintf('>>> ⚠️ 判定破碎! Time:%.4fs | Force: %.2fN > Limit: %.2fN | Z: %.3fm\n', ...
    %         % current_time, f_mag, F_break, p_node(3));
    % 
    % % imc.is_broken = true;  % <--- 暂时注释掉这一行！
    % % Fc(:) = 0; ...         % <--- 暂时注释掉这一行！
    % % return;                % <--- 暂时注释掉这一行！
    %         end
    %         % ====================================================
    % 
    %         % 如果没断，正常施加力
    %         Fc(idx) = Fc(idx) + f_mag * normal;
    %         J_node = -k_c * (normal * normal');
    %         Jc(idx, idx) = Jc(idx, idx) + J_node;
    % 
            if imc.compute_friction
                v_rel = v_node - V_ice;
                v_tan = v_rel - dot(v_rel, normal) * normal;
                if norm(v_tan) > 1e-6
                    tangent_dir = v_tan / norm(v_tan);
                    f_fr = -mu_k * f_mag * tangent_dir;
                    Ffr(idx) = Ffr(idx) + f_fr;
                end
            end
        end
    end
    imc.C = []; 
end