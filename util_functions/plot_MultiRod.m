function plot_MultiRod(MultiRod, ctime, sim_params, environment, imc)

    n_nodes = MultiRod.n_nodes;
    q = MultiRod.q;
    edges = MultiRod.Edges;
    n_edges = MultiRod.n_edges;
    n_edges_dof = MultiRod.n_edges_dof;

    % 1. 准备绳索数据
    x1 = q(1:3:3*n_nodes);
    x2 = q(2:3:3*n_nodes);
    x3 = q(3:3:3*n_nodes);

    L = sum(sqrt(diff(x1).^2 + diff(x2).^2 + diff(x3).^2));
    scale_factor = 0.1 * L;
    
    a1 = scale_factor * MultiRod.a1;
    a2 = scale_factor * MultiRod.a2;
    m1 = scale_factor * MultiRod.m1;
    m2 = scale_factor * MultiRod.m2;

    % 2. 准备绘图窗口
    if isempty(get(groot,'CurrentFigure'))
        figure(2);
    end
    clf; 
    hold on;

    % 设置坐标轴范围
    if isfield(sim_params, 'plot_x'), xlim(sim_params.plot_x); end
    if isfield(sim_params, 'plot_y'), ylim(sim_params.plot_y); end
    if isfield(sim_params, 'plot_z'), zlim(sim_params.plot_z); end

    %% 3. 绘制冰柱组
    if isfield(imc, 'num_ice') || isfield(imc, 'ice_positions')
        
        use_custom_positions = isfield(imc, 'ice_positions');
        if use_custom_positions
            num_ice = size(imc.ice_positions, 1);
        else
            num_ice = imc.num_ice;
            R_array = imc.array_radius;
            L_center = imc.array_center_dist;
        end
        
        R_ice = imc.ice_radius;
        if isfield(imc, 'z_root'), z_root = imc.z_root; else, z_root = 0.07; end

        % --- [计算公转角度] ---
        if isfield(imc, 'theta_accumulated')
            theta_orb = -imc.theta_accumulated;
        else
            theta_orb = -imc.omega_mag * ctime;
        end
        rot_mat = [cos(theta_orb), -sin(theta_orb); sin(theta_orb), cos(theta_orb)];
        
        % --- [计算自转角度 (仅用于环形阵列)] ---
        if isfield(imc, 'omega_spin') && ~use_custom_positions
            theta_spin = imc.omega_spin * ctime; 
        else
            theta_spin = 0;
        end

        % 生成圆柱体网格
        [X_cyl, Y_cyl, Z_cyl] = cylinder(R_ice, 20);
        Z_cyl = Z_cyl * z_root; 

        % 遍历并绘制每一根冰柱
        for j = 1:num_ice
            if ~isfield(imc, 'is_broken') || ~imc.is_broken(j)
                
                % 【核心修改】：提取当前冰柱的初始 (x, y) 坐标
                if use_custom_positions
                    P0_j_xy = imc.ice_positions(j, :)';
                else
                    phi_j = 2 * pi * (j - 1) / num_ice + theta_spin; 
                    P0_j_xy = [L_center + R_array * cos(phi_j); R_array * sin(phi_j)];
                end
                
                % 绕主轴公转
                P_ice_xy = rot_mat * P0_j_xy;

                % 平移网格数据
                X_plot = X_cyl + P_ice_xy(1);
                Y_plot = Y_cyl + P_ice_xy(2);

                % 绘制 3D 表面 
                surf(X_plot, Y_plot, Z_cyl, 'FaceColor', [0 1 1], ...
                    'EdgeColor', 'k', 'EdgeAlpha', 0.2, 'FaceAlpha', 0.6); 
            end
        end
    end
    
    %% 4. 绘制绳索
    for i = 1:n_edges
        n1 = edges(i,1);
        n2 = edges(i,2);
        n1pos = q(3*n1-2:3*n1);
        n2pos = q(3*n2-2:3*n2);
        
        plot3([n1pos(1); n2pos(1)], [n1pos(2); n2pos(2)], [n1pos(3); n2pos(3)], 'ko-', 'LineWidth', 1.5, 'MarkerSize', 3);

        % 绘制固定点 (红色)
        if ~isempty(MultiRod.fixed_nodes)
            if ismember(n1, MultiRod.fixed_nodes), plot3(n1pos(1), n1pos(2), n1pos(3), 'ro', 'MarkerFaceColor', 'r'); end
            if ismember(n2, MultiRod.fixed_nodes), plot3(n2pos(1), n2pos(2), n2pos(3), 'ro', 'MarkerFaceColor', 'r'); end
        end
    end

    % 5. 绘制 Bishop Frame (可选)
    if sim_params.showFrames
        for c = 1:n_edges_dof
            n1 = edges(c,1); n2 = edges(c,2);
            xp = (q(3*n1-2:3*n1) + q(3*n2-2:3*n2)) / 2;
            plot3([xp(1), xp(1)+m1(c,1)], [xp(2), xp(2)+m1(c,2)], [xp(3), xp(3)+m1(c,3)], 'r-');
            plot3([xp(1), xp(1)+m2(c,1)], [xp(2), xp(2)+m2(c,2)], [xp(3), xp(3)+m2(c,3)], 'g-');
        end
    end

    % 6. 收尾
    hold off;
    title(sprintf('Time: %.4f s (Rotating Frame)', ctime));
    
    if isfield(sim_params, 'view')
        if sim_params.view == "xy", view(2);
        elseif sim_params.view == "xz", view(0, 90);
        else, view(3); end
    else
        view(3);
    end
    
    xlabel('X_R (m)'); ylabel('Y_R (m)'); zlabel('Z_R (m)');
    grid on; axis equal;
    drawnow;
end