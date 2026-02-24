function plot_MultiRod(MultiRod, ctime, sim_params, environment, imc)

    n_nodes = MultiRod.n_nodes;
    q = MultiRod.q;
    edges = MultiRod.Edges;
    n_edges = MultiRod.n_edges;
    n_edges_dof = MultiRod.n_edges_dof;

    % 1. 准备绳索数据 (在随动坐标系下，q即为相对位置，直接绘制即可)
    x1 = q(1:3:3*n_nodes);
    x2 = q(2:3:3*n_nodes);
    x3 = q(3:3:3*n_nodes);

    % 计算总长用于缩放箭头 (视觉优化)
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
    clf; % 清除上一帧
    hold on;

    % 设置坐标轴范围
    if isfield(sim_params, 'plot_x'), xlim(sim_params.plot_x); end
    if isfield(sim_params, 'plot_y'), ylim(sim_params.plot_y); end
    if isfield(sim_params, 'plot_z'), zlim(sim_params.plot_z); end

    % %% --- 3. [核心新增] 绘制旋转的冰柱 (随动系视角) ---
    % % 逻辑：绳子不动，冰柱以 -omega 速度绕中心公转
    % 
    % % A. 读取几何参数 (优先 imc，其次 environment)
    % r_ice = 0.004; d_ice = 0.15; % 默认值
    % if isfield(imc, 'ice_radius'), r_ice = imc.ice_radius; 
    % elseif isfield(environment, 'contact_params'), r_ice = environment.contact_params.ice_radius; end
    % 
    % if isfield(imc, 'ice_center_dist'), d_ice = imc.ice_center_dist; 
    % elseif isfield(environment, 'contact_params'), d_ice = environment.contact_params.ice_center_dist; end
    % 
    % % B. 计算旋转角度
    % omega_val = 0;
    % if isfield(sim_params, 'omega'), omega_val = norm(sim_params.omega); end
    % % 随动系下冰柱反向旋转
    % ice_angle = -omega_val * ctime; 
    % 
    % % C. 生成并旋转圆柱体网格
    % [xc, yc, zc] = cylinder(r_ice, 20);
    % zc = zc * 0.2 - 0.1; % 高度 [-0.1, 0.1]
    % 
    % % 初始位置：X轴正方向 d_ice 处
    % xc = xc + d_ice; 
    % 
    % % 应用旋转矩阵 (绕 Z 轴旋转 ice_angle)
    % % x' = x*cos(t) - y*sin(t)
    % % y' = x*sin(t) + y*cos(t)
    % xc_rot = xc * cos(ice_angle) - yc * sin(ice_angle);
    % yc_rot = xc * sin(ice_angle) + yc * cos(ice_angle);
    % 
    % % D. 绘制冰柱
    % surf(xc_rot, yc_rot, zc, 'FaceColor', [0 1 1], 'EdgeColor', 'none', 'FaceAlpha', 0.6);

%% --- 3. [核心新增] 绘制旋转的冰柱 (随动系视角) ---
    % 逻辑：绳子不动，冰柱以 -omega 速度绕中心公转
    
    % 【新增判断】：只有在当前帧的时间 (ctime) 小于断裂时间时，才绘制冰柱
    is_broken_now = isfield(imc, 'breaking_time') && (ctime >= imc.breaking_time);
    
    if ~is_broken_now
        % A. 读取几何参数 (优先 imc，其次 environment)
        r_ice = 0.004; d_ice = 0.15; % 默认值
        if isfield(imc, 'ice_radius'), r_ice = imc.ice_radius; 
        elseif isfield(environment, 'contact_params'), r_ice = environment.contact_params.ice_radius; end
        
        if isfield(imc, 'ice_center_dist'), d_ice = imc.ice_center_dist; 
        elseif isfield(environment, 'contact_params'), d_ice = environment.contact_params.ice_center_dist; end

        % B. 计算旋转角度
        omega_val = 0;
        if isfield(sim_params, 'omega'), omega_val = norm(sim_params.omega); end
        % 随动系下冰柱反向旋转
        ice_angle = -omega_val * ctime; 

        % C. 生成并旋转圆柱体网格
        [xc, yc, zc] = cylinder(r_ice, 20);
        zc = zc * 0.2 - 0.1; % 高度 [-0.1, 0.1]
        
        % 初始位置：X轴正方向 d_ice 处
        xc = xc + d_ice; 
        
        % 应用旋转矩阵 (绕 Z 轴旋转 ice_angle)
        % x' = x*cos(t) - y*sin(t)
        % y' = x*sin(t) + y*cos(t)
        xc_rot = xc * cos(ice_angle) - yc * sin(ice_angle);
        yc_rot = xc * sin(ice_angle) + yc * cos(ice_angle);

        % D. 绘制冰柱
        surf(xc_rot, yc_rot, zc, 'FaceColor', [0 1 1], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    end

    %% 4. 绘制绳索
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
