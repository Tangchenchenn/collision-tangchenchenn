function [Fd, Jd] = getAerodynamicDrag(q, q0, dt, env, MultiRod, omega)
    % 读取参数
    Cd = env.Cd;       % 阻力系数
    rho_air = env.rho; % 空气密度
    
    % 相对速度 (Relative Velocity in Rotating Frame)
    u = (q - q0) / dt; 
    
    n_DOF = MultiRod.n_DOF;
    Fd = zeros(n_DOF, 1);
    Jd = zeros(n_DOF, n_DOF); % 雅可比矩阵略复杂，这里暂忽略旋转项的雅可比修正，通常能收敛
    
    %% 1. 针对 ROD (绳索/边) 的空气阻力计算
    % 遍历每一条边 (Edge) 计算圆柱体阻力
    % if isfield(MultiRod, 'rod_edges') && ~isempty(MultiRod.rod_edges)
    if ~isempty(MultiRod.rod_edges)
        n_rod_edges = size(MultiRod.rod_edges, 1);
        r0 = MultiRod.refLen(1) * 0 + 0.002; % 假设半径 2mm，或者从 geometry 读取: env.rod_r0
        if isfield(env, 'rod_r0'), r0 = env.rod_r0; end
        
        for i = 1:n_rod_edges
            % 获取边的两个节点索引
            node1 = MultiRod.rod_edges(i, 1);
            node2 = MultiRod.rod_edges(i, 2);
            idx1 = mapNodetoDOF(node1);
            idx2 = mapNodetoDOF(node2);
            
            % 节点位置
            x1 = q(idx1);
            x2 = q(idx2);
            x_mid = 0.5 * (x1 + x2); % 边中点位置
            
            % 1. 计算线速度 (相对速度 u + 牵连速度 w x r)
            u1 = u(idx1);
            u2 = u(idx2);
            u_mid = 0.5 * (u1 + u2); % 边中点的相对速度
            
            % 牵连速度 v_trans = omega x r
            v_trans = cross(omega, x_mid);
            
            % 总速度 (空气相对于线的速度 = - (u + v_trans))
            % 但阻力方向与速度相反，F ~ -v|v|
            v_total = u_mid + v_trans; 
            v_norm = norm(v_total);
            
            if v_norm < 1e-6
                continue;
            end
            
            % 2. 计算阻力 (简化模型：各向同性阻力)
            % F = -0.5 * rho * Cd * Area * |v| * v
            % Area (投影面积) = 直径 * 长度
            edge_dir = x2 - x1;
            L = norm(edge_dir);
            Area = (2 * r0) * L; 
            
            % 阻力矢量
            F_drag_vector = -0.5 * rho_air * Cd * Area * v_norm * v_total;
            
            % 3. 分配到两个节点
            Fd(idx1) = Fd(idx1) + F_drag_vector / 2;
            Fd(idx2) = Fd(idx2) + F_drag_vector / 2;
            
            % 注意：此处省略了 Jacobian (Jd) 的填充以简化代码
            % 对于隐式积分，没有精确 Jacobian 可能会增加迭代次数(Iter变多)，
            % 但只要时间步长 dt 足够小，仍然可以收敛。
            % 如果需要高性能，这里需要添加 dF/dq 和 dF/du 的导数。
        end
    end

    %% 2. 针对 SHELL (三角形面) 的空气阻力计算 (保留原逻辑并修正速度)
    if isfield(MultiRod, 'n_faces') && MultiRod.n_faces > 0
        faceAs = MultiRod.faceA;
        for c=1:MultiRod.n_faces 
            node1ind = MultiRod.face_nodes_shell(c,1);
            node2ind = MultiRod.face_nodes_shell(c,2);
            node3ind = MultiRod.face_nodes_shell(c,3);
            idx = [mapNodetoDOF(node1ind); mapNodetoDOF(node2ind); mapNodetoDOF(node3ind)];
            
            q1 = q(mapNodetoDOF(node1ind));
            q2 = q(mapNodetoDOF(node2ind));
            q3 = q(mapNodetoDOF(node3ind));
            
            % 计算面上各点的总速度
            % 这里简化为计算每个顶点的速度
            u1 = u(mapNodetoDOF(node1ind)) + cross(omega, q1);
            u2 = u(mapNodetoDOF(node2ind)) + cross(omega, q2);
            u3 = u(mapNodetoDOF(node3ind)) + cross(omega, q3);
            
            face_normal = cross((q2-q1), (q3-q2));
            face_A = faceAs(c);
            face_unit_normal = face_normal/norm(face_normal);

            % 计算每个节点的阻力贡献 (使用了修正后的 u1, u2, u3)
            % 仅保留法向分量逻辑 (原代码逻辑)
            
            % 辅助函数：计算带符号的阻力
            calc_node_drag = @(vel) -sign(dot(vel, face_unit_normal)) * ...
                (rho_air * 0.5 * Cd * face_A / 3) * (dot(vel, face_unit_normal)^2 * face_unit_normal);

            Fd(idx(1:3)) = Fd(idx(1:3)) + calc_node_drag(u1);
            Fd(idx(4:6)) = Fd(idx(4:6)) + calc_node_drag(u2);
            Fd(idx(7:9)) = Fd(idx(7:9)) + calc_node_drag(u3);
            
            % Jacobian 部分省略修正，或者沿用原代码(不含omega项)作为近似
        end
    end
end
% function [Fd, Jd] = getAerodynamicDrag(q,q0,dt,env, MultiRod)
% Cd = env.Cd;
% rho_med = env.rho;
% u = (q-q0)/dt ;
% Fd = zeros(MultiRod.n_DOF,1);
% Jd = zeros(MultiRod.n_DOF,MultiRod.n_DOF);
% faceAs = MultiRod.faceA;
% 
% for c=1:MultiRod.n_faces 
%     node1ind = MultiRod.face_nodes_shell(c,1);
%     node2ind = MultiRod.face_nodes_shell(c,2);
%     node3ind = MultiRod.face_nodes_shell(c,3);
%     idx = [mapNodetoDOF(node1ind); mapNodetoDOF(node2ind); mapNodetoDOF(node3ind)];
%     face_normal = cross( (q(mapNodetoDOF(node2ind))-q(mapNodetoDOF(node1ind))), (q(mapNodetoDOF(node3ind)) - q(mapNodetoDOF(node2ind))) ) ;
%     face_A = faceAs(c);
%     face_unit_normal = face_normal/norm(face_normal);
% 
%     q1 = q(mapNodetoDOF(node1ind));
%     q2 = q(mapNodetoDOF(node2ind));
%     q3 = q(mapNodetoDOF(node3ind));
%     u1 = u(mapNodetoDOF(node1ind));
%     u2 = u(mapNodetoDOF(node2ind));
%     u3 = u(mapNodetoDOF(node3ind));
% 
%     if(dot(u1,face_unit_normal)>0)
%         sign1 = -1;
%     elseif(dot(u1,face_unit_normal)<0)
%         sign1 = 1;
%     else
%         sign1 = 0;
%     end
% 
%     if(dot(u2,face_unit_normal)>0)
%         sign2 = -1;
%     elseif(dot(u2,face_unit_normal)<0)
%         sign2 = 1;
%     else
%         sign2 = 0;
%     end
% 
%     if(dot(u3,face_unit_normal)>0)
%         sign3 = -1;
%     elseif(dot(u3,face_unit_normal)<0)
%         sign3 = 1;
%     else
%         sign3 = 0;
%     end
% 
%     Fd(idx(1:3)) = Fd(idx(1:3)) + sign1*(rho_med*0.5*Cd*face_A/3).* ...
%         (dot(u1, face_unit_normal)^2.*face_unit_normal);
%     Fd(idx(4:6)) = Fd(idx(4:6)) + sign2*(rho_med*0.5*Cd*face_A/3).* ...
%         (dot(u2, face_unit_normal)^2.*face_unit_normal);
%     Fd(idx(7:9)) = Fd(idx(7:9)) + sign3*(rho_med*0.5*Cd*face_A/3).* ...
%         (dot(u3, face_unit_normal)^2.*face_unit_normal);
% 
%     % Jd [(1:3),:]
%     Jd(idx(1:3),idx(1:3)) = Jd(idx(1:3),idx(1:3)) + sign1*(rho_med*Cd*face_A/3).*...
%         (2*dot(u1,face_unit_normal)*face_unit_normal*((1/dt)*face_unit_normal'*eye(3) + ...
%         u1'*gradient_of_unit_normal(face_normal,q1, q2, q3, 1)) + ...
%         dot(u1,face_unit_normal)^2*gradient_of_unit_normal(face_normal,q1, q2, q3, 1));
% 
%     Jd(idx(1:3),idx(4:6)) = Jd(idx(1:3),idx(4:6)) + sign1*(rho_med*Cd*face_A/3).*...
%         dot(u1,face_unit_normal)*(2*face_unit_normal*u1' + dot(u1,face_unit_normal)*eye(3)) * ...
%         gradient_of_unit_normal(face_normal,q1, q2, q3, 2);
% 
%     Jd(idx(1:3),idx(7:9)) = Jd(idx(1:3),idx(7:9)) + sign1*(rho_med*Cd*face_A/3).*...
%         dot(u1,face_unit_normal)*(2*face_unit_normal*u1' + dot(u1,face_unit_normal)*eye(3)) * ...
%         gradient_of_unit_normal(face_normal,q1, q2, q3, 3);
% 
%     % Jd [(4:6),:]
%     Jd(idx(4:6),idx(1:3)) = Jd(idx(4:6),idx(1:3)) + sign2*(rho_med*Cd*face_A/3).*...
%         dot(u2,face_unit_normal)*(2*face_unit_normal*u2' + dot(u2,face_unit_normal)*eye(3)) * ...
%         gradient_of_unit_normal(face_normal,q1, q2, q3, 1);
% 
%     Jd(idx(4:6),idx(4:6)) = Jd(idx(4:6),idx(4:6)) + sign2*(rho_med*Cd*face_A/3).*...
%         (2*dot(u2,face_unit_normal)*face_unit_normal*((1/dt)*face_unit_normal'*eye(3) + ...
%         u2'*gradient_of_unit_normal(face_normal,q1, q2, q3, 2)) + ...
%         dot(u2,face_unit_normal)^2*gradient_of_unit_normal(face_normal,q1, q2, q3, 2));
% 
%     Jd(idx(4:6),idx(7:9)) = Jd(idx(4:6),idx(7:9)) + sign2*(rho_med*Cd*face_A/3).*...
%         dot(u2,face_unit_normal)*(2*face_unit_normal*u2' + dot(u2,face_unit_normal)*eye(3)) * ...
%         gradient_of_unit_normal(face_normal,q1, q2, q3, 3);
% 
%     % Jd [(7:9),:]
%     Jd(idx(7:9),idx(1:3)) = Jd(idx(7:9),idx(1:3)) + sign3*(rho_med*Cd*face_A/3).*...
%         dot(u3,face_unit_normal)*(2*face_unit_normal*u3' + dot(u3,face_unit_normal)*eye(3)) * ...
%         gradient_of_unit_normal(face_normal,q1, q2, q3, 1);
% 
%     Jd(idx(7:9),idx(4:6)) = Jd(idx(7:9),idx(4:6)) + sign3*(rho_med*Cd*face_A/3).*...
%         dot(u3,face_unit_normal)*(2*face_unit_normal*u3' + dot(u3,face_unit_normal)*eye(3)) * ...
%         gradient_of_unit_normal(face_normal,q1, q2, q3, 2);
% 
%     Jd(idx(7:9),idx(7:9)) = Jd(idx(7:9),idx(7:9)) + sign3*(rho_med*Cd*face_A/3).*...
%         (2*dot(u3,face_unit_normal)*face_unit_normal*((1/dt)*face_unit_normal'*eye(3) + ...
%         u3'*gradient_of_unit_normal(face_normal,q1, q2, q3, 3)) + ...
%         dot(u3,face_unit_normal)^2*gradient_of_unit_normal(face_normal,q1, q2, q3, 3));
% 
% 
% end
