function [Fcor, Jcor] = getCoriolisForce(MultiRod, q, q0, dt, omega)
    n_nodes = MultiRod.n_nodes;
    n_DOF = MultiRod.n_DOF;
    Fcor = zeros(n_DOF, 1);
    Jcor = zeros(n_DOF, n_DOF);

    Omega = [0, -omega(3), omega(2); omega(3), 0, -omega(1); -omega(2), omega(1), 0];
    v_rel = (q - q0) / dt; % 相对随动系的速度

    for i = 1:n_nodes
        idx = (3*i-2):(3*i);
        m_node = MultiRod.massVec(idx(1));

        % 科氏力 F = -2 * m * (omega x v_rel)
        F_node = -2 * m_node * (Omega * v_rel(idx));
        Fcor(idx) = F_node;

        % 雅可比项 J = dF/dq = dF/dv * dv/dq = (-2 * m * Omega) * (1/dt)
        Jcor(idx, idx) = -2 * m_node * Omega / dt;
    end
end
% function [Fcor, Jcor] = getCoriolisForce(MultiRod, q, q0, dt, omega)
%     n_nodes = MultiRod.n_nodes;
%     n_DOF = MultiRod.n_DOF;
% 
%     Omega = [0, -omega(3), omega(2); omega(3), 0, -omega(1); -omega(2), omega(1), 0];
% 
%     % 1. 向量化速度
%     v_rel = (q(1:3*n_nodes) - q0(1:3*n_nodes)) / dt;
%     v_reshaped = reshape(v_rel, 3, n_nodes);
% 
%     % 2. 质量
%     m_nodes = MultiRod.massVec(1:3:end);
% 
%     % 3. 计算力 F = -2 * m * (omega x v)
%     % omega x v 等价于 Omega * v
%     F_matrix = -2 * (Omega * v_reshaped) .* m_nodes';
% 
%     Fcor = zeros(n_DOF, 1);
%     Fcor(1:3*n_nodes) = F_matrix(:);
% 
%     % 4. Jacobian
%     % J = -2 * m * Omega * (1/dt)
%     factor = -2 / dt;
%     Jcor_small = factor * kron(spdiags(m_nodes, 0, n_nodes, n_nodes), Omega);
% 
%     Jcor = sparse(n_DOF, n_DOF);
%     Jcor(1:3*n_nodes, 1:3*n_nodes) = Jcor_small;
% end