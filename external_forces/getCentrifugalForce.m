function [Fcen, Jcen] = getCentrifugalForce(MultiRod, q, omega)
    n_nodes = MultiRod.n_nodes;
    n_DOF = MultiRod.n_DOF;
    Fcen = zeros(n_DOF, 1);
    Jcen = zeros(n_DOF, n_DOF);

    % 旋转算子矩阵 (omega x r) 对应的反对称矩阵 Omega
    Omega = [0, -omega(3), omega(2); omega(3), 0, -omega(1); -omega(2), omega(1), 0];
    Omega2 = Omega * Omega; % 对应 omega x (omega x r)

    for i = 1:n_nodes
        idx = (3*i-2):(3*i);
        m_node = MultiRod.massVec(idx(1)); % 获取该节点的平移质量
        pos = q(idx); % 当前节点在随动系下的坐标

        % 离心力 F = -m * Omega^2 * pos
        F_node = - m_node * (Omega2 * pos);
        Fcen(idx) = F_node;

        % 离心刚度项 (Jacobian) J = dF/dq = -m * Omega^2
        % 注意：在隐式步进器中 J 是减去的，所以这里符号需与 F 对齐
        Jcen(idx, idx) = - m_node * Omega2;
    end
end
% function [Fcen, Jcen] = getCentrifugalForce(MultiRod, q, omega)
%     n_nodes = MultiRod.n_nodes;
%     n_DOF = MultiRod.n_DOF;
% 
%     % 1. 预计算 Omega平方矩阵 (3x3)
%     % Omega = [0, -omega(3), omega(2); omega(3), 0, -omega(1); -omega(2), omega(1), 0];
%     % Omega2 = Omega * Omega; 
%     % 手动展开计算 Omega^2 通常更快且更精确
%     w1 = omega(1); w2 = omega(2); w3 = omega(3);
%     Omega2 = [ -w2^2-w3^2,   w1*w2,      w1*w3;
%                 w1*w2,      -w1^2-w3^2,  w2*w3;
%                 w1*w3,       w2*w3,     -w1^2-w2^2 ];
% 
%     % 2. 提取所有节点的质量 (假设每个节点3个DOF质量相同，取每隔3个的值)
%     % MultiRod.massVec 是 n_DOF x 1 向量
%     m_nodes = MultiRod.massVec(1:3:end); 
% 
%     % 3. 向量化计算力 F = -m * Omega^2 * r
%     % 将 q 重塑为 [3, n_nodes]
%     q_reshaped = reshape(q(1:3*n_nodes), 3, n_nodes);
% 
%     % 计算每个节点的力 (3 x n_nodes)
%     % 利用广播机制: F_matrix = -Omega2 * q * m
%     % 注意：m_nodes 需要转置匹配乘法，或者逐列相乘
%     F_matrix = -Omega2 * q_reshaped .* m_nodes'; 
% 
%     % 重塑回列向量并填充到 Fcen
%     Fcen = zeros(n_DOF, 1);
%     Fcen(1:3*n_nodes) = F_matrix(:);
% 
%     % 4. 向量化构建稀疏 Jacobian
%     % Jacobian 是块对角矩阵，每个块是 -m_i * Omega2
% 
%     % 构建重复的 Omega2 序列
%     % 使用 kron 将质量向量扩展到对角块
%     % 方法：利用 Kronecker 积构建稀疏矩阵
% 
%     % 这是一个高效构建块对角矩阵的方法：
%     % J_blocks = -m_i * Omega2
%     % 我们构建一个大矩阵，实际上是 mass_matrix * (I kron Omega2) 的变体
% 
%     % 简单做法：利用 sparse 和 kron
%     % J = - kron(diag(m_nodes), Omega2);
%     Jcen_small = -kron(spdiags(m_nodes, 0, n_nodes, n_nodes), Omega2);
% 
%     Jcen = sparse(n_DOF, n_DOF);
%     Jcen(1:3*n_nodes, 1:3*n_nodes) = Jcen_small;
% end