function idx = getBoundary1Ring(f, v)
% GETBOUNDARY1RING 获取边界点及其1邻域点
%
% 输入:
%   f   - M×3 三角形面片顶点索引
%   v   - N×3 或 N×2 顶点坐标
%
% 输出:
%   idx - n×1 顶点编号（边界点 + 1邻域点）

    N = size(v, 1);

    % 1. 提取全部三角形边
    E = [f(:,[1 2]);
         f(:,[2 3]);
         f(:,[3 1])];

    E = sort(E, 2);

    % 2. 统计每条边出现的次数
    [uniqueE, ~, ic] = unique(E, 'rows');
    counts = accumarray(ic, 1);

    % 3. 只出现一次的边为边界边
    boundaryEdges = uniqueE(counts == 1, :);

    % 4. 提取边界顶点
    boundaryVerts = unique(boundaryEdges(:));

    % 5. 构造顶点邻接矩阵
    A = sparse([E(:,1); E(:,2)], ...
               [E(:,2); E(:,1)], ...
               1, N, N);

    % 6. 找出边界顶点的1邻域顶点
    neighborVerts = find(any(A(:, boundaryVerts), 2));

    % 7. 合并、去重、转为列向量
    idx = unique([boundaryVerts; neighborVerts]);

end