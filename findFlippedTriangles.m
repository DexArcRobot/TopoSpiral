function [flipped, degenerate] = findFlippedTriangles(V, F)
    A = V(F(:,1), :);
    B = V(F(:,2), :);
    C = V(F(:,3), :);

    % 有向面积的两倍：正值为逆时针，负值为顺时针
    area2 = (B(:,1)-A(:,1)).*(C(:,2)-A(:,2)) ...
          - (B(:,2)-A(:,2)).*(C(:,1)-A(:,1));

    flipped    = find(area2 < 0);
    degenerate = find(area2 == 0);
end