function V_rot = rotatePoints(V, angleDeg)
% rotatePoints 将二维点绕原点旋转
%   V        : N×2 点坐标 [x, y]
%   angleDeg : 旋转角度（度），正值逆时针，负值顺时针
%   V_rot    : 旋转后的 N×2 点坐标

    validateattributes(V, {'numeric'}, {'2d', 'ncols', 2});
    validateattributes(angleDeg, {'numeric'}, ...
                       {'scalar', 'real', 'finite'});

    R = [cosd(angleDeg), -sind(angleDeg);
         sind(angleDeg),  cosd(angleDeg)];

    V_rot = double(V) * R.';
end