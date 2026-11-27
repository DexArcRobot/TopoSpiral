function [gradientVec, gradientTangent, normalizedNormVec, triangleCenter, A_t] = compute_gradientGPU(vertex, face, scalar)
%%see https://blog.csdn.net/weixin_40211518/article/details/89469816
%     vertex = gpuArray(vertex);
%     face = gpuArray(face);
% 
% 
%     node1 = vertex(face(:,1), :);
%     node2 = vertex(face(:,2), :);
%     node3 = vertex(face(:,3), :);
%     A_t = vecnorm(cross(node2- node1, node3- node1, 2),2,2);
%     normalizedNormVec = cross(node2- node1, node3- node1, 2)./vecnorm(cross(node2- node1, node3- node1, 2),2,2);
%     gradientTangent = (scalar(face(:,2)) - scalar(face(:,1))).* (node1- node3)' ./ (2*A_t') + (scalar(face(:,3)) - scalar(face(:,1))).* (node2- node1)' ./ (2*A_t');
%     normalizedNormVec = dlarray(gpuArray(normalizedNormVec'),'CB');
%     gradientVec = cross(normalizedNormVec, gradientTangent, 1);
% 
%     triangleCenter = (node1 + node2 + node3)' / 3;

    % ---- vertex/face in GPU ----
% vertex = gpuArray(vertex);
% face   = gpuArray(face);

node1 = vertex(face(:,1), :);
node2 = vertex(face(:,2), :);
node3 = vertex(face(:,3), :);

triangleCenter = (node1 + node2 + node3)' / 3;

% ---- compute normals ----
Nvec = cross(node2 - node1, node3 - node1, 2);   % N×3
A_t  = vecnorm(Nvec,2,2);                        % N×1
normalizedNormVec = Nvec ./ A_t;                 % N×3

% ---- compute tangent gradient (contains Fpred) ----
gradientTangent = ...
    (scalar(face(:,2)) - scalar(face(:,1))) .* (node1' - node3') ./ (2*A_t') + ...
    (scalar(face(:,3)) - scalar(face(:,1))) .* (node2' - node1') ./ (2*A_t');

% ---- convert both to dlarray + correct shape (3×N) ----
%normalizedNormVec = dlarray(normalizedNormVec','CB');  % 3×N
%gradientTangent   = dlarray(gradientTangent,'CB');    % 3×N

% ---- cross ----
gradientVec = cross(normalizedNormVec', gradientTangent, 1);  % 3×N

end
