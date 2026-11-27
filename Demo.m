close all
clear
clc
global frame;

k_c = 1/3;
alpha = 0.001;
EsMode = [5, 2];%[5, 2]:symmetrized Itakura–Saito-type divergence. shen 2025 topology;  [1, 2]:zou qiang CAD 2014

[v, f, n, objname] = stlRead("Input\headlight.stl");
v = v(:,[1,3,2]);

[freeboundsList, longestBoundIdx] = getMeshBounds(f, v);


% inner_bd_id = 5;
% [tau] = annularSlitMapping(v, f, inner_bd_id);
% vCSM = exp(tau).';
% vFlatten = [real(vCSM)', imag(vCSM)'];


[vCSM, info] = harmonicSlitMap(v,f,4,1);
vFlatten = [real(vCSM), imag(vCSM)];


for idx = 1:length(freeboundsList)
    bdPt = freeboundsList{idx};
    radOnBd = vecnorm(vFlatten(bdPt, :), 2, 2);
    avgRadOnBd = mean(radOnBd);
    vFlatten(bdPt, :) = vFlatten(bdPt, :) ./ (radOnBd + eps) * avgRadOnBd;
end

vCSM = (vFlatten(:, 1) + 1i * vFlatten(:, 2)).';
%%
if 1 %headlight.stl
[pt,dist] = dsearchn([real(vCSM)', imag(vCSM)'],[0.796683,-0.326572]);
angleDeg = -0.1;
V_rot = rotatePoints([real(vCSM(pt))', imag(vCSM(pt))'], angleDeg)
vCSM(pt) =  V_rot(1) + 1i*V_rot(2);

[pt,dist] = dsearchn([real(vCSM)', imag(vCSM)'],[0.48716,0.631216]);
angleDeg = -0.01;
V_rot = rotatePoints([real(vCSM(pt))', imag(vCSM(pt))'], angleDeg)
vCSM(pt) =  V_rot(1) + 1i*V_rot(2);
end
if 0 %headlight1500.stl
[pt,dist] = dsearchn([real(vCSM)', imag(vCSM)'],[-0.350982,-0.862266]);
angleDeg = 0.5;
V_rot = rotatePoints([real(vCSM(pt))', imag(vCSM(pt))'], angleDeg)
vCSM(pt) =  V_rot(1) + 1i*V_rot(2);

[pt,dist] = dsearchn([real(vCSM)', imag(vCSM)'],[0.799955,-0.318322]);
angleDeg = -1;
V_rot = rotatePoints([real(vCSM(pt))', imag(vCSM(pt))'], angleDeg)
vCSM(pt) =  V_rot(1) + 1i*V_rot(2);
end
%%
TR1 = triangulation(f, real(vCSM)', imag(vCSM)'); %conformal slit mapping mesh
TR2 = triangulation(f, v(:,1), v(:,2), v(:,3)); %3D surface mesh
TR3 = triangulation(f, vFlatten(:,1), vFlatten(:,2)); %flattened surface mesh

figure(1)
hold on
axis equal
triplot(TR1,'color',[0.9 0.9 0.9])
[flipped, degenerate] = findFlippedTriangles([real(vCSM)', imag(vCSM)'], f);

if ~isempty(flipped)
    patch('Vertices', [real(vCSM)', imag(vCSM)'], 'Faces', f(flipped,:), ...
          'FaceColor', 'r', 'EdgeColor', 'r','lineWidth',3);
end

figure(2)
hold on
axis equal
trimesh(TR2)



figure(3)
hold on
axis equal
triplot(TR3)

ptOnBound = [];
for idx = 1:length(freeboundsList)
    ptOnBound = [ptOnBound; freeboundsList{idx}];
end
ptNotOnBound = setdiff((1:length(vCSM))',ptOnBound);

%%
FV.faces = f;
FV.vertices = v;
[FaceSFM, Bp, Tp] = GetSecondBasicForm( FV );
FaceSFMMatrix = zeros(length(FaceSFM),3);
for idx = 1:length(FaceSFM)
    curFaceSFM = cell2mat(FaceSFM(idx));
    FaceSFMMatrix(idx,:) = [curFaceSFM(1,1), curFaceSFM(1,2), curFaceSFM(2,2)];
end

%%

itStep = 0.1;

boundPtIdx = freeBoundary(TR1);
scalar = abs(vCSM);
loss = getScalarFieldEnergy(scalar, v, f,  ...
    Bp, Tp, FaceSFMMatrix, k_c, EsMode, alpha, boundPtIdx);
loss

%% ordered_ring{} 包含顶点的三角组成的网格边界上的点，逆时针排列

num_vertices = length(vCSM);
vertex_triangles = cell(num_vertices, 1);

for vId = 1:length(vertex_triangles)
    [m,n] = find(f == vId);
    vertex_triangles{vId} = unique(m);
end


ordered_ring = {};
for vIdx = 1:num_vertices
    f_temp = f(vertex_triangles{vIdx}, :);
    [ringBounds, ~] = getMeshBounds(f_temp, [real(vCSM)', imag(vCSM)']);
    vBoundId = find(ringBounds{1} == vIdx);
    if ~isempty(vBoundId)
        ringBounds{1} = circshift(ringBounds{1}, 1 - vBoundId);%在边上的中心点放队首
    end
    ordered_ring{vIdx} = ringBounds{1};
end



%%
EE= [];
iterTimes = 250;
colrr = cool(iterTimes);
close all
if ~exist('t_vCSM')
    t_vCSM = vCSM*30;
end
vCSM = t_vCSM;    

figure(1)
hold on
axis equal
triplot(TR1,'color',[0 0 0],'linewidth',3)
%vCSM([100 150  200 300]) = vCSM([100 150  200 300])*1.05;

scalar = abs(vCSM);%scalar = imag(vCSM); %debug
[~, ~, normalizedNormVec, ~, ~] = compute_gradient(v, f, scalar');
[~, Avertex, ~, ~, ~] = CalcVertexNormals(FV, normalizedNormVec);
allE = [];
lambda = [1 100 0]; % [Ew Ek Ec]

k_c = 1;
vertexEnergyGradient = zeros(length(vCSM),2);
adaptionGrad = vertexEnergyGradient;

bd1Ring = getBoundary1Ring(f, v);

for idx = 1:iterTimes

    if idx == 200
        itStep = 0.01;
    end

    x = real(vCSM)';
    y = imag(vCSM)';
    z = x*0;
% 
%     maxNorm = max(vecnorm(vertexEnergyGradient, 2, 2));
    figure(1)
    hold on
    axis equal
    TRM2  = triangulation(f,[x,y,z]);
    trimesh(TRM2,'EdgeColor',colrr(idx,:),'FaceColor','none');
    %plot(real(vCSM(freeboundsList{1}(cornerIdx{1}))),imag(vCSM(freeboundsList{1}(cornerIdx{1}))),'k*')

    %
    scalar = abs(vCSM);%scalar = imag(vCSM); %debug
    [curvature2Gradient, gradientVec, triangleCenter, A_t] = getCurvature2Gradient(v, f, scalar', Bp, Tp, FaceSFMMatrix);
    
    energyW = calBasicEnergy(gradientVec, k_c, curvature2Gradient, A_t, EsMode(1));
    energyK = getEk(v, f, scalar', Bp, Tp, FaceSFMMatrix, boundPtIdx);
    %
    vertexEnergy = 0;
    vertexEnergyGradient = 0;
    loss
    if 0
        [vertexEnergy, vertexEnergyGradient] = variationalDescent20250904(f, v, vCSM, lambda, FaceSFMMatrix, Tp, Bp, EsMode, k_c, vertex_triangles, ordered_ring, Avertex);
    else
        scalar = dlarray(scalar);
        [loss_t, gradients] = dlfeval(@modelGradients, scalar, ...
            v, f, Bp, Tp, FaceSFMMatrix, k_c, EsMode, alpha, boundPtIdx);
        
        % 对每个点半径的负偏导，统一为列向量
        negativeGradient = -extractdata(gradients(:));
        
        % 每个点自己的径向单位向量，N×2
        vCSM_XY = [real(vCSM(:)), imag(vCSM(:))];
        radialUnit = vCSM_XY ./ (vecnorm(vCSM_XY,2,2) + eps);
        
        % 转成二维下降方向，N×2
        vertexEnergyGradient = negativeGradient .* radialUnit;
        loss = extractdata(loss_t);
    end
    
    [dequeIntersecCell] = pointConvexIterArea(f, vCSM, ordered_ring);
    adaptionGrad = gradientAdaptation(vertexEnergyGradient, dequeIntersecCell, vCSM, freeboundsList, Avertex);
    %adaptionGrad = vertexEnergyGradient;

    kk = (1:length(vCSM))';
    %kk = [100 150  200 300];
%     itV = 0.05;
    % for idx2 = 5
    %     adaptionGrad(freeboundsList{idx2},:) = 0;
    %     vertexEnergyGradient(freeboundsList{idx2},:) = 0;
    % end

    if mod(idx,20) == 1

        figure()
        subplot(2,1,1)
        TR4 = triangulation(f, real(vCSM)', imag(vCSM)'); %conformal slit mapping mesh
        hold on
        axis equal
%         xlim([-max(abs(vCSM))*1.1,max(abs(vCSM))*1.1])
%         ylim([-max(abs(vCSM))*1.1,max(abs(vCSM))*1.1])
        triplot(TR4,'color',[0.3 0.9 0.3])
        for bdIdx = 1:length(freeboundsList)
            plot(vCSM([freeboundsList{bdIdx};freeboundsList{bdIdx}(1)]),'LineWidth',2);
        end
        quiver(real(vCSM)', imag(vCSM)',vertexEnergyGradient(:,1),vertexEnergyGradient(:,2),0,'b');
        quiver(real(vCSM)', imag(vCSM)',adaptionGrad(:,1),adaptionGrad(:,2),0,'r');
%         quiver([real(vCSM(freeboundsList{1}))';real(vCSM(freeboundsList{2}))';real(vCSM(freeboundsList{3}))';real(vCSM(freeboundsList{4}))'], ...
%        [imag(vCSM(freeboundsList{1}))';imag(vCSM(freeboundsList{2}))';imag(vCSM(freeboundsList{3}))';imag(vCSM(freeboundsList{4}))'], ...
%        [adaptionGrad((freeboundsList{1}),1);adaptionGrad((freeboundsList{2}),1);adaptionGrad((freeboundsList{3}),1);adaptionGrad((freeboundsList{4}),1)], ...
%        [adaptionGrad((freeboundsList{1}),2);adaptionGrad((freeboundsList{2}),2);adaptionGrad((freeboundsList{3}),2);adaptionGrad((freeboundsList{4}),2)],'r','AutoScale','off');
        if 0
            subplot(2,1,2)
            hold on
            axis equal
            triplot(TR2);
            myPlot(v(freeboundsList{1},:),'k-')
            r = abs(vCSM(:));
            rmin = min(r);
            rmax = max(r);
            rLevels = linspace(rmin, rmax, 101);
            rLevels = rLevels(2:end-1);
            [C, h] = tricontour(f, v(:,1), v(:,2), r, rLevels);
        else
            subplot(2,1,2)
            hold on
            axis equal
            triplot(TR2);
            myPlot(v(freeboundsList{1},:),'k-')
    %         xlim([-max(abs(vCSM))*1.1,max(abs(vCSM))*1.1])
    %         ylim([-max(abs(vCSM))*1.1,max(abs(vCSM))*1.1])
            for idx2 = 1:1:99
                r = min(abs(vCSM)) + (max(abs(vCSM))-min(abs(vCSM)))/100 * idx2;
                t = 0:pi/500:2*pi;
                circle = [cos(t)', sin(t)'] * r;
                w2D = inverseMapping(circle, TR4, TR2);
                hold on
                plot(w2D(:,1), w2D(:,2),'r');
            end
        end

    end
    %adaptionGrad(bd1Ring,:) = 0;
    adaptionGrad(freeboundsList{4},:) = 0;
    %var = 1e-2;
    %vertexEnergyGradient = vertexEnergyGradient./max(vecnorm(vertexEnergyGradient,2,2));
    vCSM(kk) = vCSM(kk) + itStep*(adaptionGrad(kk,1) + 1i*adaptionGrad(kk,2)).'*1^idx;
    %vCSM(kk) = vCSM(kk) + var*(vertexEnergyGradient(kk,1) + 1i*vertexEnergyGradient(kk,2)).' * idx;
    %vCSM = vCSM + 0.1*(adaptionGrad(:,1) + 1i*adaptionGrad(:,2)).';
    %vCSM(ptOnBound) = vCSM(ptOnBound) + 0.004*(adaptionGrad(ptOnBound,1) + 1i*adaptionGrad(ptOnBound,2)).'*1^idx;

% %    kk = (1:length(vCSM))';  
% %    vCSM(kk) = vCSM(kk) + 0.01/maxNorm*(vertexEnergyGradient(kk,1).' + 1i*vertexEnergyGradient(kk,2).');
     allE = [allE;sum(vertexEnergy), energyW, energyK];

    [sum(vertexEnergy), energyW, energyK]
    EE= [EE; energyW, energyK,loss];
    idx
    [flipped, degenerate] = findFlippedTriangles([real(vCSM)', imag(vCSM)'], f);

    if ~isempty(flipped)
        TRInvalid = triangulation(f, real(vCSM)', imag(vCSM)'); %conformal slit mapping mesh
        figure()
        hold on
        axis equal
        triplot(TRInvalid,'color',[0.9 0.9 0.9])
        patch('Vertices', [real(vCSM)', imag(vCSM)'], 'Faces', f(flipped,:), ...
              'FaceColor', 'r', 'EdgeColor', 'r','lineWidth',3);
        warning("invalid tri");
        break;
    end
end

figure()
plot(EE(:,3))

figure()
subplot(2,1,1)
TR4 = triangulation(f, real(vCSM)', imag(vCSM)'); %conformal slit mapping mesh
hold on
axis equal
        xlim([-max(abs(vCSM))*1.1,max(abs(vCSM))*1.1])
        ylim([-max(abs(vCSM))*1.1,max(abs(vCSM))*1.1])
triplot(TR4,'color',[0.3 0.9 0.3])
for bdIdx = 1:length(freeboundsList)
    plot(vCSM([freeboundsList{bdIdx};freeboundsList{bdIdx}(1)]),'LineWidth',2);
end
quiver(real(vCSM)', imag(vCSM)',vertexEnergyGradient(:,1),vertexEnergyGradient(:,2),0,'b');
quiver(real(vCSM)', imag(vCSM)',adaptionGrad(:,1),adaptionGrad(:,2),0,'r');

subplot(2,1,2)
hold on
axis equal
triplot(TR2);
myPlot(v(freeboundsList{1},:),'k-')
for idx2 = 1:1:99
    r = min(abs(vCSM)) + (max(abs(vCSM))-min(abs(vCSM)))/100 * idx2;
    t = 0:pi/5000:2*pi;
    circle = [cos(t)', sin(t)'] * r;
    w2D = inverseMapping(circle, TR4, TR2);
    hold on
    plot(w2D(:,1), w2D(:,2),'r');
end
%% ------------------------
%  自定义前向 + loss 函数
% ------------------------
function [loss, gradients] = modelGradients(scalar, ...
    v, f, Bp, Tp, FaceSFMMatrix, k_c, EsMode, alpha, boundPtIdx)

    %scalar = scalar(:);  % 使用线性插值

    loss = getScalarFieldEnergy(scalar, v, f,  ...
        Bp, Tp, FaceSFMMatrix, k_c, EsMode, alpha, boundPtIdx);
    gradients = dlgradient(loss, scalar);

    %loss = dlarray(loss);          
    % 自动求梯度
    isgpuarray(loss)
    % gradients = dlgradient(loss, net.Learnables);
    % if any(isnan(extractdata(loss)))
    %     disp("Loss is NaN");
    % end

end

