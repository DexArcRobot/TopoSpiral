function[E, energyW, energyK, lossNeuman] = getScalarFieldEnergy(scalar, v, f, Bp, Tp, FaceSFMMatrix, k_c, EsMode, loss_alpha, boundPtIdx)
    energyW = 0;    
    energyK = 0;    

    [curvature2Gradient, gradientVec, triangleCenter, A_t] = getCurvature2Gradient(v, f, scalar', Bp, Tp, FaceSFMMatrix);
    
    if loss_alpha ~= 1
        energyW = calBasicEnergy(gradientVec, k_c, curvature2Gradient, A_t, EsMode(1));
    end
    if loss_alpha ~= 0
        %energyK = getEk(v, f, scalar', Bp, Tp, FaceSFMMatrix, boundPtIdx);
        energyK = getEkFast(v, f, scalar', Bp, Tp, FaceSFMMatrix, boundPtIdx);
    end

     lossNeuman = 0;
    % if EsMode(3) ~= 0 %矩形场景
    %     eLen = v(boundPtIdx, :) - v(circshift(boundPtIdx,1), :);
    %     eLen = eLen(bdNum, :);
    %     [gradientVec, ~, ~, ~, ~] = compute_gradient(v, f, scalar');
    %     %gradientTangent = gradientTangent';
    %     %gradientTangent(:, bdFIdx) = gradientTangent(:, bdFIdx)./vecnorm(gradientTangent(:, bdFIdx),2,1);
    %     gradientVec(bdFIdx, :) = gradientVec(bdFIdx, :)./vecnorm(gradientVec(bdFIdx, :),2,2);
    %     %gradientVec(bdFIdx, :) = gradientVec(bdFIdx, :)./vecnorm(gradientVec(bdFIdx, :),2,1);
    % 
    %     lossNeuman = cross(gradientVec(bdFIdx, :), eLen, 2);
    %     lossNeuman = sum(vecnorm(lossNeuman, 2, 2));
    %     %gradientTangent = gradientTangent';
    %     %triangleCenter = triangleCenter';
    %     %quiver3(triangleCenter(bdFIdx, 1),triangleCenter(bdFIdx, 2),triangleCenter(bdFIdx, 3),gradientVec(bdFIdx, 1),gradientVec(bdFIdx, 2),gradientVec(bdFIdx, 3));
    %     %quiver3(triangleCenter(1,bdFIdx),triangleCenter(2,bdFIdx),triangleCenter(3,bdFIdx),gradientTangent(1, bdFIdx),gradientTangent(2, bdFIdx),gradientTangent(3, bdFIdx));
    % end
    E = (1-loss_alpha)*energyW + loss_alpha*energyK; %+ lossNeuman * EsMode(3);

end

