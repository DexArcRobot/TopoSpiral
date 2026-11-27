function energyK = getEkFast(v, f, scalar, Bp, Tp, FaceSFMMatrix, boundPtIdx)
%GETEK2 Vectorized equivalent of getEk for a fixed, numeric surface mesh.
%   Inputs and energy/boundary conventions match getEk. SCALAR may be a
%   row/column vector or an unformatted dlarray. Geometry is not differentiated.
%   No vertex/face loops or vertex-normal calculations are needed.
%   As in getEk, degenerate faces and zero scalar gradients are not regularized.

scalar = reshape(scalar, 1, []);
[gradientVec, gradientTangent, ~, ~, A_t] = ...
    compute_gradientGPU(v, f, scalar);
normGradientTangent = gradientTangent ./ sqrt(sum(gradientTangent.^2, 1));
normGradientVec = gradientVec ./ sqrt(sum(gradientVec.^2, 1));

% Normal curvature in the level-set tangent direction.
uComponent = sum(Tp.' .* normGradientTangent, 1);
vComponent = sum(Bp.' .* normGradientTangent, 1);
Kn = uComponent.^2 .* FaceSFMMatrix(:,1).' + ...
    2*uComponent.*vComponent .* FaceSFMMatrix(:,2).' + ...
    vComponent.^2 .* FaceSFMMatrix(:,3).';
energyKn = sum(Kn.^2 .* A_t.');

% Cyclic edges: e0 = c-b, e1 = a-c, e2 = b-a.
a = v(f(:,1),:);
b = v(f(:,2),:);
c = v(f(:,3),:);
e0 = c-b;
e1 = a-c;
e2 = b-a;

% One sparse scatter matrix sums the three corner blocks at their vertices.
nFaces = size(f,1);
nVertices = size(v,1);
cornerToVertex = sparse(f(:), (1:3*nFaces).', ones(3*nFaces,1), ...
    nVertices, 3*nFaces);
Acorner = cornerAreas(e0, e1, e2);
Avertex = cornerToVertex * Acorner(:);

% Cotangents at a, b, c. Preserve the original denominator's eps.
denom = 2*A_t + eps;
cotA = -sum(e1.*e2,2) ./ denom;
cotB = -sum(e2.*e0,2) ./ denom;
cotC = -sum(e0.*e1,2) ./ denom;

% Original per-corner contributions, evaluated for every face at once.
X = normGradientVec.';
contributionA = sum((-cotB.*e1 + cotC.*e2).*X, 2);
contributionB = sum((-cotC.*e2 + cotA.*e0).*X, 2);
contributionC = sum((-cotA.*e0 + cotB.*e1).*X, 2);
Kg = cornerToVertex * [contributionA; contributionB; contributionC];
Kg = -Kg ./ (2*Avertex);
Kg(boundPtIdx) = 0;
energyKg = sum(Kg.^2 .* Avertex);
energyK = energyKn + energyKg;
end

function Acorner = cornerAreas(e0, e1, e2)
% Vectorize the exact area branches in CalcVertexNormals, including its
% nonacute-triangle convention (do not substitute barycentric areas).
de0 = sqrt(sum(e0.^2,2));
de1 = sqrt(sum(e1.^2,2));
de2 = sqrt(sum(e2.^2,2));
l2 = [de0.^2, de1.^2, de2.^2];
ew = l2 .* (l2(:,[2 3 1]) + l2(:,[3 1 2]) - l2);
s = (de0+de1+de2)/2;
Af = sqrt(s.*(s-de0).*(s-de1).*(s-de2));

% Preserve the if/elseif priority of the original implementation.
atA = ew(:,1) <= 0;
atB = ~atA & ew(:,2) <= 0;
atC = ~atA & ~atB & ew(:,3) <= 0;
acute = ~(atA | atB | atC);
Acorner = zeros(size(l2), 'like', l2);

ewscale = 0.5*Af(acute) ./ sum(ew(acute,:),2);
Acorner(acute,:) = ewscale .* ...
    (ew(acute,[2 3 1]) + ew(acute,[3 1 2]));

dot02 = sum(e0.*e2,2);
dot01 = sum(e0.*e1,2);
dot12 = sum(e1.*e2,2);
Acorner(atA,2) = -0.25*l2(atA,3).*Af(atA)./dot02(atA);
Acorner(atA,3) = -0.25*l2(atA,2).*Af(atA)./dot01(atA);
Acorner(atA,1) = Af(atA)-Acorner(atA,2)-Acorner(atA,3);

Acorner(atB,3) = -0.25*l2(atB,1).*Af(atB)./dot01(atB);
Acorner(atB,1) = -0.25*l2(atB,3).*Af(atB)./dot12(atB);
Acorner(atB,2) = Af(atB)-Acorner(atB,1)-Acorner(atB,3);

Acorner(atC,1) = -0.25*l2(atC,2).*Af(atC)./dot12(atC);
Acorner(atC,2) = -0.25*l2(atC,1).*Af(atC)./dot02(atC);
Acorner(atC,3) = Af(atC)-Acorner(atC,1)-Acorner(atC,2);
end
