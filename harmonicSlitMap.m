function [z, info] = harmonicSlitMap(V,F,innerBoundary,outerBoundary)
%HARMONICSLITMAP Approximate circular-slit conformal map using harmonic forms.
% [z,info] = harmonicSlitMap(V,F,innerBoundary,outerBoundary)
% V: N-by-2 or N-by-3 real coordinates; F: M-by-3 consistently oriented faces.
% Surface must be connected, orientable, manifold, genus zero, with >=2
% boundary loops. No unused vertices, duplicate faces, or degenerate faces.
% innerBoundary / outerBoundary: 1-by-1 positive integer boundary LOOP IDs.
% Loops are numbered from 1 in ascending order of their minimum vertex index.
% IDs must be distinct and <= the number of boundary loops.
% info.boundaryLoops{k} contains the ordered vertices of boundary k.
% Other boundaries become concentric arc slits. Orientation determines the
% complex orientation; planar faces should be counterclockwise.
% z: N-by-1 complex coordinates, outer radius 1; no explicit mesh cut needed.
% info.theta is a spanning-tree lift, NOT a continuous global log coordinate.
% FEM harmonic measures + Whitney-form constrained Hodge projection.
% Finite meshes only approximate conformality; inspect info diagnostics.
% Requires MATLAB with graph/conncomp (R2015b+). No extra toolbox.
% Call harmonicSlitMap() for a self-contained annulus demo.
if nargin==0
    nt=96; nr=16; r=exp(linspace(log(0.3),0,nr+1));
    t=(0:nt-1)'*2*pi/nt;
    V=[reshape(cos(t)*r,[],1),reshape(sin(t)*r,[],1)];
    F=zeros(2*nt*nr,3); k=0;
    for j=1:nr
        for i=1:nt
            a=(j-1)*nt+i; b=(j-1)*nt+mod(i,nt)+1;
            c=b+nt; d=a+nt;
            F(k+1,:)=[a d c]; F(k+2,:)=[a c b]; k=k+2;
        end
    end
    [z,info]=harmonicSlitMap(V,F,1,2);
    figure; subplot(1,2,1); triplot(F,V(:,1),V(:,2)); axis equal;
    title('Input annulus'); subplot(1,2,2);
    triplot(F,real(z),imag(z)); axis equal; title('Circular-slit map');
    fprintf('Inner radius: %.8f (analytic: 0.3)\n',info.radii(info.innerLoop));
    fprintf('Closure %.3e, period %.3e, CR error %.3e\n', ...
        info.closedResidual,info.periodResidual,info.relativeCRError);
    return
end
narginchk(4,4);
validateattributes(V,{'numeric'},{'2d','real','finite','nonempty'});
validateattributes(F,{'numeric'},{'2d','real','finite','integer','positive'});
assert(any(size(V,2)==[2 3]) && size(F,2)==3,'Invalid V or F dimensions.');
V=double(V); F=double(F); n=size(V,1); nf=size(F,1);
assert(~isempty(F) && max(F(:))<=n,'Invalid face indices.');
assert(numel(unique(F(:)))==n,'Remove unused vertices.');
assert(size(unique(sort(F,2),'rows'),1)==nf,'Duplicate faces.');
if size(V,2)==2, V(:,3)=0; end
% Global edges use increasing vertex indices; local edges follow each face.
raw=[F(:,[1 2]);F(:,[2 3]);F(:,[3 1])];
[E,~,eid]=unique(sort(raw,2),'rows'); ne=size(E,1);
sgn=2*(raw(:,2)>raw(:,1))-1;
count=accumarray(eid,1,[ne 1]);
assert(all(count<=2),'Nonmanifold edge.');
os=accumarray(eid,sgn,[ne 1]);
assert(all(os(count==2)==0),'Faces are not consistently oriented.');
assert(max(conncomp(graph(E(:,1),E(:,2),[],n)))==1,'Mesh is disconnected.');
isBoundary=count(eid)==1; be=raw(isBoundary,:);
assert(~isempty(be),'Surface has no boundary.');
bv=unique(be(:));
assert(all(accumarray(be(:,1),1,[n 1])==ismember((1:n)',bv)) && ...
       all(accumarray(be(:,2),1,[n 1])==ismember((1:n)',bv)), ...
       'Boundary is not a disjoint union of manifold loops.');
next=zeros(n,1); next(be(:,1))=be(:,2);
seen=false(n,1); loops={};
for s=bv'
    if seen(s), continue; end
    q=s; loop=[];
    while ~seen(q)
        seen(q)=true; loop(end+1,1)=q; q=next(q); %#ok<AGROW>
    end
    assert(q==s,'Invalid boundary loop.');
    loops{end+1}=loop; %#ok<AGROW>
end
nb=numel(loops);
assert(nb>=2 && n-ne+nf==2-nb,'Requires genus zero and >=2 boundaries.');
validateattributes(innerBoundary,{'numeric'}, ...
    {'real','finite','scalar','integer','>=',1,'<=',nb},mfilename,'innerBoundary',3);
validateattributes(outerBoundary,{'numeric'}, ...
    {'real','finite','scalar','integer','>=',1,'<=',nb},mfilename,'outerBoundary',4);
inner=double(innerBoundary); outer=double(outerBoundary);
assert(inner~=outer,'Inner and outer boundary IDs must be distinct.');
outerSeed=loops{outer}(1); % Internal integration anchor, not an input.
holes=setdiff(1:nb,outer,'stable'); nh=nb-1;
D0=sparse([(1:ne)';(1:ne)'],[E(:,1);E(:,2)], ...
    [-ones(ne,1);ones(ne,1)],ne,n);
D1=sparse(repmat((1:nf)',3,1),eid,sgn,nf,ne);
le=reshape(eid,nf,3); ls=reshape(sgn,nf,3);
% Face gradients and stiffness matrix; Whitney one-form mass matrix.
G=zeros(3,3,nf); normals=zeros(nf,3); areas=zeros(nf,1);
ki=zeros(9*nf,1); kj=ki; kv=ki; mi=ki; mj=ki; mv=ki;
lp=[1 2;2 3;3 1];
for f=1:nf
    ids=F(f,:); X=V(ids,:); nn=cross(X(2,:)-X(1,:),X(3,:)-X(1,:));
    a2=norm(nn); scale=max(sum((X-X([2 3 1],:)).^2,2));
    assert(a2>1e-14*scale && a2>0,'Degenerate or extremely thin face.');
    nn=nn/a2; A=a2/2; normals(f,:)=nn; areas(f)=A;
    g=[cross(nn,X(3,:)-X(2,:));cross(nn,X(1,:)-X(3,:)); ...
       cross(nn,X(2,:)-X(1,:))]/a2;
    G(:,:,f)=g; Kf=A*(g*g'); Mf=zeros(3);
    Q=A/12*(ones(3)+eye(3)); % integral lambda_i lambda_j
    for a=1:3
        i=lp(a,1); j=lp(a,2);
        for b=1:3
            k=lp(b,1); l=lp(b,2);
            Mf(a,b)=Q(i,k)*dot(g(j,:),g(l,:))-Q(i,l)*dot(g(j,:),g(k,:)) ...
                -Q(j,k)*dot(g(i,:),g(l,:))+Q(j,l)*dot(g(i,:),g(k,:));
        end
    end
    Mf=Mf.*(ls(f,:)'*ls(f,:));
    ix=(f-1)*9+(1:9); [ii,jj]=ndgrid(ids,ids);
    ki(ix)=ii(:); kj(ix)=jj(:); kv(ix)=Kf(:);
    [ii,jj]=ndgrid(le(f,:),le(f,:));
    mi(ix)=ii(:); mj(ix)=jj(:); mv(ix)=Mf(:);
end
K=sparse(ki,kj,kv,n,n); M=sparse(mi,mj,mv,ne,ne);
% Harmonic measures h_j=1 on hole j, zero on every other boundary.
H=zeros(n,nh);
for j=1:nh, H(loops{holes(j)},j)=1; end
interior=setdiff((1:n)',bv);
if ~isempty(interior)
    H(interior,:)=-K(interior,interior)\(K(interior,bv)*H(bv,:));
end
C=H'*K*H; C=(C+C')/2;
assert(rcond(C)>1e-12,'Ill-conditioned period matrix; improve the mesh.');
% Induced boundary orientation: outer CCW and holes CW in the plane.
% Integral_boundary *du = outward flux. The chosen hole has period -2*pi.
p=zeros(nh,1); p(holes==inner)=-2*pi;
c=C\p; u=H*c;
% Project the pointwise Hodge rotation onto CLOSED Whitney forms with the
% prescribed periods. Merely averaging *du on edges does not ensure closure.
B=sparse(nh,ne);
for j=1:nh
    loop=loops{holes(j)}; pairs=[loop loop([2:end 1])];
    [ok,idx]=ismember(sort(pairs,2),E,'rows'); assert(all(ok));
    B(j,:)=sparse(ones(numel(idx),1),idx, ...
        2*(pairs(:,2)>pairs(:,1))-1,1,ne);
end
rhs=zeros(ne,1); targetEnergy=0;
for f=1:nf
    g=G(:,:,f); grad=u(F(f,:))'*g;
    rot=cross(normals(f,:),grad); A=areas(f);
    targetEnergy=targetEnergy+A*dot(rot,rot);
    for a=1:3
        i=lp(a,1); j=lp(a,2);
        rhs(le(f,a))=rhs(le(f,a))+ls(f,a)*A/3*dot(g(j,:)-g(i,:),rot);
    end
end
Q=[D1;B]; target=[zeros(nf,1);p]; nq=size(Q,1);
sol=[M Q';Q sparse(nq,nq)]\[rhs;target]; beta=sol(1:ne);
% Integrate on a spanning tree. Closed forms with 2*pi-integer periods give
% path-independent exp(i*theta), even though theta itself is multivalued.
adj=sparse([E(:,1);E(:,2)],[E(:,2);E(:,1)], ...
    [(1:ne)';-(1:ne)'],n,n);
theta=zeros(n,1); visited=false(n,1); visited(outerSeed)=true;
queue=zeros(n,1); queue(1)=outerSeed; head=1; tail=1;
while head<=tail
    a=queue(head); head=head+1; [~,js,es]=find(adj(a,:));
    for k=1:numel(js)
        b=js(k);
        if ~visited(b)
            theta(b)=theta(a)+sign(es(k))*beta(abs(es(k)));
            visited(b)=true; tail=tail+1; queue(tail)=b;
        end
    end
end
z=exp(u+1i*theta);
info.boundaryLoops=loops; info.outerLoop=outer; info.innerLoop=inner;
info.radii=cellfun(@(x) exp(mean(u(x))),loops);
info.harmonicMeasures=H; info.coefficients=c; info.periodMatrix=C;
info.edges=E; info.exactForm=D0*u; info.conjugateForm=beta;
info.u=u; info.theta=theta; info.periods=B*beta;
info.closedResidual=norm(D1*beta,inf);
info.periodResidual=norm(B*beta-p,inf);
info.phaseResidual=norm(angle(exp(1i*(D0*theta-beta))),inf);
info.relativeCRError=sqrt(max(0,beta'*M*beta-2*beta'*rhs+targetEnergy) ...
    /max(targetEnergy,eps));
w1=z(F(:,2))-z(F(:,1)); w2=z(F(:,3))-z(F(:,1));
signedArea=imag(conj(w1).*w2)/2;
info.flippedFaces=find(signedArea<0);
info.degenerateFaces=find(abs(signedArea)<=1e-12*max(abs(signedArea)));
if ~isempty(info.flippedFaces)
    warning('harmonicSlitMap:Flips','%d flipped triangles; refine/improve mesh.', ...
        numel(info.flippedFaces));
end
end
