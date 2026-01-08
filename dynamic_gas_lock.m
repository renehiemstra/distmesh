clear; close all; clc;
do_plot = true;
s = 'Mesh: Structured, graded mesh on dynamic gas lock';

n = 8;
m = 2*n;
M = 2*m;
N = M;
d = 0.0;
split = false;


xa = 0.0; ya = 0.0;
xb = 0.9; yb = 0.45;
xc = 1.9; yc = 0.9;
xd = 1.7427; yd = -13.1;
xe = 6.3427; ye = 16.7;
xf = 10.3427;
xg = -3.45208;
delta = d * (ye - yd) / 2;

Q = [xa         ya;     % 1
     xb         ya;     % 2
     xb         yb;     % 3
     xc         yb;     % 4
     xe         yd;     % 10
     xf         yd;     % 9
     xf         ye;     % 8
     xg         ye;     % 7
     xd         yc;     % 5
     xa         yc;     % 6
     xa         ya];    % 1

% plot contour of geometry
% plot(P(:,1), P(:,2), 'Color', 'b', 'LineWidth', 2);
% hold on;

% mapping
mapx = @(x) pchip([0.0 0.5 1.0], [0.0 0.74 1.0], x);
mapy = @(y) pchip([0.0 1/3 2/3 1.0], [0.0 0.2 0.8 1.0], y);
mapyd = @(y) pchip([0.0 1/3 2/3 1.0], [0.0 0.75 0.94 1.0], y);
mapyu = @(y) pchip([0.0 1/3 2/3 1.0], [0.0 0.06 0.25 1.0], y);

% create triangulation
[T1, x1, y1] = triangulatepatch([xa xa; xb xb], [ya yb; ya yb], linspace(0,1,m+1)', linspace(0,1,n+1)', split);
[T2, x2, y2] = triangulatepatch([xa xa; xb xb], [yb yc; yb yc], linspace(0,1,m+1)', mapy(linspace(0,1,n+1)'), split);
[T3, x3, y3] = triangulatepatch([xb xb; xc xd], [yb yc; yb yc], linspace(0,1,m+1)', mapy(linspace(0,1,n+1)'), split);

% trimesh(T1, x1, y1, 'Color', 'k', 'LineWidth', 1);
% trimesh(T2, x2, y2, 'Color', 'k', 'LineWidth', 1);
% trimesh(T3, x3, y3, 'Color', 'k', 'LineWidth', 1);
% axis equal;

[T4, x4, y4] = triangulatepatch([xe xc; xf xf], [yd yb; yd yb-delta], mapx(linspace(0,1,M+1)'), mapyu(linspace(0,1,N+1)'), split);
% trimesh(T4, x4, y4, 'Color', 'k', 'LineWidth', 1);

[T5, x5, y5] = triangulatepatch([xc xd; xf xf], [yb yc; yb-delta yc+delta], mapx(linspace(0,1,M+1)'), mapy(linspace(0,1,n+1)'), split);
% trimesh(T5, x5, y5, 'Color', 'k', 'LineWidth', 1);

[T6, x6, y6] = triangulatepatch([xd xg; xf xf], [yc ye; yc+delta ye], mapx(linspace(0,1,M+1)'), mapyd(linspace(0,1,N+1)'), split);
% trimesh(T6, x6, y6, 'Color', 'k', 'LineWidth', 1);

%% combine mesh

% mesh coordinates
P = [x1 y1;
     x2 y2;
     x3 y3;
     x4 y4
     x5 y5
     x6 y6];
P = [P zeros(size(P,1),1)];

% mesh elements
n = [length(x1); length(x2); length(x3); length(x4); length(x5); length(x6);];
T = [T1; T2+n(1); T3+sum(n(1:2)); T4+sum(n(1:3)); T5+sum(n(1:4)); T6+sum(n(1:5))];

% mesh elements
[P,~,ic] = unique(P, 'rows');
T = ic(T(:,:));

% export mesh
mesh.vertices = P;
mesh.faces = T;

% plot mesh
if( do_plot )
  clf
  c = [.9, .9, .9];
  patch( 'vertices', mesh.vertices, 'faces', mesh.faces, 'facecolor', c );
  title(s);
  axis tight;
  axis equal;
  %axis off;
end

%% connectivity arrays

% number of elements
npts  = size(P, 1);
nelms = size(T, 1);

% halfedges
he = zeros(nelms*3, 2);
e  = zeros(ceil(1.2*ceil(nelms*3 / 2)), 2);

% check normals
for k=1:nelms
    a = P(T(k,1), :);
    b = P(T(k,2), :);
    c = P(T(k,3), :);
    n = cross(b-a, c-a);
    assert(n(1)==0 & n(2)==0 & n(3)>0);
end

% assign halfedges
for k=1:nelms
    i = 3*(k-1);
    he(i+1, 1) = T(k,1); he(i+1, 2) = T(k,2); % halfedge 1
    he(i+2, 1) = T(k,2); he(i+2, 2) = T(k,3); % halfedge 2
    he(i+3, 1) = T(k,3); he(i+3, 2) = T(k,1); % halfedge 3
end

% active list of halfedges
active = ones(nelms*3, 1);

ecounter = 1;
while true
    this = find(active==1, 1, 'first');
    %fprintf('%d\n', this);
    if isscalar(this)==false
        break % exit of list of 'active' halfedges is empty
    end
    other = find(he(:,2) == he(this, 1) & he(:,1) == he(this, 2));
    e(ecounter, 1) = this;
    active(this) = 0;
    if isscalar(other)
        e(ecounter, 2) = other;
        active(other) = 0;
    end
    ecounter = ecounter + 1;
end

e = e(1:ecounter-1, :);

% extract boundary edges
be = find(e(:,2) == 0);

bpts = he(e(be,1),:);
[bpts, ia, ic] = unique(bpts(:));

% plot the boundary points
hold on;
plot(mesh.vertices(bpts,1), mesh.vertices(bpts,2), 'o');

%% compute incidense matrix between edges and nodes
mesh.vertices = P;
mesh.faces = T;
mesh.edges = e;
mesh.halfedges = he;
mesh.incidence = incidensemat(mesh);
mesh.boundary = bpts;

% add physical objects
mesh.physobj(1).name = 'volume';
mesh.physobj(1).dim = 2;
mesh.physobj(1).tag = 1;

mesh.physobj(2).name = 'symmetry';
mesh.physobj(2).dim = 1;
mesh.physobj(2).tag = 2;

mesh.physobj(3).name = 'reservoir-top';
mesh.physobj(3).dim = 1;
mesh.physobj(3).tag = 3;

mesh.physobj(4).name = 'reservoir-bottom';
mesh.physobj(4).dim = 1;
mesh.physobj(4).tag = 4;

mesh.physobj(5).name = 'reservoir-left-top';
mesh.physobj(5).dim = 1;
mesh.physobj(5).tag = 5;

mesh.physobj(6).name = 'reservoir-left-bottom';
mesh.physobj(6).dim = 1;
mesh.physobj(6).tag = 6;

mesh.physobj(7).name = 'channel-top';
mesh.physobj(7).dim = 1;
mesh.physobj(7).tag = 7;

mesh.physobj(8).name = 'channel-bottom';
mesh.physobj(8).dim = 1;
mesh.physobj(8).tag = 8;

mesh.physobj(9).name = 'small-reservoir-left';
mesh.physobj(9).dim = 1;
mesh.physobj(9).tag = 9;

mesh.physobj(10).name = 'small-reservoir-bottom';
mesh.physobj(10).dim = 1;
mesh.physobj(10).tag = 10;

mesh.physobj(11).name = 'small-reservoir-right';
mesh.physobj(11).dim = 1;
mesh.physobj(11).tag = 11;

% add mesh entities
% mesh points
for k=1:10
    mesh.entity(k).dim = 0;
    mesh.entity(k).data = [Q(k,1) Q(k,2) 0];
end

% mesh lines
tol = 1e-4;

mesh.entity(11).dim = 1;
mesh.entity(11).data = extractboundary(mesh, Q(6,:)', Q(7,:)', tol);
mesh.entity(11).physobj = 'symmetry';

mesh.entity(12).dim = 1;
mesh.entity(12).data = extractboundary(mesh, Q(7,:)', Q(8,:)', tol);
mesh.entity(12).physobj = 'reservoir-top';

mesh.entity(13).dim = 1;
mesh.entity(13).data = extractboundary(mesh, Q(5,:)', Q(6,:)', tol);
mesh.entity(13).physobj = 'reservoir-bottom';

mesh.entity(14).dim = 1;
mesh.entity(14).data = extractboundary(mesh, Q(8,:)', Q(9,:)', tol);
mesh.entity(14).physobj = 'reservoir-left-top';

mesh.entity(15).dim = 1;
mesh.entity(15).data = extractboundary(mesh, Q(4,:)', Q(5,:)', tol);
mesh.entity(15).physobj = 'reservoir-left-bottom';

mesh.entity(16).dim = 1;
mesh.entity(16).data = extractboundary(mesh, Q(9,:)', Q(10,:)', tol);
mesh.entity(16).physobj = 'channel-top';

mesh.entity(17).dim = 1;
mesh.entity(17).data = extractboundary(mesh, Q(3,:)', Q(4,:)', tol);
mesh.entity(17).physobj = 'channel-bottom';

mesh.entity(18).dim = 1;
mesh.entity(18).data = extractboundary(mesh, Q(1,:)', Q(2,:)', tol);
mesh.entity(18).physobj = 'small-reservoir-bottom';

mesh.entity(19).dim = 1;
mesh.entity(19).data = extractboundary(mesh, Q(2,:)', Q(3,:)', tol);
mesh.entity(19).physobj = 'small-reservoir-right';

mesh.entity(20).dim = 1;
mesh.entity(20).data = extractboundary(mesh, Q(10,:)', Q(1,:)', tol);
mesh.entity(20).physobj = 'small-reservoir-left';


% mesh surfaces
k = 21;
mesh.entity(k).dim = 2;
mesh.entity(k).data = T;
mesh.entity(k).bctag = 1:6;
mesh.entity(k).physobj = 'volume';

write_msh(mesh, 'dynamic-gas-lock-boundary-layer-mesh-fine.msh');

%% plot boundaries
colors = [0 0.4470 0.7410; 
          0.8500 0.3250 0.0980;
          0.9290 0.6940 0.1250; 
          0.4940 0.1840 0.5560;
          0.4660 0.6740 0.1880;
          0.3010 0.7450 0.9330;
          0.6350 0.0780 0.1840;
          0 0 1
          0 1 0;
          1 0 0];
k = 1;
for i=1:length(mesh.entity)
    if mesh.entity(i).dim ==1
        v = mesh.entity(i).data;
        plot(P(v(:), 1), P(v(:), 2), 'Color', colors(k,:), 'LineWidth',4);
        k = k + 1;
    end
end
