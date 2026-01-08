clear; close all; clc;
do_plot = true;
s = 'Mesh: Structured, graded mesh on dynamic gas lock';

n = 3;
m = 2*n;
M = 2*m;
N = M;
d = 0.0;

xa = 0; ya = 0;
xb = 1; yb = 0;
xc = 1; yc = 2;
xd = -1; yd = 2;


Q = [xa        ya;     % 1
     xb        yb;     % 2
     xc        yc;     % 3
     xd        yd];    % 4
     

% plot contour of geometry
% plot(P(:,1), P(:,2), 'Color', 'b', 'LineWidth', 2);
% hold on;

% mapping
mapx = @(x) pchip([0.0 0.5 1.0], [0.0 0.5 1.0], x);
mapy = @(y) pchip([0.0 0.5 1.0], [0.0 0.7 1.0], y);

[T, x, y] = triangulatepatch([xa xd; xb xc], [ya yd; yb yc], mapx(linspace(0,1,M+1)'), mapy(linspace(0,1,N+1)'), false);
% trimesh(T1, x4, y4, 'Color', 'k', 'LineWidth', 1);

%% combine mesh

% mesh coordinates;
P = [x y];
P = [P zeros(size(P,1),1)];

% mesh elements
n = [length(x)];

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

mesh.physobj(5).name = 'reservoir-left';
mesh.physobj(5).dim = 1;
mesh.physobj(5).tag = 5;

% add mesh entities
% mesh points
for k=1:4
    mesh.entity(k).dim = 0;
    mesh.entity(k).data = [Q(k,1) Q(k,2) 0];
end

% mesh lines
tol = 1e-4;

k = 5;
mesh.entity(k).dim = 1;
mesh.entity(k).data = extractboundary(mesh, Q(2,:)', Q(3,:)', tol);
mesh.entity(k).physobj = 'symmetry';
k = k + 1;

mesh.entity(k).dim = 1;
mesh.entity(k).data = extractboundary(mesh, Q(3,:)', Q(4,:)', tol);
mesh.entity(k).physobj = 'reservoir-top';
k = k + 1;

mesh.entity(k).dim = 1;
mesh.entity(k).data = extractboundary(mesh, Q(1,:)', Q(2,:)', tol);
mesh.entity(k).physobj = 'reservoir-bottom';
k = k + 1;

mesh.entity(k).dim = 1;
mesh.entity(k).data = extractboundary(mesh, Q(4,:)', Q(1,:)', tol);
mesh.entity(k).physobj = 'reservoir-left';
k = k + 1;

% mesh surfaces
mesh.entity(k).dim = 2;
mesh.entity(k).data = T;
mesh.entity(k).bctag = 1:2;
mesh.entity(k).physobj = 'volume';

write_msh(mesh, 'squaremesh-fitted-corners.msh');

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