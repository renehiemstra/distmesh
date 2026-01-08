clear; close all; clc;
do_plot = true;
s = 'Mesh: Structured, graded mesh on dynamic gas lock';

n = 2;
m = 4;
n1 = 2*n; n2 = n; n3 = 8*n; 
m1 = m; m2 = 3*m; m3 = 10*m; m4 = 5*m; m5 = 6*m;
M = 2*m;
N = M;
d = 0.0;
alpha = 0.8;
beta = 0.85;
gamma = 0.9;
split = false;

xa =   0.000;   ya =  0.0000;
xb =   1.253;   yb =  0.4411;
xc =   3.840;   yc =  1.9430;
xd =  13.840;   yd =  4.1120;
xe = 113.750;   ye = 43.1400;
xf = (1-alpha)* xe + alpha * xd;

yf = (1-alpha)* ya + alpha * yc;
yg = (1-alpha)* ya + alpha * yb;
yh = (1-alpha)* ya + alpha * yd;
yi = (1-beta)* ye + beta * ya;
yk = (1-gamma)* ye + gamma * ya;

Q = [xa         ya;     % 1
     xb         ya;     % 2
     xc         ya;     % 3
     xd         ya;     % 4
     xf         ya;     % 5
     xe         ya;     % 6
     xa         yf;     % 7
     xb         yf;     % 8
     xc         yg;     % 9
     xd         yh;     % 10
     xf         yk;     % 11
     xe         yk;     % 12
     xa         yc;     % 13         
     xb         yc;     % 14
     xc         yb;     % 15
     xd         yd;     % 16
     xf         yi;     % 17    
     xe         yi;     % 18
     xd         ye;     % 19
     xf         ye;     % 20
     xe         ye;];    % 21
        
% plot contour of geometry
% plot(P(:,1), P(:,2), 'Color', 'b', 'LineWidth', 2);
% hold on;

% mapping
mapx2 = @(x) pchip([0.0 1/3 2/3 1.0], [0.0 0.2 0.5 1.0], x);
mapx3 = @(x) pchip([0.0 1/3 2/3 1.0], [0.0 0.4 0.75 1.0], x);
mapx4 = @(x) pchip([0.0 1/3 2/3 1.0], [0.0 0.6 0.85 1.0], x);
mapx5 = @(x) pchip([0.0 1/3 2/3 1.0], [0.0 0.4 0.725 1.0], x);


xx = [0.0 1/3 2/3 1.0]; yy = [0.0 0.6 0.85 1.0];
mapy = @(y) pchip(xx, yy, y);
mapy2 = @(y) pchip(xx, 1-fliplr(yy), y);

% create triangulation
[TT{1},  x{1},  y{1}] = triangulatepatch([Q(1,1) Q(7,1); Q(2,1) Q(8,1)], [Q(1,2) Q(7,2); Q(2,2) Q(8,2)], linspace(0,1,m1+1)', linspace(0,1,n1+1)', split);
[TT{2},  x{2},  y{2}] = triangulatepatch([Q(2,1) Q(8,1); Q(3,1) Q(9,1)], [Q(2,2) Q(8,2); Q(3,2) Q(9,2)], mapx2(linspace(0,1,m2+1))', linspace(0,1,n1+1)', split);
[TT{3},  x{3},  y{3}] = triangulatepatch([Q(3,1) Q(9,1); Q(4,1) Q(10,1)], [Q(3,2) Q(9,2); Q(4,2) Q(10,2)], mapx3(linspace(0,1,m3+1))', linspace(0,1,n1+1)', split);
[TT{4},  x{4},  y{4}] = triangulatepatch([Q(4,1) Q(10,1); Q(5,1) Q(11,1)], [Q(4,2) Q(10,2); Q(5,2) Q(11,2)], mapx4(linspace(0,1,m4+1))', linspace(0,1,n1+1)', split);
[TT{5},  x{5},  y{5}] = triangulatepatch([Q(5,1) Q(11,1); Q(6,1) Q(12,1)], [Q(5,2) Q(11,2); Q(6,2) Q(12,2)], mapx5(linspace(0,1,m5+1))', linspace(0,1,n1+1)', split);


[TT{6},  x{6},  y{6}] = triangulatepatch([Q(7,1) Q(13,1); Q(8,1) Q(14,1)], [Q(7,2) Q(13,2); Q(8,2) Q(14,2)], linspace(0,1,m1+1)', linspace(0,1,n2+1)', split);
[TT{7},  x{7},  y{7}] = triangulatepatch([Q(8,1) Q(14,1); Q(9,1) Q(15,1)], [Q(8,2) Q(14,2); Q(9,2) Q(15,2)], mapx2(linspace(0,1,m2+1))', linspace(0,1,n2+1)', split);
[TT{8},  x{8},  y{8}] = triangulatepatch([Q(9,1) Q(15,1); Q(10,1) Q(16,1)], [Q(9,2) Q(15,2); Q(10,2) Q(16,2)], mapx3(linspace(0,1,m3+1))', linspace(0,1,n2+1)', split);
[TT{9},  x{9},  y{9}] = triangulatepatch([Q(10,1) Q(16,1); Q(11,1) Q(17,1)], [Q(10,2) Q(16,2); Q(11,2) Q(17,2)], mapx4(linspace(0,1,m4+1))', linspace(0,1,n2+1)', split);
[TT{10}, x{10}, y{10}] = triangulatepatch([Q(11,1) Q(17,1); Q(12,1) Q(18,1)], [Q(11,2) Q(17,2); Q(12,2) Q(18,2)], mapx5(linspace(0,1,m5+1))', linspace(0,1,n2+1)', split);

[TT{11}, x{11}, y{11}] = triangulatepatch([Q(16,1) Q(19,1); Q(17,1) Q(20,1)], [Q(16,2) Q(19,2); Q(17,2) Q(20,2)], mapx4(linspace(0,1,m4+1))', mapy(linspace(0,1,n3+1))', split);
[TT{12}, x{12}, y{12}] = triangulatepatch([Q(17,1) Q(20,1); Q(18,1) Q(21,1)], [Q(17,2) Q(20,2); Q(18,2) Q(21,2)], mapx5(linspace(0,1,m5+1))', mapy(linspace(0,1,n3+1))', split);

%trimesh(TT{12}, x{12}, y{12}, 'Color', 'k', 'LineWidth', 1);
%% combine mesh

% mesh coordinates
ndigits = 8;
ncomp = length(x);
count = 0;
P = [];
T = [];
for k=1:ncomp
    P = [P; round(x{k},ndigits) round(y{k},ndigits) zeros(length(y{k}), 1)];
    T = [T; TT{k}+count];
    count = count + length(y{k});
end

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

mesh.physobj(4).name = 'reservoir-left';
mesh.physobj(4).dim = 1;
mesh.physobj(4).tag = 4;

mesh.physobj(5).name = 'reservoir-right';
mesh.physobj(5).dim = 1;
mesh.physobj(5).tag = 5;

mesh.physobj(6).name = 'channel-top-right';
mesh.physobj(6).dim = 1;
mesh.physobj(6).tag = 6;

mesh.physobj(7).name = 'channel-top-center';
mesh.physobj(7).dim = 1;
mesh.physobj(7).tag = 7;

mesh.physobj(8).name = 'channel-top-left';
mesh.physobj(8).dim = 1;
mesh.physobj(8).tag = 8;

mesh.physobj(9).name = 'channel-inflow';
mesh.physobj(9).dim = 1;
mesh.physobj(9).tag = 9;

% add mesh entities
% mesh points
% bq = [1,2,3,4,5,6,12,18,21,20,19,16,15,14,13];
bq = [1,6,21,19,16,15,14,13];
k = 0;
for i=bq
    k = k+1;
    mesh.entity(k).dim = 0;
    mesh.entity(k).data = [Q(i,1) Q(i,2) 0];
end

% mesh lines
tol = 1e-5;

mesh.entity(k+1).dim = 1;
mesh.entity(k+1).data = extractboundary(mesh, Q(1,:)', Q(6,:)', tol);
mesh.entity(k+1).physobj = 'symmetry';

mesh.entity(k+2).dim = 1;
mesh.entity(k+2).data = extractboundary(mesh, Q(6,:)', Q(21,:)', tol);
mesh.entity(k+2).physobj = 'reservoir-right';

mesh.entity(k+3).dim = 1;
mesh.entity(k+3).data = extractboundary(mesh, Q(21,:)', Q(19,:)', tol);
mesh.entity(k+3).physobj = 'reservoir-top';

mesh.entity(k+4).dim = 1;
mesh.entity(k+4).data = extractboundary(mesh, Q(19,:)', Q(16,:)', tol);
mesh.entity(k+4).physobj = 'reservoir-left';

mesh.entity(k+5).dim = 1;
mesh.entity(k+5).data = extractboundary(mesh, Q(16,:)', Q(15,:)', tol);
mesh.entity(k+5).physobj = 'channel-top-right';

mesh.entity(k+6).dim = 1;
mesh.entity(k+6).data = extractboundary(mesh, Q(15,:)', Q(14,:)', tol);
mesh.entity(k+6).physobj = 'channel-top-center';

mesh.entity(k+7).dim = 1;
mesh.entity(k+7).data = extractboundary(mesh, Q(14,:)', Q(13,:)', tol);
mesh.entity(k+7).physobj = 'channel-top-left';

mesh.entity(k+8).dim = 1;
mesh.entity(k+8).data = extractboundary(mesh, Q(13,:)', Q(1,:)', tol);
mesh.entity(k+8).physobj = 'channel-inflow';

% mesh surfaces
k = k+9;
mesh.entity(k).dim = 2;
mesh.entity(k).data = T;
mesh.entity(k).bctag = 1:8;
mesh.entity(k).physobj = 'volume';

write_msh(mesh, 'nozzleflow.msh');

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
z = 1;
for i=1:length(mesh.entity)
    if mesh.entity(i).dim == 1
        v = mesh.entity(i).data;
        plot(P(v(:), 1), P(v(:), 2), 'Color', colors(z,:), 'LineWidth',4);
        z = z + 1;
    end
end
