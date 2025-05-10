clear; close all; clc;
do_plot = true;
s = 'Mesh: structured mesh on Reden testcase';

xa = 10;
ya = 5;
o = 20;
m = 70;
n = 10;

% mapping channel - part 1
xx = linspace(0,1,13);
yy = [0.0 0.2 0.5 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0];
map1x = @(x) pchip(xx, yy, x);
map1y = @(y) pchip([0.0 0.5 1.0], [0.0 0.35 0.5], y);
[T1, x, y] = mapsquare(m, n);
x1 = map1x(x); y1 = map1y(y);

% mapping channel - part 2
map2x = @(x) pchip(xx, 20-flip(yy), x);
map2y = @(y) pchip([0.0 0.5 1.0], [0.0 0.35 0.5], y);
[T2, x, y] = mapsquare(m, n);
x2 = map2x(x); y2 = map2y(y);

% mapping channel inlet
map3x = @(x) pchip([0.0 1/4 1/2 3/4 1.0], [-xa -2 -1 -0.25 0.0], x);
map3y = @(y) pchip([0.0 0.5 1.0], [0.0 0.35 0.5], y);
[T3, x, y] = mapsquare(o, n);
x3 = map3x(x); y3 = map3y(y);

map5x = @(x) pchip([0.0 1/4 1/2 3/4 1.0], 20 -flip([-xa -2 -1 -0.25 0.0]), x);
map5y = @(y) pchip([0.0 0.5 1.0], [0.0 0.35 0.5], y);
[T5, x, y] = mapsquare(o, n);
x5 = map5x(x); y5 = map5y(y);

% mapping channel above inlet
map4x = map3x;
map4y = @(y) pchip([0.0 0.5 1.0], [0.5 0.65 ya], y);
[T4, x, y] = mapsquare(o, n);
x4 = map4x(x); y4 = map4y(y);

map6x = @(x) pchip([0.0 1/4 1/2 3/4 1.0], 20-flip([-xa -2 -1 -0.25 0.0]), x);
map6y = @(y) pchip([0.0 0.5 1.0], [0.5 0.65 ya], y);
[T6, x, y] = mapsquare(o, n);
x6 = map6x(x); y6 = map6y(y);

%% combine mesh

% mesh coordinates
P = [x1 y1;
     x2 y2;
     x3 y3;
     x4 y4
     x5 y5;
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
  % axis off;
end

%% connectivity arrays

% number of elements
npts  = size(P, 1);
nelms = size(T, 1);

% halfedges
he = zeros(nelms*3, 2);
e  = zeros(1.2*ceil(nelms*3 / 2), 2);

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
    fprintf('%d\n', this);
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
P(bpts,:)

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

mesh.physobj(2).name = 'inlet';
mesh.physobj(2).dim = 1;
mesh.physobj(2).tag = 2;

mesh.physobj(3).name = 'outlet';
mesh.physobj(3).dim = 1;
mesh.physobj(3).tag = 3;

mesh.physobj(4).name = 'symmetry';
mesh.physobj(4).dim = 1;
mesh.physobj(4).tag = 4;

mesh.physobj(5).name = 'chamber1';
mesh.physobj(5).dim = 1;
mesh.physobj(5).tag = 5;

mesh.physobj(6).name = 'chamber2';
mesh.physobj(6).dim = 1;
mesh.physobj(6).tag = 6;

mesh.physobj(7).name = 'chamber3';
mesh.physobj(7).dim = 1;
mesh.physobj(7).tag = 7;

mesh.physobj(8).name = 'chamber4';
mesh.physobj(8).dim = 1;
mesh.physobj(8).tag = 8;

mesh.physobj(9).name = 'channel';
mesh.physobj(9).dim = 1;
mesh.physobj(9).tag = 9;

% add mesh entities
% mesh points
mesh.entity(1).dim = 0;
mesh.entity(1).data = [-xa 0 0];
mesh.entity(2).dim = 0;
mesh.entity(2).data = [o+xa 0 0];
mesh.entity(3).dim = 0;
mesh.entity(3).data = [-xa ya 0];
mesh.entity(4).dim = 0;
mesh.entity(4).data = [o+xa ya 0];
mesh.entity(5).dim = 0;
mesh.entity(5).data = [0 ya 0];
mesh.entity(6).dim = 0;
mesh.entity(6).data = [o ya 0];
mesh.entity(7).dim = 0;
mesh.entity(7).data = [0 0.5 0];
mesh.entity(8).dim = 0;
mesh.entity(8).data = [o 0.5 0];
% mesh lines
tol = 1e-4;
mesh.entity(9).dim = 1;
mesh.entity(9).data = extractboundary(mesh, @(x) abs(x+xa)<tol, @(y) y<ya+tol & y>-tol);
mesh.entity(9).physobj = 'inlet';
mesh.entity(10).dim = 1;
mesh.entity(10).data = extractboundary(mesh, @(x) x<o+xa+tol & x>-xa-tol, @(y) abs(y)<tol);
mesh.entity(10).physobj = 'symmetry';
mesh.entity(11).dim = 1;
mesh.entity(11).data = extractboundary(mesh, @(x) abs(x-(o+xa))<tol, @(y) y<ya+tol & y>-0.01);
mesh.entity(11).physobj = 'outlet';
mesh.entity(12).dim = 1;
mesh.entity(12).data = extractboundary(mesh, @(x) x<tol & x>-xa-tol, @(y) abs(y-ya)<tol);
mesh.entity(12).physobj = 'chamber4';
mesh.entity(13).dim = 1;
mesh.entity(13).data = extractboundary(mesh, @(x) x<o+xa+tol & x>o-tol, @(y) abs(y-ya)<tol);
mesh.entity(13).physobj = 'chamber1';
mesh.entity(14).dim = 1;
mesh.entity(14).data = extractboundary(mesh, @(x) abs(x)<tol, @(y) y<ya+tol & y>0.5-tol);
mesh.entity(14).physobj = 'chamber3';
mesh.entity(15).dim = 1;
mesh.entity(15).data = extractboundary(mesh, @(x) abs(x-o)<tol, @(y) y<ya+tol & y>0.5-tol);
mesh.entity(15).physobj = 'chamber2';
mesh.entity(16).dim = 1;
mesh.entity(16).data = extractboundary(mesh, @(x) x<o+tol & x>0-tol, @(y) abs(y-0.5)<tol);
mesh.entity(16).physobj = 'channel';

% mesh surfaces
mesh.entity(17).dim = 2;
mesh.entity(17).data = T;
mesh.entity(17).bctag = 1:8;
mesh.entity(17).physobj = 'volume';





write_msh(mesh, 'reden.msh');


%% plot boundaries
colors = [0 0.4470 0.7410; 
          0.8500 0.3250 0.0980;
          0.9290 0.6940 0.1250; 
          0.4940 0.1840 0.5560;
          0.4660 0.6740 0.1880;
          0.3010 0.7450 0.9330;
          0.6350 0.0780 0.1840;
          1 0 0];
k = 1;
for i=1:length(mesh.entity)
    if mesh.entity(i).dim ==1
        v = mesh.entity(i).data;
        plot(P(v(:), 1), P(v(:), 2), 'Color', colors(k,:), 'LineWidth',2);
        k = k + 1;
    end
end
