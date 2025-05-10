clear; close all; clc;
do_plot = true;
s = 'Mesh: Uniform mesh on Reden testcase';

xa = -1;
ya = 1;
o = 20;
m = 10;
n = 10;

% mapping channel
map1x = @(x) pchip([0.0 0.5 1.0], [0.5 2.0 10.0], x);
map1y = @(y) pchip([0.0 0.5 1.0], [0.0 0.35 0.5], y);

% create triangulation
[T1, x, y] = mapsquare(m, n);
x1 = map1x(x); y1 = map1y(y);
% trimesh(T1, map1x(x1), map1y(y1), 'Color', 'k', 'LineWidth', 2);
% hold on;

% mapping channel inlet
map2x = @(x) pchip([0.0 1.0], [0.0 0.5], x);
map2y = @(y) pchip([0.0 0.5 1.0], [0.0 0.35 0.5], y);

% create triangulation
[T2, x, y] = mapsquare(o, n);
x2 = map2x(x); y2 = map2y(y);
% trimesh(T2, map2x(x2), map2y(y2), 'Color', 'g', 'LineWidth', 2);

% mapping channel inlet
map3x = @(x) pchip([0.0 1.0], [xa 0.0], x);
map3y = @(y) pchip([0.0 0.5 1.0], [0.0 0.35 0.5], y);

% create triangulation
[T3, x, y] = mapsquare(o, n);
x3 = map3x(x); y3 = map3y(y);
% trimesh(T3, map3x(x3), map3y(y3), 'Color', 'r', 'LineWidth', 2);

% mapping channel above inlet
map4x = @(x) pchip([0.0 1.0], [xa 0.0], x);
map4y = @(y) pchip([0.0 0.5 1.0], [0.5 0.65 ya], y);

% create triangulation
[T4, x, y] = mapsquare(o, n);
x4 = map4x(x); y4 = map4y(y);
% trimesh(T4, map4x(x4), map4y(y4), 'Color', 'b', 'LineWidth', 2);

axis equal;

%% compute height function

% mesh density at cornernodes a-f
ha = 0.5;
hb = 1.0;
hc = 0.5;
hd = 0.2;
he = 0.1;
hf = 0.2;

% density functions at the boundary components
g{1} = @(x,y) pchip([0.0 0.5], [ha hb], y);
g{2} = @(x,y) pchip([0.5 1.0], [hb ha], y);
g{3} = @(x,y) pchip([xa 0], [ha ha], x);
g{4} = @(x,y) pchip([1 5], [hc hd], y);
g{5} = @(x,y) pchip([-5 0], [he hd], x);
g{6} = @(x,y) pchip([0 5], [hf he], y);
g{7} = @(x,y) pchip([-5 xa], [hf ha], x);

m = 10;
% boundary coordinates
p = [repmat(xa,1,m)'      linspace(0.0,0.5,m)'; 
     repmat(xa,1,m)'      linspace(0.5,ya,m)';
     linspace(xa,0,2*m)'    repmat(ya,1,2*m)';
     repmat(0,1,4*m)'       linspace(ya,5.0,4*m)';
     linspace(-5,0,5*m)'    repmat(5.0,1,5*m)';
     repmat(-5,1,5*m)'      linspace(0.0,5.0,5*m)';
     linspace(-5,-1,4*m)'   repmat(0.0,1,4*m)'];

% x- and y-coordinates
x = p(:,1); y = p(:,2); v = zeros(length(x), 1);

% boundary 1
for k=1:7
    ii = boundarycomponent(x, y, xa, ya, k);
    v(ii) = g{k}(x(ii), y(ii));
end

[X, Y] = meshgrid(linspace(-5, 0, 10), linspace(0, 5, 10));
V = griddata(x, y, v, X, Y);
mesh(X,Y,V); hold on;
plot3(x, y, v, 'o');


%% Unstructured part 5
do_plot = true;
s = 'Mesh: Uniform mesh on Reden testcase';

% fixed points
P =[-5                      0;
    -5                      5;
    0                       5;
    repmat(xa,1,n+1)'         map3y(linspace(0,1,n+1))';
    repmat(xa,1,n+1)'         map4y(linspace(0,1,n+1))';
    map4x(linspace(0,1,o+1))' repmat(ya,1,o+1)'];

pv = [-5.0  0.0;
       xa   0.0;
       xa   ya;
       0.0  ya;
       0.0  5.0;
      -5.0  5.0;
      -5.0  0.0];

disp(s)
fd = { 'l_dpolygon', [], pv };
fh = @(p) 1 ./ griddata(x, y, v, p(:,1), p(:,2));
[p,tri] = distmesh( fd, fh, 0.03, [-5,0; 0,5], P, [], 2000);

if( do_plot )
  clf
  c = [.9, .9, .9];
  patch( 'vertices', p, 'faces', tri, 'facecolor', c )
  hold on;
  patch( 'vertices', [x1,y1], 'faces', T1, 'facecolor', c )
  patch( 'vertices', [x2 y2], 'faces', T2, 'facecolor', c )
  patch( 'vertices', [x3 y3], 'faces', T3, 'facecolor', c )
  patch( 'vertices', [x4 y4], 'faces', T4, 'facecolor', c )
  title(s)
  axis tight
  axis equal
end
hold on;
plot(P(:,1),P(:,2),'o')
size(tri, 1) + size(T1, 1) + size(T2, 1) + size(T3, 1) + size(T4, 1)


