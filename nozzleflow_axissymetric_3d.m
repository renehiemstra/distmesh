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

% number of mesh vertices
N_pts = size(mesh.vertices, 1);

% --- 1. DEFINE BOUNDARY TOPOLOGY ---
% Format: {Name, Start Point, End Point, is_symmetry_axis}
bnd_defs = {
    {'symmetry',           Q(1,:)',  Q(6,:)',  true},
    {'reservoir-right',    Q(6,:)',  Q(21,:)', false},
    {'reservoir-top',      Q(21,:)', Q(19,:)', false},
    {'reservoir-left',     Q(19,:)', Q(16,:)', false},
    {'channel-top-right',  Q(16,:)', Q(15,:)', false},
    {'channel-top-center', Q(15,:)', Q(14,:)', false},
    {'channel-top-left',   Q(14,:)', Q(13,:)', false},
    {'channel-inflow',     Q(13,:)', Q(1,:)',  false}
};

angle_deg = 5;
angle_rad = angle_deg * pi / 180;
N_pts = size(mesh.vertices, 1);
tol = 1e-5;

mesh.entity = [];
k = 0; % Entity counter

% --- 2. GENERATE 0D POINTS ---
bq = [1,6,21,19,16,15,14,13];
for i = bq
    % 0-degree point
    k = k+1; 
    mesh.entity(k).dim = 0; 
    mesh.entity(k).data = [Q(i,1), Q(i,2), 0]; 
    mesh.entity(k).physobj = sprintf('pt-%d-0deg', i);
    
    % 5-degree point
    k = k+1; 
    mesh.entity(k).dim = 0; 
    mesh.entity(k).data = [Q(i,1), Q(i,2)*cos(angle_rad), Q(i,2)*sin(angle_rad)]; 
    mesh.entity(k).physobj = sprintf('pt-%d-5deg', i);
end

% --- 3. GENERATE 1D LINES AND 2D SURFACES ---
for i = 1:length(bnd_defs)
    name = bnd_defs{i}{1}; pA = bnd_defs{i}{2}; pB = bnd_defs{i}{3}; is_sym = bnd_defs{i}{4};
    
    L_0deg = extractboundary(mesh, pA, pB, tol);
    L_5deg = L_0deg + N_pts;
    
    % 1D line on 0-deg plane
    k = k+1; mesh.entity(k).dim = 1; mesh.entity(k).data = L_0deg; mesh.entity(k).physobj = [name '-line-0deg'];
    
    if ~is_sym
        % 1D line on 5-deg plane
        k = k+1; mesh.entity(k).dim = 1; mesh.entity(k).data = L_5deg; mesh.entity(k).physobj = [name '-line-5deg'];
        % 2D Extruded Surface connecting them
        k = k+1; mesh.entity(k).dim = 2; mesh.entity(k).data = extrude_boundary_edges(L_0deg, N_pts); mesh.entity(k).physobj = name;
    end
end

% --- 4. 3D VERTICES & TETRAHEDRA (With orientation preservation) ---
P1 = mesh.vertices; 
P2 = zeros(size(P1)); P2(:,1) = P1(:,1); P2(:,2) = P1(:,2)*cos(angle_rad); P2(:,3) = P1(:,2)*sin(angle_rad);    
P_3D_raw = [P1; P2];

% Deduplicate vertices along the axis
ndigits = 8;
[mesh.vertices, ~, ic] = unique(round(P_3D_raw, ndigits), 'rows', 'stable');

nelms_2D = size(T, 1);
T_3D_raw = zeros(nelms_2D * 3, 4);
tet_count = 1;
for i = 1:nelms_2D
    base_nodes = T(i, :);
    [~, min_idx] = min(base_nodes);
    if min_idx == 1,     n1 = base_nodes(1); n2 = base_nodes(2); n3 = base_nodes(3);
    elseif min_idx == 2, n1 = base_nodes(2); n2 = base_nodes(3); n3 = base_nodes(1);
    else,                n1 = base_nodes(3); n2 = base_nodes(1); n3 = base_nodes(2); end
    
    n1_t = n1 + N_pts; n2_t = n2 + N_pts; n3_t = n3 + N_pts;
    if n2 < n3
        T_3D_raw(tet_count,:) = [n1, n2, n3, n1_t]; T_3D_raw(tet_count+1,:) = [n2, n3, n1_t, n2_t]; T_3D_raw(tet_count+2,:) = [n3, n1_t, n2_t, n3_t];
    else
        T_3D_raw(tet_count,:) = [n1, n2, n3, n1_t]; T_3D_raw(tet_count+1,:) = [n2, n3, n1_t, n3_t]; T_3D_raw(tet_count+2,:) = [n2, n1_t, n2_t, n3_t];
    end
    tet_count = tet_count + 3;
end

% --- 5. ADD SYMMETRY PLANES & VOLUME ---
k = k+1; mesh.entity(k).dim = 2; mesh.entity(k).data = [T(:,1), T(:,3), T(:,2)]; mesh.entity(k).physobj = 'sym-0deg';
k = k+1; mesh.entity(k).dim = 2; mesh.entity(k).data = T + N_pts; mesh.entity(k).physobj = 'sym-5deg';
k = k+1; mesh.entity(k).dim = 3; mesh.entity(k).data = T_3D_raw; mesh.entity(k).physobj = 'volume';

% --- 6. APPLY DEDUPLICATION MAPPING AND FILTER DEGENERATES ---
valid_entities = true(length(mesh.entity), 1);
for i = 1:length(mesh.entity)
    if mesh.entity(i).dim == 0
        % Snap to deduplicated coordinates
        mesh.entity(i).data = round(mesh.entity(i).data, ndigits);
        % If this is a 5-deg point and it lies on the axis (y=0), drop it to avoid duplicates
        if abs(mesh.entity(i).data(2)) < 1e-10 && contains(mesh.entity(i).physobj, '5deg')
            valid_entities(i) = false; 
        end
    else
        mapped_data = ic(mesh.entity(i).data);
        if mesh.entity(i).dim == 1
            valid = mapped_data(:,1) ~= mapped_data(:,2); % Drop lines that collapsed to points
            mesh.entity(i).data = mapped_data(valid, :);
        else
            valid = sum(diff(sort(mapped_data, 2), 1, 2) == 0, 2) == 0; % Drop surfaces/volumes that collapsed to lines/surfaces
            mesh.entity(i).data = mapped_data(valid, :);
        end
        if isempty(mesh.entity(i).data)
            valid_entities(i) = false;
        end
    end
    mesh.entity(i).bctag = []; 
end
mesh.entity = mesh.entity(valid_entities); % Remove empty entities

% --- 7. DYNAMICALLY BUILD PHYSICAL OBJECTS ---
% This cleanly prevents Gmsh element dimension mismatch errors forever!
mesh.physobj = [];
p_tag = 1;
for i = 1:length(mesh.entity)
    if isfield(mesh.entity(i), 'physobj') && ~isempty(mesh.entity(i).physobj)
        phys_name = mesh.entity(i).physobj;
        % Add it only if we haven't already defined this physical group
        if isempty(mesh.physobj) || ~any(strcmp({mesh.physobj.name}, phys_name))
            mesh.physobj(p_tag).name = phys_name;
            mesh.physobj(p_tag).dim = mesh.entity(i).dim;
            mesh.physobj(p_tag).tag = p_tag;
            p_tag = p_tag + 1;
        end
    end
end

% Write the Gmsh file
write_msh(mesh, 'nozzleflow_3D.msh');

%% 8. 3D VISUALIZATION IN MATLAB
figure('Name', '3D Full Topology Verification', 'Color', 'w');
hold on; grid on; view(3);

colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 
          0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880; 0.3010 0.7450 0.9330;
          0.6350 0.0780 0.1840; 0 0 1; 0 1 0; 1 0 0];
z_color = 1;

for i = 1:length(mesh.entity)
    % Plot 2D Surfaces
    if mesh.entity(i).dim == 2
        faces = mesh.entity(i).data;
        alpha_val = 0.9;
        if contains(mesh.entity(i).physobj, 'sym')
            alpha_val = 0.10; % Extremely transparent sym planes so we can see inside
        end
        patch('Vertices', mesh.vertices, 'Faces', faces, 'FaceColor', colors(z_color,:), ...
              'EdgeColor', 'none', 'FaceAlpha', alpha_val, 'DisplayName', mesh.entity(i).physobj);
        z_color = mod(z_color, size(colors,1)) + 1;
        
    % Plot 1D Lines (Boundaries)
    elseif mesh.entity(i).dim == 1
        edges = mesh.entity(i).data;
        X = [mesh.vertices(edges(:,1), 1), mesh.vertices(edges(:,2), 1)]';
        Y = [mesh.vertices(edges(:,1), 2), mesh.vertices(edges(:,2), 2)]';
        Z = [mesh.vertices(edges(:,1), 3), mesh.vertices(edges(:,2), 3)]';
        line(X, Y, Z, 'Color', 'k', 'LineWidth', 2, 'HandleVisibility', 'off');
        
    % Plot 0D Points (Corners)
    elseif mesh.entity(i).dim == 0
        pt = mesh.entity(i).data;
        plot3(pt(1), pt(2), pt(3), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'y', 'HandleVisibility', 'off');
    end
end

axis equal; xlabel('X'); ylabel('Y'); zlabel('Z');
title('Full 3D Extrusion (Surfaces, Lines, and Points)');
legend('show', 'Location', 'eastoutside');

function V_3D = extrude_boundary_edges(V_2D, N_pts)
    m = size(V_2D, 1);
    V_3D = zeros(m * 2, 3);
    n_a = V_2D(:, 1); n_b = V_2D(:, 2);
    n_a_top = n_a + N_pts; n_b_top = n_b + N_pts;
    for i = 1:m
        if n_a(i) < n_b(i)
            V_3D(2*i-1, :) = [n_a(i), n_b(i), n_a_top(i)];
            V_3D(2*i, :)   = [n_b(i), n_b_top(i), n_a_top(i)];
        else
            V_3D(2*i-1, :) = [n_a(i), n_b(i), n_b_top(i)];
            V_3D(2*i, :)   = [n_a(i), n_b_top(i), n_a_top(i)];
        end
    end
end