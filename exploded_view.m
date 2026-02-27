% ---------------------------------------------------------
% 14-Sector Velocity Space Decomposition (Exploded View)
% ---------------------------------------------------------

% 1. Define the 4 normal vectors of a regular spatial tetrahedron
% These represent the 4 zero-flux boundary planes: v * n_i = 0
N = [ 1,  1,  1;
      1, -1, -1;
     -1,  1, -1;
     -1, -1,  1]' / sqrt(3);

% 2. Generate a dense point cloud on the surface of a sphere (v_max = 1)
% We use 10,000 points to ensure the curved spherical caps look smooth.
[az, el] = meshgrid(linspace(0, 2*pi, 100), linspace(-pi/2, pi/2, 100));
pts = [cos(el(:)).*cos(az(:)), cos(el(:)).*sin(az(:)), sin(el(:))];

% 3. Calculate the "signature" of each point
% The signature is the sign of the dot product with each of the 4 normals.
dots = pts * N;
sig = sign(dots);
sig(sig == 0) = 1; % Safeguard for floating point exact zeros

% Convert the 4-part signature (-1 or 1) into a unique integer ID
% Because sum(N) = 0, this will mathematically yield exactly 14 unique IDs.
sig_id = (sig(:,1)>0)*8 + (sig(:,2)>0)*4 + (sig(:,3)>0)*2 + (sig(:,4)>0)*1;
unique_ids = unique(sig_id);

% 4. Prepare the figure
figure('Color', 'w', 'Name', '14-Sector Velocity Space', 'Position', [100, 100, 800, 800]);
hold on; view(3); grid on; axis equal; axis off;
cameratoolbar; % Allows you to easily click and drag to rotate the 3D view

% Aesthetic settings
explode_factor = 0.3; % Distance to pull the pieces apart radially
colors = lines(length(unique_ids)); % 14 distinct colors

% 5. Build and plot each of the 14 polyhedral cones
for i = 1:length(unique_ids)
    id = unique_ids(i);
    
    % Extract the surface points belonging strictly to this specific cone
    cone_surf_pts = pts(sig_id == id, :);
    
    if size(cone_surf_pts, 1) < 4
        continue; % Skip if mathematically degenerate (safeguard)
    end
    
    % Calculate the center of the spherical cap to dictate the explosion direction
    dir = mean(cone_surf_pts, 1);
    dir = dir / norm(dir);
    
    % Shift the surface points outward
    shifted_surf = cone_surf_pts + explode_factor * dir;
    
    % The apex of every pyramid is the origin [0,0,0]. Shift it outward too.
    shifted_origin = [0, 0, 0] + explode_factor * dir;
    
    % Combine the apex and the spherical cap points
    region_pts = [shifted_origin; shifted_surf];
    
    % Generate the solid 3D volume using a convex hull
    K = convhull(region_pts(:,1), region_pts(:,2), region_pts(:,3));
    
    % Plot the polyhedral cone
    trisurf(K, region_pts(:,1), region_pts(:,2), region_pts(:,3), ...
        'FaceColor', colors(i,:), 'FaceAlpha', 0.85, ...
        'EdgeColor', 'k', 'EdgeAlpha', 0.2, 'LineWidth', 0.5);
end

% 6. Lighting and annotations
camlight right;
camlight left;
lighting gouraud;
title('Exploded View of the 14 Velocity Sectors', 'FontSize', 16, 'FontWeight', 'bold');
subtitle('8 Triangular Pyramids (Face/Vertex)  |  6 Quadrangular Pyramids (Edge)', 'FontSize', 12);