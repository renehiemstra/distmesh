function V = extractboundary(mesh, a, b, tol)
    I = getnodesonaxis(mesh.vertices, a, b, tol);
    assert(size(I,1)>0, 'error: selection is empty.')
    A = mesh.incidence(I,:);
    v = full(sum(A, 1));
    J = find(v == 2);
    % node number of boundary component
    j = extractboundarycomponent(mesh.halfedges(mesh.edges(J,1),:));
    V = [j(1:end-1) j(2:end)];
end

function i = getnodesonaxis(P, a, b, tol)
    m = size(P,1);          % number of physical points
    X = P(:,1:2) - repmat(a',m,1); % relative coordinates of physical points
    u = (b-a)./norm(b-a);   % unit vector a -> b
    U = X * u;              % coordinates on line a -> b
    v = X - U * u';         % normal vector from points to line a -> b 
    V = sqrt(v(:,1).^2 + v(:,2).^2); % shortest distance from points to line a -> b
    U = U ./ norm(b-a);     % normalize coordinates in u-dir
    i = find(V < tol & U > -tol & U < 1.0+tol); % indices that satisfy path condition
end

function v = extractboundarycomponent(p)
    % size of numer of halfedges
    n = size(p, 1);
    assert(n>0, 'error: selection is empty.')
    % setup a linked list to extract the boundary curve
    list = java.util.LinkedList;
    % start in the middle and work forwards
    i = 1; j = 0;
    while true
        list.addLast(int32(p(i,1)));
        nextn = p(i,2);
        j = find(p(:,1) == nextn);
        if isscalar(j)==false
            break
        end
        i = j;
    end
    % we need to add one more point
    list.addLast(int32(p(i,2)));
    % start in the middle and work backwards
    nextn = list.getFirst();
    while true
        i = find(p(:,2) == nextn);
        if isscalar(i)==false
            break
        end
        list.addFirst(int32(p(i,1)));
        nextn = p(i,1);
    end
    % fill linked list into matlab array
    v = zeros(n,1);
    count = 1;
    while list.size() > 0
        v(count) = list.pop();
        count = count + 1;
    end
end


















function V_3D = extrude_boundary_edges(V_2D, N_pts)
    % V_2D: [m x 2] array of 2D edge node indices (e.g., from extractboundary)
    % N_pts: Total number of vertices in the 2D plane (offset for the 5-deg plane)
    % V_3D: [2*m x 3] array of 3D triangle node indices
    
    m = size(V_2D, 1);
    V_3D = zeros(m * 2, 3);
    
    % Base nodes
    n1 = V_2D(:, 1);
    n2 = V_2D(:, 2);
    
    % Top nodes (5-degree plane)
    n1_top = n1 + N_pts;
    n2_top = n2 + N_pts;
    
    % Split the extruded quadrilateral into two triangles
    % Triangle 1: [n1, n2, n1_top]
    V_3D(1:2:end, :) = [n1, n2, n1_top];
    
    % Triangle 2: [n2, n2_top, n1_top]
    V_3D(2:2:end, :) = [n2, n2_top, n1_top];
end