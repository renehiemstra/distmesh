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