function [T, x, y] = triangulatepatch(Px, Py, u, v, split)

    % mesh nodes
    X = bilinearmap(Px, u, v);
    Y = bilinearmap(Py, u, v);
    
    if split
        [T,x,y] = face_vertex_map_with_split(X, Y);
    else
        [T,x,y] = face_vertex_map(X, Y);
    end
end

function [T,x,y] = face_vertex_map(X, Y)
% get sizes of mesh in each dimension
    M = size(X,1); % number of vertices in x-dir
    N = size(X,2); % number of vertices in y-dir
    m = M-1; % number of faces in x-dir
    n = N-1; % number of faces in y-dir

    x = X(:);
    y = Y(:);
    
    % assign node numbers to triangles
    T = zeros(m*n, 3);
    count = 1;
    for j = 1:n
        for i = 1:m
            
            % corners of square i*j
            a = i + (j-1)*(m+1);
            b = a + 1;
            d = a + m+1;
            c = d + 1;
            
            % assign nodenumbers to triangle T1
            T(count, 1) = a; T(count, 2) = b; T(count, 3) = c;
            count = count + 1;
            % assign nodenumbers to triangle T2
            T(count, 1) = a; T(count, 2) = c; T(count, 3) = d;
            count = count + 1;
            
        end
    end
end

function [T,x,y] = face_vertex_map_with_split(X, Y)

    % get sizes of mesh in each dimension
    M = size(X,1); % number of vertices in x-dir
    N = size(X,2); % number of vertices in y-dir
    m = M-1; % number of faces in x-dir
    n = N-1; % number of faces in y-dir

    x = [X(:); 0; 0];
    y = [Y(:); 0; 0];
    
    % assign node numbers to triangles
    T = zeros(m*n+4, 3);
    count = 1;
    for j = 1:n
        for i = 1:m
            % corners of square i*j
            a = i + (j-1)*(m+1);
            b = a + 1;
            d = a + m+1;
            c = d + 1;
            % special case where T1 needs to be broken up into
            % several triangles to make sure a triangle only touches
            % one boundary
            if (i==m) && (j==1)
                % compute center of T1 and add to x-y-vertices
                p = M*N+1;
                x(p) = (x(a) + x(b) + x(c)) / 3;
                y(p) = (y(a) + y(b) + y(c)) / 3;
                % add subtriangles of triangle T1
                T(count, 1) = a; T(count, 2) = b; T(count, 3) = p;
                count = count + 1;
                T(count, 1) = b; T(count, 2) = c; T(count, 3) = p;
                count = count + 1;
                T(count, 1) = c; T(count, 2) = a; T(count, 3) = p;
                count = count + 1;
                % assign nodenumbers to triangle T2
                T(count, 1) = a; T(count, 2) = c; T(count, 3) = d;
                count = count + 1;

            % special case where T2 needs to be broken up into
            % several triangles to make sure a triangle only touches
            % one boundary
            elseif (i==1) && (j==n)
                % compute center of T1 and add to x-y-vertices
                p = M*N+2;
                x(p) = (x(a) + x(c) + x(d)) / 3;
                y(p) = (y(a) + y(c) + y(d)) / 3;
                % assign nodenumbers to triangle T1
                T(count, 1) = a; T(count, 2) = b; T(count, 3) = c;
                count = count + 1;
                % add subtriangles of triangle T2
                T(count, 1) = a; T(count, 2) = c; T(count, 3) = p;
                count = count + 1;
                T(count, 1) = c; T(count, 2) = d; T(count, 3) = p;
                count = count + 1;
                T(count, 1) = d; T(count, 2) = a; T(count, 3) = p;
                count = count + 1;
            else
                % assign nodenumbers to triangle T1
                T(count, 1) = a; T(count, 2) = b; T(count, 3) = c;
                count = count + 1;
                % assign nodenumbers to triangle T2
                T(count, 1) = a; T(count, 2) = c; T(count, 3) = d;
                count = count + 1;
            end
        end
    end
end

function x = bilinearmap(X, u, v)
    x = [u 1-u] * X * [v 1-v]';
end