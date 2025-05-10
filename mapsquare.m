function [T, x, y] = mapsquare(m, n)
    % mesh nodes
    [Y, X] = meshgrid(linspace(0, 1, n+1), linspace(0, 1, m+1));
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
            % assign nodenumbers to triangle 1
            T(count, 1) = a; T(count, 2) = b; T(count, 3) = c;
            count = count + 1;
            % assign nodenumbers to triangle 2
            T(count, 1) = a; T(count, 2) = c; T(count, 3) = d;
            count = count + 1;
        end
    end
end