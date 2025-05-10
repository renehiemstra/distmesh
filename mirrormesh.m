function [P, T] = mirrormesh(P, T)
    P(:, 1) = -P(:,1);
end