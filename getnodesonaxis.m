function i = getnodesonaxis(P, fx, fy)
    i = find(fx(P(:,1)) & fy(P(:,2)));
end