function ii = gethalfedgesonaxis(he, P, v)
    epsilon = 1e-4;
    ii = find(abs(P(:,1)-v) < epsilon);
end