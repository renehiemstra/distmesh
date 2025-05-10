function V = extractboundary(mesh, fx, fy)
    I = getnodesonaxis(mesh.vertices, fx, fy);
    assert(size(I,1)>0, 'error: selection is empty.')
    A = mesh.incidence(I,:);
    v = full(sum(A, 1));
    J = find(v == 2);
    % node number of boundary component
    j = extractboundarycomponent(mesh.halfedges(mesh.edges(J,1),:));
    V = [j(1:end-1) j(2:end)];
end