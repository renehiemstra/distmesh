function D = incidensemat(mesh)
    ne = size(mesh.edges, 1);
    nh = size(mesh.halfedges, 1);
    i = mesh.halfedges(mesh.edges(:,1),:)'; i = i(:);
    j = repmat(1:ne, 2, 1); j = j(:)';
    v = [repmat(1, 1, ne); repmat(1, 1, ne)]; v = v(:);
    D = sparse(i,j,v);
end