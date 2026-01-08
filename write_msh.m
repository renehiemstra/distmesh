function write_msh(mesh,filename)
    fid = fopen(filename,'w');
    % start writing
    write_mesh_format(fid);
    write_physical_names(fid, mesh.physobj);
    mesh.entity = write_entities(fid, mesh.vertices, mesh.entity, mesh.physobj);
    write_vertices(fid, mesh.vertices, mesh.entity, mesh.boundary);
    write_elements(fid, mesh.entity);
    % stop writing
    fclose(fid);
end

function write_physical_names(fid, objects)
    m = length(objects);
    fprintf(fid,'%s\n', '$PhysicalNames');
    fprintf(fid,'%d\n', m);
    for i=1:m
        fprintf(fid,'%d %d \"%s\"\n', objects(i).dim, objects(i).tag, objects(i).name);
    end
    fprintf(fid,'%s\n', '$EndPhysicalNames');
end

function entity = write_entities(fid, vertices, entity, physobj)
    m = length(entity);
    fprintf(fid,'%s\n', '$Entities');
    d = getentitystats(entity);
    fprintf(fid,'%d %d %d %d\n', d(1), d(2), d(3), d(4));
    pointtag = 1;
    for i=1:m
        if entity(i).dim == 0
            write_point_entity(fid, vertices, entity(i).data, pointtag);
            entity(i).tag = pointtag;
            pointtag = pointtag + 1;
        end
    end
    linetag = 1;
    for i=1:m
        if entity(i).dim == 1
            physobjtag = getphysicalobjectid(physobj, entity(i).physobj);
            tag_a = getpointid(entity, vertices(entity(i).data(1), :));
            tag_b = getpointid(entity, vertices(entity(i).data(end), :));
            write_line_entity(fid, vertices, entity(i).data, linetag, physobjtag, tag_a, tag_b);
            entity(i).tag = linetag;
            linetag = linetag + 1;
        end
    end
    surfacetag = 1;
    for i=1:m
        if entity(i).dim == 2
            physobjtag = getphysicalobjectid(physobj, entity(i).physobj);
            write_surface_entity(fid, vertices, entity(i).data, entity(i).bctag, surfacetag, physobjtag);
            entity(i).tag = surfacetag;
            surfacetag = surfacetag + 1;
        end
    end
    fprintf(fid,'%s\n', '$EndEntities');
end

function id = getphysicalobjectid(physobj, name)
    id = 0;
    for i=1:length(physobj)
        if strcmp(physobj(i).name, name)
            id = physobj(i).tag;
        end
    end
    assert(id~=0, 'not a valid physical object name.')
end

function id = getpointid(entity, p)
    id = 0;
    for i=1:length(entity)
        if entity(i).dim == 0 && norm(entity(i).data - p) < 1e-5
            id = entity(i).tag;
            break
        end
    end
    assert(id~=0, sprintf('not a valid point entity with these coordinates: (%0.3f, %0.3f)', p(1), p(2)));
end

function write_point_entity(fid, vertices, data, tag)
    assert(size(data,1)==1 && size(data,2)==3, 'error: size of entity data should be (1,3).');
    getnodeid(vertices, data, 1e-5); % sanity check
    s = sprintf('%d', tag);
    p1 = sprintf('%0.5f', data(1));
    p2 = sprintf('%0.5f', data(2));
    p3 = sprintf('%0.5f', data(3));
    fprintf(fid,'%s%*s %s%*s %s%*s %s%*s %d\n',...
        s, 8-length(s), '', ...
        p1, 12-length(p1), '', ...
        p2, 12-length(p2), '', ...
        p3, 12-length(p3), '', ...
        0);
end

function write_line_entity(fid, vertices, data, tag, physobjtag, patag, pbtag)
    a = min(vertices(data, :));
    b = max(vertices(data, :));
    s = sprintf('%d', tag);
    a1 = sprintf('%0.5f', a(1));
    a2 = sprintf('%0.5f', a(2));
    a3 = sprintf('%0.5f', a(3));
    b1 = sprintf('%0.5f', b(1));
    b2 = sprintf('%0.5f', b(2));
    b3 = sprintf('%0.5f', b(3));
    fprintf(fid,'%s%*s %s%*s %s%*s %s%*s %s%*s %s%*s %s%*s %d %d %d %d %d\n',...
        s, 8-length(s), '', ...
        a1, 12-length(a1), '', ...
        a2, 12-length(a2), '', ...
        a3, 12-length(a3), '', ...
        b1, 12-length(b1), '', ...
        b2, 12-length(b2), '', ...
        b3, 12-length(b3), '', ...
        1, ...
        physobjtag, ...
        2,...
        -patag, ...
        pbtag);
end

function write_surface_entity(fid, vertices, data, bctag, tag, physobjtag)
    a = min(vertices(data(:), :));
    b = max(vertices(data(:), :));
    a1 = sprintf('%0.5f', a(1));
    a2 = sprintf('%0.5f', a(2));
    a3 = sprintf('%0.5f', a(3));
    b1 = sprintf('%0.5f', b(1));
    b2 = sprintf('%0.5f', b(2));
    b3 = sprintf('%0.5f', b(3));
    s = sprintf('%d', tag);
    fprintf(fid,'%s%*s %s%*s %s%*s %s%*s %s%*s %s%*s %s%*s %d %d %d %s\n',...
        s, 8-length(s), '', ...
        a1, 12-length(a1), '', ...
        a2, 12-length(a2), '', ...
        a3, 12-length(a3), '', ...
        b1, 12-length(b1), '', ...
        b2, 12-length(b2), '', ...
        b3, 12-length(b3), '', ...
        1, ...
        physobjtag, ...
        length(bctag),...
        num2str(bctag));
end

function b = comp(P, p, tol)
    b = abs(P-p)<tol;
end

function write_mesh_format(fid)
    fprintf(fid,'%s\n', '$MeshFormat');
    fprintf(fid,'%s\n', '4.1 0 8');
    fprintf(fid,'%s\n', '$EndMeshFormat');
end

function write_vertices(fid, vertices, entity, boundary)
    fprintf(fid,'%s\n', '$Nodes');
    n = size(vertices,1);
    fprintf(fid,'%d %d %d %d\n', length(entity), n, 1, n);
    for i=1:length(entity)
        if entity(i).dim == 0
            data = entity(i).data;
            fprintf(fid,'%d %d %d %d\n', entity(i).dim, entity(i).tag, 0, size(data,1));
            fprintf(fid, '%d\n', getnodeid(vertices, data, 1e-5));
            fprintf(fid,'%5.5f %5.5f %5.5f\n', data(1), data(2), data(3));
        elseif entity(i).dim == 1
            data = entity(i).data;
            fprintf(fid,'%d %d %d %d\n', entity(i).dim, entity(i).tag, 0, size(data,1)-1);
            for j=2:size(data,1)
                fprintf(fid,'%d\n', data(j,1));
            end
            for j=2:size(data,1)
                k = data(j,1);
                fprintf(fid,'%5.5f %5.5f %5.5f\n', vertices(k,1), vertices(k,2), vertices(k,3));
            end
        elseif entity(i).dim == 2
            fprintf(fid,'%d %d %d %d\n', entity(i).dim, entity(i).tag, 0, n - length(boundary));
            % node tags
            for j=1:size(vertices,1)
                if ismember(j, boundary)==false
                    fprintf(fid,'%d\n', j);
                end
            end
            % % node coordinates
            for j=1:size(vertices,1)
                if ismember(j, boundary)==false
                    fprintf(fid,'%5.5f %5.5f %5.5f\n', vertices(j,1), vertices(j,2), vertices(j,3));
                end
            end
        end
    end
    fprintf(fid,'%s\n', '$EndNodes');
end

function nodeid = getnodeid(vertices, data, tol)
    assert(size(data,1)==1 && size(data,2)==3);
    nodeid = find(comp(vertices(:,1), data(1), tol) & comp(vertices(:,2), data(2), tol) & comp(vertices(:,3), data(3), tol));
    assert(isscalar(nodeid) && int32(nodeid) == nodeid, sprintf('error: entity (%0.3f, %0.3f, %0.3f) does not correspond to a valid mesh vertex.', data(1), data(2), data(3)));
end

function write_elements(fid, entity)
    fprintf(fid,'%s\n', '$Elements');
    [nobjs, nelms] = getobjectstats(entity);
    fprintf(fid,'%d %d %d %d\n', nobjs, nelms, 1, nelms);
    etag = 1;
    for i=1:length(entity)
        if entity(i).dim > 0
            etag = write_element_group(fid, entity(i), etag);
        end
    end
    fprintf(fid,'%s\n', '$EndElements');
end

function [nentities, nelements] = getobjectstats(entity)
    nentities = 0;
    nelements = 0;
    for i=1:length(entity)
        if entity(i).dim>0
            nentities = nentities + 1;
            nelements = nelements + size(entity(i).data, 1);
        end
    end
end

function [dimobj] = getentitystats(entity)
    dimobj = zeros(4,1);
    for i=1:length(entity)
        d = entity(i).dim;
        dimobj(d+1) = dimobj(d+1) + 1; 
    end
end

function tag = write_element_group(fid, entity, tag)
    [m, n] = size(entity.data);
    fprintf(fid,'%d %d %d %d\n', entity.dim, entity.tag, n-1, m);
    data = entity.data;
    switch entity.dim
        case 1
            for i=1:m
                fprintf(fid,'%d %d %d\n', tag, data(i,1), data(i,2));
                tag = tag + 1;
            end
        case 2
            for i=1:m
                fprintf(fid,'%d %d %d %d\n', tag, data(i,1), data(i,2), data(i,3));
                tag = tag + 1;
            end
        otherwise
            error('not implemented.');
    end
end