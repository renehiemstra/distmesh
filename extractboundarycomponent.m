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