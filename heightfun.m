% mesh density function
function [h] = heightfun(p)
    h = zeros(size(p,1),1);
    for i=1:size(p,1)
        h(i) = evalheightfun(p(i,1),p(i,2));
    end
end


function [h] = evalheightfun(x,y)
    hmax = 2; hmin = 0.2; scal = 0.2;
    if (x > 0) && (x < 10)
        h = hmin;
    else
        h = min(hmin + scal * distfrominlet(x, y), hmax);
    end
end

% compute distance from line (10,0) - (10,0.5)
function d = distfrominlet(x, y)
    dx = x - 10;
    dy = min(y - 0.5, y - 0.0);
    d = sqrt(dx^2 + dy^2);
end
