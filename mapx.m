function y = mapx(x, xa, xb, xc, xd, ya, yb, yc, yd, t1, t2)
    n = length(x);
    y = zeros(n, 1);
    s1 = (yb-ya) / (xb-xa);
    s2 = (yc-yb) / (xc-xb);
    s3 = (yd-yc) / (xd-xc);
    for k=1:n
        y(k) = v(x(k), xa, xb, xc, xd, s1, s2, s3, ya, yb, yc, yd);
    end
end


function y = v(x, xa, xb, xc, xd, s1, s2, s3, ya, yb, yc, yd)
    if x >= xa && x < xb
        y = ya + s1 * (x-xa);
    elseif x >= xb && x < xc
        y = yb + s2 * (x-xb); 
    elseif x >= xc && x <= xd
        y = yc + s3 * (x-xc);
    else
        error("Invalid coordinate.")
    end
end

function y = f1(x, xa, xb, ya, yb, ta, tb)
    t = (x-xa) / (xb-xa);
    b = bernsteinMatrix(3, (x-xa) / (xb-xa));
    (yb - ya)
    [ya; ; ;yb]
end
