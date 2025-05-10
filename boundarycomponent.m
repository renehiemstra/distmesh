function b = boundarycomponent(x, y, xa, ya, k)
    epsilon = 0.001;
    if k==1 
        b = abs(x-xa) < epsilon & y>=0-epsilon & y<=0.5;
    elseif k==2
        b = abs(x-xa) < epsilon & y>=0.5 & y<=ya+epsilon;
    elseif k==3
        b = x>=xa-epsilon & x<=0.0+epsilon & abs(y-ya) < epsilon;
    elseif k==4
        b = abs(x) < epsilon & y>=1.0 & y<=5.0;
    elseif k==5
        b = x>=-5.0 & x<=0.0 & abs(y-5.0) < epsilon;
    elseif k==6
        b = abs(x+5.0) < epsilon & y>=0 & y<=5.0;
    elseif k==7
        b = x>=-5.0 & x<=-1 & abs(y) < epsilon;
    end
end
