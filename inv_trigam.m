%Numerical approximation of the inverse trigamma function
function y = inv_trigam(x)
    lb = .001;
    ub = 100;
    y = mean([lb,ub]);
    z = psi(1,y);
    eps = .001;
    
    uz = psi(1,ub);
    while uz > x
        ub = ub*2;
        uz = psi(1,ub);
    end
    
    while abs((x-z)/x) > eps
        y = mean([lb,ub]);
        z = psi(1,y);
        if z < x
            ub = y;
        else
            lb = y;
        end
    end
    
end