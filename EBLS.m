function [x,alpha,nf,ng]=EBLS(fun,x0,d,alpha_start)
    c1=.01; c2=0.3;
    l=0; u=1000000; % Idk what u should be
    alpha=alpha_start;
    
    fs=numf; gs=numg;

    f0=feval(fun,x0,1);
    g0=feval(fun,x0,2);
    k=0;
    while true
        % Armijo check
        if feval(fun,x0+alpha*d,1) > f0 + c1*alpha*dot(g0,d)
            u=alpha; alpha=(l+u)/2;
            k=k+1;
        else
            % Weak Wolfe check
            g=feval(fun,x0+alpha*d,2);
            if dot(d,g) < c2*dot(g0,d)
                l=alpha;
                if u > 10^5
                    alpha=2*l;
                else
                    alpha=(l+u)/2;
                end
                k=k+1;
            else 
                % Otherwise, accept alpha.
                break;
            end
        end
        if k>100 % If # adjustments to alpha is greater than 100, break.
            break;
        end
    end
    x=x0+alpha*d;
    nf=numf-fs; ng=numg-gs;
end