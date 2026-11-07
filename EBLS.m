function [x,alpha,nf,ng,lsflag]=EBLS(fun,x0,d,alpha_start)
    lsflag=false;
    global numf numg
    c1=.5; c2=0.5;
    l=0; u=Inf;
    alpha=alpha_start;
    
    fs=numf; gs=numg;

    f0=feval(fun,x0,1);
    g0=feval(fun,x0,2);
    k=0;
    while true
        x=x0+alpha*d;
        % Armijo check
        if feval(fun,x,1) > f0 + c1*alpha*dot(g0,d)
            u=alpha; alpha=(l+u)/2;
            k=k+1;
        else
            % Weak Wolfe check
            g=feval(fun,x,2);
            if dot(d,g) < c2*dot(g0,d)
                l=alpha;
                if u == Inf
                    alpha=2*l;
                else
                    alpha=(l+u)/2;
                end
                k=k+1;
            else 
                % Otherwise, accept alpha.
                lsflag=true;
                break;
            end
        end
        if k>100 % If # adjustments to alpha is greater than 100, break.
            fprintf('Max alpha adjustments reached.\n');
            break;
        end
    end
    nf=numf-fs; ng=numg-gs;
end