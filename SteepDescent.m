function [inform,x] = SteepDescent(fun,x,sdparams)
    [f,g]=feval(fun,x,3);
    err = max(abs(g));
    k=0;
    while err > 10^-5*(1+abs(f))
        err=max(abs(g));
        alpha_start=1;
        [x,alpha,nf,ng]=EBLS(fun,x,d,alpha_start);
        
        if numf > 100000
            fprintf('Function evaluations reached 100,000.\n');
            break
        end
        k=k+1;
       % fprintf('Steepest descent, current step: (%.5g,%.5g)\n',x(1),x(2));
    end
    fprintf('Steepest descent converged in %d steps.\n',k);
    
    fprintf('Total function evaluations: %d. Gradient evaluations: %d.\n',...
        numf,ng);

    inform=[, ,numf,numg,feval(fun,x,1),err,x];
end