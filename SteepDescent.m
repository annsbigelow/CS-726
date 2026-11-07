function [inform,x] = SteepDescent(fun,x,sdparams)
    global numf numg

    [f,d,err,a0]=sdparams{:};
    
    success_flag = true;

    k=0;
    while err > 10^-5*(1+abs(f))
        alpha_start=a0;
        [x,alpha,nf,ng,lsflag]=EBLS(fun,x,d,alpha_start);
        if lsflag==false
            success_flag=false;
            break
        end
        
        [f,g] = feval(fun,x,3);
        d=-g;
        err=max(abs(g));

        k=k+1;
        if k > 100000
            success_flag=false;
            break
        end
    end

    inform=[success_flag,k,feval(fun,x,1),err];
end