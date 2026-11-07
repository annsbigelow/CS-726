function [inform,x] = CG_PRplus(fun,x,PRparams)
    global numf numg
    
    [f,d,err,a0]=PRparams{:};
    
    success_flag = true;
    g=-d;
    k=0;
    while err > 10^-5*(1+abs(f))
        alpha_start=a0;
        g0 = g;
        [x,alpha,nf,ng,lsflag]=EBLS(fun,x,d,alpha_start);
        if lsflag==false
            success_flag=false;
            break
        end
    
        [f,g] = feval(fun,x,3); err=max(abs(g));
        beta_pr = (g'*(g-g0))/(g0'*g0);
        beta = max([beta_pr,0]);
        d=-g+beta*d;
        
        k=k+1;
        if k > 100000
            success_flag=false;
            break
        end
    end
    
    inform=[success_flag,k,feval(fun,x,1),err];
end