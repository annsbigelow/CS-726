clear
c1=.001; c2=0.5;
fun = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
grad = @(x) [-400*x(1)*(x(2)-x(1)^2)-2*(1-x(1)), 200*(x(2)-x(1)^2)];
tol = .0001;

for rho=[.1 .3 .5 .7 .9]
    linesearch = hw1_LinesearchClass(c1,c2,rho,0,0);
    x = [-.5 .5];
    strat=1;
    g=grad(x); linesearch.geval=1;
    f=fun(x); linesearch.eval=1;
    
    err = norm(g);
    k=0;
    while err > tol
        alpha_start=1;
        [x,f,g,alpha]=linesearch.backtrack(fun,grad,x,f,-g,alpha_start);
        err=norm(g);
        if linesearch.eval > 100000
            fprintf('Function evaluations reached 100,000.\n');
            break
        end
        k=k+1;
    end
    
    fprintf('beta=%g, steps=%d, total function evaluations = %d, gradient evaluations = %d.\n',...
        rho,k,linesearch.eval,linesearch.geval);
end