clear
c1=.001; c2=0.5; rho=.5;
linesearch = hw1_LinesearchClass(c1,c2,rho,0,0);

fun = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
grad = @(x) [-400*x(1)*(x(2)-x(1)^2)-2*(1-x(1)), 200*(x(2)-x(1)^2)];

tol = .0001;

%% Strategy 1 for choosing alpha_start:
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
   % fprintf('Steepest descent, current step: (%.5g,%.5g)\n',x(1),x(2));
end
fprintf('Steepest descent (backtrack) converged in %d steps.\n',k);

fprintf('Total function evaluations: %d. Gradient evaluations: %d.\n',...
    linesearch.eval,linesearch.geval);
fprintf('Strategy %d. beta=%g.\n',strat,rho);
%% Strategy 2 for choosing alpha_start:
x = [-.5 .5];
strat=2;
g=grad(x); linesearch.geval=1;
f=fun(x); linesearch.eval=1;

err = norm(g);
eta=2.0;
k=0;
while err > tol
    if k==0
        alpha_start=1;
    else
        alpha_start=eta*alpha;
    end
    [x,f,g,alpha]=linesearch.backtrack(fun,grad,x,f,-g,alpha_start);
    err=norm(g);
    if linesearch.eval > 100000
        fprintf('Function evaluations reached 100,000.\n');
        break
    end
    k=k+1;
   % fprintf('Steepest descent, current step: (%.5g,%.5g)\n',x(1),x(2));
end
fprintf('Steepest descent (backtrack) converged in %d steps.\n',k);

fprintf('Total function evaluations: %d. Gradient evaluations: %d.\n',...
    linesearch.eval,linesearch.geval);
fprintf('Strategy %d, beta=%g, eta=%g\n',strat,rho,eta);