clear 
global numf numg

n=4; m=n/4;
x0=zeros(n,1);
for i=0:m-1
    x0(4*i+1)=3; x0(4*i+2)=-2;
    x0(4*i+3)=0; x0(4*i+4)=1;
end

numf=0; numg=0;

fun = @xpowsing;
%% Steepest Descent
[f0,g0]=feval(fun,x0,3);
d0=-g0;
sdparams={f0,d0,max(abs(g0)),0.1};

[informsd,x] = SteepDescent(fun,x0,sdparams);

success_flag=informsd(1); k=informsd(2); 
final_f=informsd(3); final_err=informsd(4);

if success_flag==0
    fprintf(['Failure: steepest descent did not converge in less than ' ...
        '100,000 iterations (or the linesearch was terminated).\n']);
else
    fprintf(['SD converged in %d iterations: number of f evaluations = %.6g, ' ...
        'number of g evaluations = %.6g, f = %.6g, error = %.6g\n'], ...
        k, numf, numg, final_f, final_err);
    % fprintf('SD solution x = [');
    % fprintf('%.6g ',x);
    % fprintf(']\n');
end
%% Nonlinear CG: Fletcher-Reeves
numf=0; numg=0;

[f0,g0]=feval(fun,x0,3);
d0=-g0;
frparams={f0,d0,max(abs(g0)),0.1};

[informFR, x] = CG_FR(fun,x0,frparams);

success_flag=informFR(1); k=informFR(2); 
final_f=informFR(3); final_err=informFR(4);

if success_flag==0
    fprintf(['Failure: CG-FR did not converge in less than ' ...
        '100,000 iterations (or the linesearch was terminated).\n']);
else
    fprintf(['CG-FR converged in %d iterations: number of f evaluations = %.6g,' ...
        ' number of g evaluations = %.6g, f = %.6g, error = %.6g\n'], ...
        k, numf, numg, final_f, final_err);
    % fprintf('CG-FR solution x = [');
    % fprintf('%.6g ',x);
    % fprintf(']\n');
end
%% Nonlinear CG: Polak-Ribiere + 
numf=0; numg=0;

[f0,g0]=feval(fun,x0,3);
d0=-g0;
PRparams={f0,d0,max(abs(g0)),0.1};

[informPR, x] = CG_PRplus(fun,x0,PRparams);

success_flag=informPR(1); k=informPR(2); 
final_f=informPR(3); final_err=informPR(4);

if success_flag==0
    fprintf(['Failure: CG-PR+ did not converge in less than ' ...
        '100,000 iterations (or the linesearch was terminated).\n']);
else
    fprintf(['CG-PR+ converged in %d iterations: number of f evaluations = %.6g, ' ...
        'number of g evaluations = %.6g, f = %.6g, error = %.6g\n'], ...
        k, numf, numg, final_f, final_err);
    % fprintf('CG-PR+ solution x = [');
    % fprintf('%.6g ',x);
    % fprintf(']\n');
end