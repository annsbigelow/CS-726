classdef hw1_LinesearchClass < handle
    % This is built-in "handle" class so that eval can be accumulated
    properties
        c1
        c2
        rho
        eval % function evaluations
        geval % gradient evaluations
    end
    methods
        % Class constructor
        function obj = hw1_LinesearchClass(c1_,c2_,rho_,eval_,geval_)
            obj.c1 = c1_;
            obj.c2 = c2_;
            obj.rho = rho_;
            obj.eval=eval_;
            obj.geval=geval_;
        end

        function [x,f,g,alpha]=backtrack(obj,fun,grad,x0,f0,d,alpha_start)
            if alpha_start <= 0 
                error('Error: initial alpha is not positive.');
            end
            alpha=alpha_start;
            % f0 passed in to avoid evaluations at the same point through GD
            %   iterations
            x = x0 + alpha*d;
            f = feval(fun,x); obj.eval=obj.eval+1;

            g=-d;
            
            %k=0;
            while f > f0 + (obj.c1*alpha*dot(g,d))
                alpha = obj.rho*alpha;
                x = x0 + alpha*d;
                f = feval(fun,x); obj.eval=obj.eval+1;
                %fprintf('Backtracking trial step: (%g,%g,%g)\n',x(1),x(2),feval(fun,x));
                %k=k+1;
            end
            g = feval(grad,x); obj.geval=obj.geval+1;
            %fprintf('Final number of backtracking iterations: %d\n',k);
        end

        function z = phi(obj,a,x0,d,fun)
            z = feval(fun,x0+a*d);
            obj.eval=obj.eval+1;
        end
        function [z,g]=dphi(obj,a,x0,d,grad)
            g=feval(grad,x0+a*d);
            z=dot(d,g);
            obj.geval=obj.geval+1;
        end

        function [x,f,g,alpha]=EBLS(obj,fun,grad,x0,f0,d,alpha_start)
            l=0; u=1000000; % Idk what u should be
            alpha=alpha_start;
            g=-d;
            
            k=0;
            while true
                % Armijo check
                p=phi(obj,alpha,x0,d,fun);
                if p > f0 + obj.c1*alpha*dot(g,d)
                    u=alpha; alpha=(l+u)/2;
                    k=k+1;
                else
                    % Weak Wolfe check
                    [dp,gt]=dphi(obj,alpha,x0,d,grad);
                    if dp < obj.c2*dot(g,d)
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
                if k>25 % If # adjustments to alpha is greater than 25, break.
                    break;
                end
            end
            x=x0+alpha*d;
            f=p; g=gt;
        end %function
    end %methods
end %classdef