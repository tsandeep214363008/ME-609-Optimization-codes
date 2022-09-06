function [x_f,f_x_f,fun_eval] = Newtons_method(n,m,R,obj_function,constraint,x_0)

x=sym('x',[1,n]); % Defining the variables

k=0;% Iteration counter
% Defining the variables
syms alph g;
syms x_k x_new x_opt;
syms del_f_k;
syms H_i_k;
syms del_f_k_new;

feval=0;
x_k=x_0;

while k<=50 % Termination condition 1
    k=k+1; % Iteration counter
    
    % Penalty function
    f=obj_function;
    for r=1:m
        % Bracket operator
        if subs(constraint(r),x,x_k)<0 
            f=f+R*(constraint(r))^2;
        end
    end
    feval=feval+m; % function/bracket evaluations 

    %derivative of the function
    del_f=sym('del_f',[1,n]);
    for i=1:n
        del_f(i)=diff(f,x(i));
    end
    
    %Hessian inverse of the function
    H_i=inv(hessian(f,x));
    
    %substitution of values in the derivative
    del_f_k=vpa(subs(del_f,x,x_k)); feval=feval+2*n;
    % Termination condition 2
    if norm(del_f_k)<=0.001
        x_opt=x_k;
        break
    end    

    %substitution of values in the hessian
     H_i_k=vpa(subs(H_i,x,x_k)); feval=feval+3*n+4*((n*n-n)/2);
%     t=eig(H_i_k);
%     for i=1:n
%         if sign(t(i))==-1||imag(t(i))~=0
%             fprintf("Hessian is not a positive semi-definite matrix in the %dth iteration\n",k);
%             fprintf("Try different values of Initial guess\n" );
%             return
%         end    
%     end  

    %performing unidirectional search
    x_new=x_k-(alph*H_i_k*del_f_k.').'/norm(H_i_k*del_f_k.');
    p=subs(f,x,x_new);
    g=matlabFunction(subs(f,x,x_new));
    
    % Finding the unidirectional search direction using Bounding phase and
    % Golden section search method
        % Bounding phase method
        x0=0.5;delta=0.1;%Initial values to be considered
        % condition for Bounding Phase Method
        if ((g(x0+delta) >= g(x0))&&(g(x0-delta) <= g(x0)))
            delta = -delta;
        elseif ((g(x0+delta) <= g(x0))&&(g(x0-delta) >= g(x0)))
            delta = delta;
        else 
        % Alternate initial values to be considered
            x0=0.2;delta=0.1;
            if ((g(x0+delta) >= g(x0))&&(g(x0-delta) <= g(x0)))
                delta = -delta;
            elseif ((g(x0+delta) <= g(x0))&&(g(x0-delta) >= g(x0)))
                delta = delta;
            else
                fprintf("Try different values of Initial guess and increment for bounding phase method\n"); 
                return 
            end
        end
        j = 0; % Iteration counter
        feval=feval+3;
        x1 = x0 + delta;
        f0 = g(x0);
        f1 = g(x1);
        % While loop for getting the optimum interval
        while f0>f1 && j<3
            j = j + 1;
            y = x0;
            x0 = x1;
            x1 = x1+(2^j)*delta;
            f0 = g(x0);
            f1 = g(x1);
            feval=feval+1;
        end
        % Assigning the Interval extreme values to a and b 
        % This assigning will help us in continuing the Golden Section search Method with (a,b) as initial interval 
        if delta>0
            a=y;
            b=x1;
            %assignin('base','a',y); 
            %assignin('base','b',x1);
        else
            a=x1;
            b=y;
            %assignin('base','a',x1);
            %assignin('base','b',y);
        end
        % Golden section search method
        z1=a;
        z2=b;
        % Normalizing the interval between 0 and 1
        aw = 0;
        bw = 1;
        Lw = bw - aw;
        w1 = aw + .618 * Lw;
        w2 = bw - .618 * Lw;
        y1 = g(w1*(b-a)+a);
        y2 = g(w2*(b-a)+a);
        feval=feval+2;
        while abs(Lw) > 0.001   % Termination condition E=0.001(assumed)
            % Golden section search method condition
            % Removing the intervals
            if y1<y2
                aw=w2;
            elseif y1>y2
                bw=w1;
            else
                aw=w2;
                bw=w1;
            end
            Lw = bw - aw;
            w1 = aw + .618 * Lw;
            w2 = bw - .618 * Lw;
            y1 = g(w1*(b-a)+a);
            y2 = g(w2*(b-a)+a);
            z1=aw*(b-a)+a;
            z2=bw*(b-a)+a;
            feval=feval+1;
        end
        alph_opt=(z1+z2)*.5;
    % Completion of Unidirectional search    
    % New point
    x_new=vpa(subs(x_new,alph,alph_opt));
    del_f_k_new=vpa(subs(del_f,x,x_new)); feval=feval+2*n;

    % Termination conidition 3 (Linear independency check) 
    if abs(dot(del_f_k_new,del_f_k))<=0.00001
        x_opt=x_new;
        break 
    end

    % Termination conidition 4
    if (norm(x_new-x_k)/norm(x_k))<=0.001
        x_opt=x_new;
        break
    end

    x_k=x_new;

end

% Optimum point and Optimum value from the Newton's method
if k==51
    x_f=x_k;
else    
    x_f=x_opt;
end
f_x_f=vpa(subs(f,x,x_f));
% Function evaluations
fun_eval=feval;

end