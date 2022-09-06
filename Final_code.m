
% Initial Guess
input_0 = inputdlg({'Enter Initial guess as space-separated numbers ='}, 'Input values', [1 50]); 
x_0=str2num(input_0{1}); 
R=0.1; % Initial penalty parameter

%Objective function1
question=1;
n=2; % No. of variables
m=6; % Total No. of total constraints
s=4; % No. of variable bounds
x=sym('x',[1,n]); % Creating symbolic variables
constraint=sym('constraint',[1,m]);
% Objective function, constraints and variable bounds 
% Also normalization is performed alongside
obj_function=(x(1)-10)^3+(x(2)-20)^3; %  obj_function=obj_function/7973;
constraint(1)=(x(1)-5)^2+(x(2)-5)^2-100;    constraint(1)=constraint(1)/100;
constraint(2)=-1*((x(1)-6)^2+(x(2)-5)^2-82.81); constraint(2)=constraint(2)/82.81; 
constraint(3)=x(1)-13;  constraint(3)=constraint(3)/13;
constraint(4)=20-x(1);  constraint(4)=constraint(4)/20;
constraint(5)=x(2);
constraint(6)=4-x(2);   constraint(6)=constraint(6)/4;

%Objective function2
% question=2;
% n=2; % No. of variables 
% m=6; % Total No. of total constraints
% s=4; % No. of variable bounds
% x=sym('x',[1,n]); % Creating symbolic variables
% constraint=sym('constraint',[1,m]);
% % Objective function, constraints and variable bounds 
% % Also normalization is done alongside
% obj_function=-1*(sin(2*pi*x(2))*(sin(2*pi*x(1)))^3)/((x(1))^3*(x(1)+x(2))); 
% constraint(1)=-1*(x(1)^2-x(2)+1);
% constraint(2)=-1*(1-x(1)+(x(2)-4)^2);
% constraint(3)=x(1);
% constraint(4)=10-x(1);  constraint(4)=constraint(4)/10;
% constraint(5)=x(2);
% constraint(6)=10-x(2);  constraint(6)=constraint(6)/10;

% Objective function3
% question=3;
% n=8; % No. of variables 
% m=22; % Total No. of total constraints
% s=16; % No. of variable bounds
% x=sym('x',[1,n]); % Creating symbolic variables
% constraint=sym('constraint',[1,m]);
% % Objective function, constraints and variable bounds 
% % Also normalization is done alongside
% obj_function=x(1)+x(2)+x(3);    obj_function=obj_function/30000;
% constraint(1)=-1*(-1+0.0025*(x(4)+x(6)));
% constraint(2)=-1*(-1+0.0025*(-x(4)+x(5)+x(7)));
% constraint(3)=-1*(-1+0.01*(-x(6)+x(8)));
% constraint(4)=-1*(100*x(1)-x(1)*x(6)+833.33252*x(4)-83333.333);  constraint(4)=constraint(4)/83333.333;
% constraint(5)=-1*(x(2)*x(4)-x(2)*x(7)-1250*x(4)+1250*x(5));
% constraint(6)=-1*(x(3)*x(5)-x(3)*x(8)-2500*x(5)+1250000);  constraint(6)=constraint(6)/1250000;
% constraint(7)=x(1)-100; constraint(7)=constraint(7)/100;
% constraint(8)=10000-x(1);   constraint(8)=constraint(8)/10000;   
% constraint(9)=x(2)-1000; constraint(9)=constraint(9)/1000;
% constraint(10)=10000-x(2);   constraint(10)=constraint(10)/10000;  
% constraint(11)=x(3)-1000; constraint(11)=constraint(11)/1000;
% constraint(12)=10000-x(3);   constraint(12)=constraint(12)/10000; 
% k=13;
% for j=4:8
%     constraint(k)=x(j)-10; constraint(k)=constraint(k)/10;
%     constraint(k+1)=1000-x(j);   constraint(k+1)=constraint(k+1)/1000; 
%     k=k+2;
% end

% Creating empty matrices to store the values later
P=[];
constraint_violation=[];
f_opt=[];
feval_iter=[];
x_k=x_0;
E=1; % Considered to start the while loop  
t=0;
feval=0;

while abs(E)>0.1  % Termination condition 1
    
    t=t+1;% Iteration counter

    % Unconstrained optimization is performed using newton's method
    [x_f,f_x_f,fun_eval] = Newtons_method(n,m,R,obj_function,constraint,x_k);  
    
    % Optimum point and Optimum value from the Newton's method
    x_k=x_f;
    P(t)=f_x_f;

    % Function evaluations
    feval=feval+fun_eval;
    feval_iter(t)=fun_eval;

    % Constraint violation calculation for every constraint 
    % and for every iteration
    for i=1:m-s
        if subs(constraint(i),x,x_k)<0
            constraint_violation(i,t)=vpa(subs(constraint(i),x,x_k));
        else
            constraint_violation(i,t)=0;
        end
    end
    
    % Comparing Penalizing functions
    if t==1
        E=1;
    else
        E=P(t)-P(t-1);
    end
    
    % Upadating Penalty parameter
    R=10*R; % c=10 is considered
    
    % Denormalizing the objective functions
    if question==1
        f_opt(t)=subs(obj_function,x,x_k);
    end
    if question==2
        f_opt(t)=-subs(obj_function,x,x_k);
    end    
    if question==3
        f_opt(t)=subs(obj_function,x,x_k)*3000;
    end

    % Termination condition 2 
    if question==1
        if R==Inf||t==15
            break
        end 
    else
        if R==Inf
            break
        end 
    end
end
  
fprintf('The optimum point of the function: [');
fprintf('%2f', x_k);
fprintf(']\n');
fprintf('The optimum value of the function: %f\n',f_opt(t));
fprintf('No. of iterations: %d\n',t);
fprintf('Total function evaluations: %d\n',feval);

% Convergence plots
% figure(1)
% plot1=plot(1:length(f_opt),f_opt);
% hold on
% scatter(1:length(f_opt),f_opt)  
% xlabel('Iterations')
% ylabel('Function Value')
% ax = plot1.Parent;  
% set(ax, 'XTick', 0:1:length(f_opt))
% 
% figure(2)
% plot2=plot(1:length(feval_iter),feval_iter);
% hold on
% scatter(1:length(feval_iter),feval_iter)  
% xlabel('Iterations')
% ylabel('Functional evaluations')
% ax = plot2.Parent;  
% set(ax, 'XTick', 0:1:length(feval_iter))



