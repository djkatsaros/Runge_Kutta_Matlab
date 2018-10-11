% Dean Katsaros
% Dec. 6, 2016
% Implementing 4th order runge kutta method for the ODE y' = ?5ty^2 + 5/t ?
% 1/t^2 , y(1) = 1 , 1 \leq t \leq T = 20

clc; close all; clear all; 
format long
  
allh = [ 0.2,0.1,0.05,0.02,0.01,0.005,0.002];     % determines number of divisions
errors_RK4 = [];                                % initialize blank vectors for error, orders, etc.
errors_BDF2 = [];                               
p_the_orders_BDF2 = [];                         
p_the_orders_RK4 = [];
mesh_sizes = 19./allh;                      % vector of all mesh_sizes for each h
F_ty = @(t,r) -5*t*r^2 + 5*t^(-1) - t^(-2);   % ODE, exact solution is y(t) = 1/t 
FDeriv = @(t,r) -10*t*r;                    % Derivative of ODE f(t,y)

for i  = 1:length(allh);
    
    h = allh(i);                            % particular h, mesh_size, N/2, etc. 
    N = mesh_sizes(i) ;
    No2 = floor(N/2);
    x = 1:h:20;                              % Mesh of t_i's from 1 to 20            
    y = zeros(1,length(x)); 
    y(1) = 1;                 % initialize, i.e., this is y_1         
    
    for j = 1:(length(x)-1)   
        
        k_1 = F_ty(x(j),y(j));                         % Calculate each K_i
        k_2 = F_ty(x(j)+0.5*h,y(j)+0.5*h*k_1);         % y(j) = y_j, x(j) = t_j
        k_3 = F_ty((x(j)+0.5*h),(y(j)+0.5*h*k_2));
        k_4 = F_ty((x(j)+h),(y(j)+h*k_3));

        y(j+1) = y(j) + (h/6)*(k_1 + 2*k_2 + 2*k_3 + k_4);  % RK4 equation
         
    end
    
    % calculate error via max_{n=N/2 ... N} |y_n - (y(t_n) = t_n^{-1})|
    pot_errors_RK4 = abs(y(No2:end) -x(No2:end).^(-1) );
    errors_RK4 = cat(2, errors_RK4, max(pot_errors_RK4));
    
    % Re - initialize for BDF 2, Backward euler or BDF 1 used to calculate
    % y_2 from y_1 and t_2, Newton's method used to find y_2
    y = zeros(1,length(x));   
    y(1) = 1; 
    r = zeros(1,100);
    r(1) = y(1);
    for l = 2:100;
        r(l) =  r(l-1) -(r(l-1) +  h*F_ty(x(2),r(l-1))-1)/(1+h*FDeriv(x(2),r(l-1)));   
        % Backward Euler formula, solved with newtons method
    end
    y(2) = r(100);
    
    % calculate $y_{k+1}$ using Newtons method, $x^(m+1) = x^m  - 
    % F(t_{k+1},x^m)/F'(t_{k+1},x^m)$, where F(t_k,x) = x
    % + (4/3)*y_{k-1} - (1/3)*y_{k-2} - (2*h/3)*f(t_{k+1},x)
    for k = 3:length(x);
        r = zeros(1,100);
        r(1) = y(k-1);   % r(l) = y_{k}^l for each k, i.e., the lth interation 
        %in newtons method for y_k. Implemented up to k = meshsize(i)
        
        for l = 2:100;
        r(l) = r(l-1)  + ((r(l-1) - (4/3)*y(k-1) + (1/3)*y(k-2) - (2*h/3)*F_ty(x(k),r(l-1)))/((2*h/3)*(-10*x(k-1)*r(l-1)) - 1 ));
        end
        
        % Solution for Newton's method iteration l
        y(k) = r(100);
        
    end
    
    % save all errors via formula as above
    pot_errors_BDF2 = abs(y(No2:end)- x(No2:end).^(-1));
    errors_BDF2 = cat(2, errors_BDF2, max(pot_errors_BDF2));
        
end

% log scale plot for the errors, both on the same plot
figure(1)
loglog( mesh_sizes, errors_RK4); hold on;
loglog(mesh_sizes, errors_BDF2)
title('log-log Scale Plots of Errors vs. the mesh size "N" ');
xlabel('Mesh Size');
ylabel('Error');

% Calculate errors 
order_BDF2 = log(errors_BDF2(1:end-1)./errors_BDF2(2:end))./log(mesh_sizes(2:end)./mesh_sizes(1:end-1));
order_RK4 = log(errors_RK4(1:end-1)./errors_RK4(2:end))./log(mesh_sizes(2:end)./mesh_sizes(1:end-1));

% plot the orders, take the last point of the lot as the order. 
figure(2)
plot( order_RK4 ); hold on;
plot( order_BDF2 );
title('Plot of the orders, calculated for each combination of n and m via p = ln(em/en)/ln(n/m)');
ylabel('order p ');
