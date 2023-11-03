% Cole Stoumbaugh
clear;clc
%% Givens
a = [10000 12000 40000];  % Semi Major axis
e = [.5 .1 .8]; % eccentricity

mu = 3.986005*10^5;

%% Calculations

% ODE45 Method
Q = numel(a);
q = 1;

while q < (Q+1)

    P = 2*pi*sqrt(a(1,q)^3/mu);  % period
    h= sqrt(mu*a(1,q)*(1 - e(1,q)^2));   % Specific angular momentum                                         
    R_p = a(1,q)*(1 - e(1,q));    % radius of perigee
    V_p = (h/R_p);  % perigee velocity
    dt = 0.01; % step size

    X(1) = R_p; % initial x-position
    X_dot(1) = 0;   % initial velocity in the x-direction 
    Y(1) = 0;   % initial y-position 
    Y_dot(1) = V_p; % initial velocity in the y-direction 

    Z_o = [X(1) Y(1) X_dot(1) Y_dot(1)];

    n = P/dt;
% Eulers    

    for t = 1:n
        R = sqrt(X(t)^2 + Y(t)^2);
        x_ddot = (-mu/R^3)*X(t);
        y_ddot = (-mu/R^3)*Y(t);
        X_dot(t+1) = X_dot(t)+x_ddot*dt;
        Y_dot(t+1) = Y_dot(t)+y_ddot*dt;
        X(t+1) = X(t) + X_dot(t)*dt;
        Y(t+1) = Y(t) + Y_dot(t)*dt;
    end

%% Plots
    t = [1, P];
    [~,Z] = ode45(@f2, t, Z_o);
    figure(q); 
    plot(Z(:,1),Z(:,2)); axis equal; grid on;
    title('ODE45 Method - Problem ' + sprintf("%d", q));
    xlabel('x-axis [km]');
    ylabel('y-axis [km]');

    % Euler's Method
    figure(q+3);
    plot(X,Y); axis equal;
    title('Eulers Method - Problem ' + sprintf("%d", q));
    xlabel('x-axis [km]');
    ylabel('y-axis [km]');

    q = q+1;
end

%% ODE Function
function Zdot = f2(~,Z)
    mu = 398600.5;
    R=sqrt(Z(1)^2 + Z(2)^2);     
    Zdot=[Z(3); Z(4); -mu*Z(1)/R^3; -mu*Z(2)/R^3;];
end
