% Cole Stoumbaugh
clear;clc

%% Given Values

R = [-6114 -1246 3373; 10000 0 0; -6998 4252 -3085; 15455 16005 0; 24912.16 0 0; 7199 9700 15940; 7199 9700 15940];
V = [-0.9811 -8.785 -3.633; 0 4.464 -4.464; 2.072 -10.319 -0.8566; -2 3 0; 0 4 0; 4.464 4.464 0; 2.464 2.464 1];

x = numel(R)/3;
y = 1;

mu = 398600.5;

k = [0 0 1]; % making k-hat matrix

%% Calculations

while y < (x+1)
    RV(y,1:3) = dot(R(y,1:3),V(y,1:3),2);
    H(y,1:3) = cross(R(y,1:3),V(y,1:3));
    N(y,1:3) = cross(k,H(y,1:3));

    r(y,1) = norm(R(y,1:3));
    vV(y,1) = norm(V(y,1:3));
    h(y,1) = norm(H(y,1:3));
    n(y,1) = norm(N(y,1:3));
    y = y+1;
end

ME = ((vV.^2)/2)-(mu./r);
E = ((cross(V,H))/mu)-(R./r);

hk = H(1:x,3);
ni = N(1:x,1);
nj = N(1:x,2);
ek = E(1:x,3);
rk = R(1:x,3);
ei = E(1:x,1);
ej = R(1:x,2);
ri = R(1:x,1);

y = 1;
while y < (x+1)
    e(y,1) = norm(E(y,1:3));
    y = y+1;
end

ER = dot(E,R,2);
NE = dot(N,E,2);
NR = dot(N,R,2);

a = abs((-mu./(2*ME))); % semi-major axis
i = acosd(hk./h); % inclination

y=1;
while y<=x %While loop to find remianing COE's, alternate COE's, and their respective half plane checks)

    raan(y,1) = acosd(ni(y,1)/n(y,1)); % raan
    if nj(y,1) < 0  % half-plane check
        raan(y,1) = 360 - raan(y,1);
    end

    rp(y,1) = abs(a(y,1)*(1-e(y,1)));

    w(y,1) = acosd(NE(y,1)/(n(y,1)*e(y,1)));  % argument of perigee
    if ek(y,1) < 0
        w(y,1) = 360 - w(y,1);
    end

    v(y,1) = acosd(ER(y,1)/(e(y,1)*r(y,1)));  % true anomaly
    if RV(y,1) < 0
        v(y,1) = 360 - v(y,1);
    end

    if i(y,1) == 0 || i(y,1) == 180 % existance check
        raan(y,1) = NaN;
        w(y,1) = NaN;
    end

    if -.001 < e(y,1) && e(y,1) < +.001
        v(y,1) = NaN;
        w(y,1) = NaN;
    end

    u(y,1) = NaN;   % Argument of Latitude

    Pi(y,1) = NaN;  % Longitude of Perigee
    
    l(y,1) = NaN;   % True longitude

    if i(y,1) < .001    % existance checks
        if -.001 < e(y,1) && e(y,1) < .001
            l(y,1) = acosd(ri(y,1)/r(y,1));
            if ej(y,1) < 0 
                l(y,1) = 360 - l(y,1);
            end
        else
            Pi(y,1) = acosd(ei(y,1)/e(y,1));
            if ej(y,1) < 0 
                Pi(y,1) = 360 - Pi(y,1);
            end
        end
    else
        if e(y,1) < .001
            u(y,1) = acosd(NR(y,1)/(n(y,1)*r(y,1)));
            if rk(y,1) < 0 
                u(y,1) = 360 - u(y,1);
            end
        end
    end
%% Orbit Types
    if e(y,1)<0.001 %Condition for circular orbit
        o = 'Orbit is circular';
        fprintf('Case # %d \nOrbit is circular\n',y)
    else
        if abs(1-e(y,1))<0.001 %Condition for parabolic orbit
            o = 'Orbit is parabolic';
            fprintf('Case # %d \nOrbit is parabolic\n',y)
        else
            if e(y,1)>1
                o = 'Orbit is hyperbola'; %Condition for hyperboolic orbit
                fprintf('Case # %d \nOrbit is hyperbolic\n',y)
            else
                o = 'Orbit is elliptical'; %Condition for elliptical orbit
                fprintf('Case # %d \nOrbit is ellipital\n',y)
            end
        end
    end

    %if e(y,1) < 0.001
        if i(y,1)==0 || i(y,1)==180 %Condition for equatorial orbit
            o = 'Orbit is equatorial';
            fprintf('and equatorial\n',y)
        end
    %end

    if rp(y,1) < 6628 % Condition for Ballistic
        o = 'Orbit is Ballistic';
        fprintf('and Ballistic\n',y)
    end
%% output
    
    Variable = categorical({'Case';'a';'e';'i';'RAAN';'Pi';'w';'u';'v';'l'});
    Value = [y;a(y,1);e(y,1);i(y,1);raan(y,1);Pi(y,1);w(y,1);u(y,1);v(y,1);l(y,1)];
    Tab = table(Variable, Value);
    writetable(Tab);
    disp(Tab);
    y = y+1;
end
