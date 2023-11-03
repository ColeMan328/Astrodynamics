% Cole Stoumbaugh
clear;clc

%% Given Value
% Site Location
LST = [45; 262.5; 294.122]*(pi/180);
L = [30; 30; 38.83]*(pi/180);
He = [0; 2; 1.839];
% Sat Data
P = [638; 600; 430];
Az = [30; 270; 201.1]*(pi/180);
El = [90; 90; 87]*(pi/180);
dP = [0; 0; .5];
dAz = [0; 0; 2]*(pi/180);
dEl = [14.1; (0.015*(180/pi)); 1]*(pi/180);
% Constants
Re = 6378.137;
we = [0; 0; (15*(pi/(180*3600)))];
ee = 0.08182;

x = numel(P); 
y = 1;
z = 1;

%% Calculations
while y < (x+1)
    X(y,1) = (Re/sqrt(1-(ee^2)*sin(L(y,1))^2) + He(y,1))*cos(L(y,1));
    Z(y,1) = ((Re*(1-ee^2))/sqrt(1-(ee^2)*sin(L(y,1))^2) + He(y,1))*sin(L(y,1));

    Rs(y,1) = X(y,1)*cos(LST(y,1));
    Rs(y,2) = X(y,1)*sin(LST(y,1));
    Rs(y,3) = Z(y,1);

    ps(1,y) = -P(y,1)*cos(El(y,1))*cos(Az(y,1));
    ps(2,y) = P(y,1)*cos(El(y,1))*sin(Az(y,1));
    ps(3,y) = P(y,1)*sin(El(y,1));

    ROT(z:(z+2),1:3) = [sin(L(y,1))*cos(LST(y,1)) -sin(LST(y,1)) cos(L(y,1))*cos(LST(y,1));...
        sin(L(y,1))*sin(LST(y,1)) cos(LST(y,1)) cos(L(y,1))*sin(LST(y,1));...
        -cos(L(y,1)) 0 sin(L(y,1))];

    dps(1,y) = -dP(y,1)*cos(El(y,1))*cos(Az(y,1)) + P(y,1)*dEl(y,1)*sin(El(y,1))*cos(Az(y,1)) + P(y,1)*dAz(y,1)*cos(El(y,1))*sin(Az(y,1));
    dps(2,y) = -dP(y,1)*cos(El(y,1))*sin(Az(y,1)) - P(y,1)*dEl(y,1)*sin(El(y,1))*sin(Az(y,1)) + P(y,1)*dAz(y,1)*cos(El(y,1))*cos(Az(y,1));
    dps(3,y) = dP(y,1)*sin(El(y,1)) + P(y,1)*dEl(y,1)*cos(El(y,1));

    p(y,1:3) = ROT(z:(z+2),1:3)*ps(1:3,y);     
    R(y,1:3) = Rs(y,1:3) + p(y,1:3);

    dp(y,1:3) = ROT(z:(z+2),1:3)*dps(1:3,y);
    V(y,1:3) = dp(y,1:3) + cross(we,R(y,1:3));

    z = z + 3;
    y = y + 1;
end

x = numel(R)/3;
y = 1;

mu = 398600.5;

k = [0 0 1]; % making k-hat matrix

% Calculations

while y < (x+1)
    RV(y,1:3) = dot(R(y,1:3),V(y,1:3),2);
    H(y,1:3) = cross(R(y,1:3),V(y,1:3));
    N(y,1:3) = cross(k,H(y,1:3));

    r(y,1) = norm(R(y,1:3));
    v(y,1) = norm(V(y,1:3));
    h(y,1) = norm(H(y,1:3));
    n(y,1) = norm(N(y,1:3));
    ME(y,1) = ((v(y,1)^2)/2)-(mu/r(y,1));
    %E(y,1:3) = ((cross(V(y,1:3),H(y,1:3)))/mu)-(R(y,1:3)/r(y,1));
    E(y,1:3) = (1/mu)*((v(y,1)^2-(mu/r(y,1)))*R(y,1:3) - dot(R(y,1:3),V(y,1:3))*V(y,1:3));
    y = y+1;
end


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

a = -mu./(2*ME); % semi-major axis
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
    if RV < 0
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
% Orbit Types
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
% output
    fprintf("R = %.2f I + %.2f J + %.2f K \n",R(y,1),R(y,2),R(y,3))
    fprintf("V = %.2f I + %.2f J + %.2f K \n",V(y,1),V(y,2),V(y,3))
    Variable = categorical({'Case';'a';'e';'i';'RAAN';'Pi';'w';'u';'v';'l'});
    Value = [y;a(y,1);e(y,1);i(y,1);raan(y,1);Pi(y,1);w(y,1);u(y,1);v(y,1);l(y,1)];
    Tab = table(Variable, Value);
    writetable(Tab);
    disp(Tab);
    y = y+1;
end
