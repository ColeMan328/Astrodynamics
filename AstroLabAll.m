% Cole Stoumbaugh
clear;clc

%% Given Values
LST = [45; 262.5; 294.122]; % Local Sidereal Time
L = [30; 30; 38.83];    % Launch Site latitude
H = [0; 2; 1.839];

P = [683; 600; 430];    % Range (distance from station to satellite)
Az = [30; 270; 201.1];  % Azimuth (angle measured from north clockwise)
El = [90; 90; 87];  % Elevation (angle measured from locat horizion)
dP = [0; 0; .5];    % Range rate
dAz = [0; 0; 2];    % Rate of change of azimuth
dEl = [14.1; 0.015; 1]; % Rate of change of elevation

Re = 6378.1;    % Radius of earth
w = [0; 0; 15.041]; % angular velocity of earth

x = numel(P);
y = 1;
z = 1;

%% Calculations
while y < (x+1)
    Rs(y,1) = -Re*cosd(L(y,1))*cosd(LST(y,1));
    Rs(y,2) = Re*cosd(L(y,1))*sind(LST(y,1));
    Rs(y,3) = Re*sind(L(y,1));

    ps(1,y) = -P(y,1)*cosd(El(y,1))*cosd(Az(y,1));
    ps(2,y) = P(y,1)*cosd(El(y,1))*sind(Az(y,1));
    ps(3,y) = P(y,1)*sind(El(y,1));

    ROT(z:(z+2),1:3) = [sind(L(y,1))*cosd(LST(y,1)) -sind(LST(y,1)) cosd(L(y,1))*cosd(LST(y,1));...
        sind(L(y,1))*sind(LST(y,1)) cosd(LST(y,1)) cosd(L(y,1))*sind(LST(y,1));...
        -cosd(L(y,1)) 0 sind(LST(y,1))];

    dps(1,y) = -dP(y,1)*cosd(El(y,1))*cosd(Az(y,1)) + P(y,1)*dEl(y,1)*sind(El(y,1))*cosd(Az(y,1)) + P(y,1)*dAz(y,1)*cosd(El(y,1))*sind(Az(y,1));
    dps(2,y) = -dP(y,1)*cosd(El(y,1))*sind(Az(y,1)) - P(y,1)*dEl(y,1)*sind(El(y,1))*sind(Az(y,1)) + P(y,1)*dAz(y,1)*cosd(El(y,1))*cosd(Az(y,1));
    dps(3,y) = dP(y,1)*sind(El(y,1)) + P(y,1)*dEl(y,1)*cosd(El(y,1));

    p(y,1:3) = ROT(z:(z+2),1:3)*ps(1:3,y);
    R(1:3,y) = Rs(y,1:3) + p(y,1:3);

    dp(1:3,y) = ROT(z:(z+2),1:3)*dps(1:3,y);
    V(y,1:3) = dp(1:3,y) + cross(w,R(1:3,y));

    z = z + 3;
    y = y + 1;
end


x = numel(R)/3;
y = 1;

mu = 398600.5;

k = zeros(x,3); % making k-hat matrix
k(1:x,3) = 1;

%% Calculations

RV = dot(R,V,2);
H = cross(R,V);
N = cross(k,H);

while y < (x+1)
    r(y,1) = norm(R(y,1:3));
    v(y,1) = norm(V(y,1:3));
    h(y,1) = norm(H(y,1:3));
    n(y,1) = norm(N(y,1:3));
    y = y+1;
end

ME = ((v.^2)/2)-(mu./r);
E = ((cross(V,H))/mu)-(R./r);

hk = H(1:x,3);
ni = N(1:x,1);
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
    if ek(y,1) > 0  % half-plane check
        raan(y,1) = 360 - raan(y,1);
    end

    rp(y,1) = abs(a(y,1)*(1-e(y,1)));

    w(y,1) = acosd(NE(y,1)/(n(y,1)*e(y,1)));  % argument of perigee
    if ek > 0
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

    if e(y,1) < 0.001
        if i(y,1)==0 %Condition for equatorial orbit
            o = 'Orbit is equatorial';
            fprintf('and equatorial\n',y)
        end
    end

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
