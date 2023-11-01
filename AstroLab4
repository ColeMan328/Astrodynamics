% Cole Stoumbaugh
clear;clc

% Given Values
LST = [45; 262.5; 294.122]*(pi/180);
L = [30; 30; 38.83]*(pi/180);
He = [0; 2; 1.839];
Long = 0;

P = [638; 600; 430];
Az = [30; 270; 201.1]*(pi/180);
El = [90; 90; 87]*(pi/180);
dP = [0; 0; .5];
dAz = [0; 0; 2]*(pi/180);
dEl = [14.1; (0.015*(180/pi)); 1]*(pi/180);

Re = 6378.1;
we = [0; 0; (7.292082578*10^(-5))];

x = numel(P); 
y = 1;
z = 1;

% Calculations
while y < (x+1)
    Rs(y,1) = (Re+He(y,1))*cos(L(y,1))*cos(LST(y,1));
    Rs(y,2) = (Re+He(y,1))*cos(L(y,1))*sin(LST(y,1));
    Rs(y,3) = (Re+He(y,1))*sin(L(y,1));

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
zz = 1;

mu = 398600.5;

k = [0 0 1]; % making k-hat matrix
ee = 0.08182;

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
    E(y,1:3) = ((cross(V(y,1:3),H(y,1:3)))/mu)-(R(y,1:3)/r(y,1));
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
    
    Variable = categorical({'Case';'a';'e';'i';'RAAN';'Pi';'w';'u';'v';'l'});
    Value = [y;a(y,1);e(y,1);i(y,1);raan(y,1);Pi(y,1);w(y,1);u(y,1);v(y,1);l(y,1)];
    Tab = table(Variable, Value);
    writetable(Tab);
    disp(Tab);

    % Kepler type 2 problem: find v_f with known TOF and v_i 
    
    Ei(y,1) = acosd((e(y,1) + cosd(v(y,1)))/(1 + e(y,1)*cosd(v(y,1)))); % degrees
    % quad check
    if v(y,1) > 180
        Ei(y,1) = 360 - Ei(y,1);
    end
    
    Ei(y,1) = Ei(y,1)*(pi/180);    % radians
    Mi(y,1) = Ei(y,1) - e(y,1)*sin(Ei(y,1)); % radians 
    
    n(y,1) = sqrt(mu/a(y,1)^3); % mean motion
    Mf(y,1) = Mi(y,1) + n(y,1)*TOF(y,1);
    k(y,1) = floor(Mf(y,1)/(2*pi));
    Mf(y,1) = Mf(y,1) - k*2*pi;    % radians
    if Mf(y,1) > 2*pi
         Mf(y,1) = Mf - k*2*pi;
    end
    
    % Use Newton's iteration method
    M0(y,1) = Mf(y,1) - e(y,1)*sin(Mf(y,1));
    E0(y,1) = Mf(y,1) + (Mf(y,1) - M0(y,1))/(1 - e(y,1)*cos(Mf(y,1)));
    Mn(y,1) = E0(y,1) - e(y,1)*sin(E0(y,1));
    En(y,1) = E0(y,1) + (Mf(y,1) - Mn(y,1))/(1 - e(y,1)*cos(E0(y,1)));
    eps(y,1) = En(y,1)-E0(y,1);
    
    while abs(eps(y,1)) >= 0.00001
        Mn(y,1) = En(y,1) - e(y,1)*sin(En(y,1));
        EnPlus1(y,1) = En(y,1) + (Mf(y,1) - Mn(y,1))/(1 - e(y,1)*cos(En(y,1)));
        
        eps(y,1) = EnPlus1(y,1)-En(y,1);
        En(y,1) = EnPlus1(y,1);
    end
    
    Ef(y,1) = En(y,1);
    vf(y,1) = acosd((cos(Ef(y,1))-e(y,1))/(1-e(y,1)*cos(Ef(y,1)))); % degrees
    % quad check
    if Ef(y,1) > pi
        vf(y,1) = 360 - vf(y,1);
    end
    fprintf(fid,'v:\t\t Final True Anomaly\t = %0.3f [degrees]\n',vf);
    
    %% Display type of orbit
    if e(y,1) < 0.001
        fprintf(fid,'\n\nThe orbit is circular');
    else
        if 0 < e(y,1) && e(y,1) < 1 && abs(1-e(y,1)) > 0.001
            fprintf(fid,'\nThe orbit is elliptical');
        else
            if abs(1-e(y,1)) < 0.001
                fprintf(fid,'\nThe orbit is parabolic');
        else
            if e(y,1) > 1
                  fprintf(fid,'\nThe orbit is hyperbolic');
            end
            end
        end
    end
    
    %% Find position and velocity vectors
    
    r_new(y,1) = (a(y,1)*(1-e(y,1)^2))/(1+e(y,1)*cosd(vf(y,1)));
    R_pqw(1:3,y) = [r_new(y,1)*cosd(vf(y,1)) r_new(y,1)*sind(vf(y,1)) 0]';
    v_pqw(1:3,y) = sqrt(mu/(a(y,1)*(1-e(y,1)^2)))*[-sind(vf(y,1)) e+cosd(vf(y,1)) 0]';
    
    ROT1i(y,1) = [1 0 0;0 cosd(-i(y,1)) -sind(-i(y,1));0 sind(-i(y,1)) cosd(-i(y,1))];
    ROT3raan(y,1) = [cosd(-raan(y,1)) -sind(-raan(y,1)) 0;sind(-raan(y,1)) cosd(-raan(y,1)) 0;0 0 1];
    ROT3omega(y,1) = [cosd(-w(y,1)) -sind(-w(y,1)) 0;sind(-w(y,1)) cosd(-w(y,1)) 0;0 0 1];
    T(y,1) = ROT3omega*ROT1i*ROT3raan;
    
    Rn(y,1:3) = inv(T(y,1))*R_pqw(y,1);
    Vn(y,1:3) = inv(T(y,1))*v_pqw(y,1);
    
    fprintf('\n\n\n(b) Updated position and velocity vectors:\n\n');
    fprintf('R = %0.3f I + %0.3f J + %0.3f K [km]\n',Rn(y,1),Rn(y,2),Rn(y,3));
    fprintf('V = %0.3f I + %0.3f J + %0.3f K [km/s]\n',v_ijk(y,1),v_ijk(y,2),v_ijk(y,3));
    
    %% Find rho_ijk and new Az, El, and range
    
    Pn(y,1:3) = Rn(y,1) - Rs(y,1:3);
    rho_ijk(y,1) = norm(Pn(y,1:3));
    Psn(y,1:3) = inv(ROT(zz:zz+2,1:3))*Pn(y,1:3);
    
    Pnew(y,1) = norm(Psn(y,1:3));
    Elnew(y,1) = asind(Psn(y,3)/Pnew(y,1));
    Aznew(y,1) = acosd(-Psn(y,1)/(Pnew(y,1)*cosd(Elnew(y,1))));
    
    if Psn(y,2) < 0
        Aznew(y,1) = 360 - Aznew(y,1);
    end
    
    
    %% Final tracking station
    % Tracking from same radar site
    
    GST_f(y,1) = GST(y,1) + TOF_hr(y,1)*15;
    LST_f(y,1) = GST_f(y,1) + Long(y,1);
    
    x(y,1) = ((Re/sqrt(1-ee(y,1)^2*(sind(L(y,1)))^2))+He(y,1))*cosd(L(y,1));
    z(y,1) = (((Re*(1-ee(y,1)^2))/sqrt(1-ee(y,1)^2*(sind(L(y,1)))^2))+He(y,1))*sind(L(y,1));
    
    Rsnew(y,1:3) = [x*cosd(LST_f), x*sind(LST_f), z];
    
    rho_ijk_new(y,1:3) = Rn(y,1:3) - Rsnew(y,1:3);

    ROTn(zz,zz+2,1:3) = [sind(L(y,1))*cosd(LST_f(y,1)), -sind(LST_f(y,1)), cosd(L(y,1))*cosd(LST_f(y,1));...
        sind(L(y,1))*sind(LST_f(y,1)), cosd(LST_f(y,1)), cosd(L(y,1))*sind(LST_f(y,1));...
        -cosd(L(y,1)), 0, sind(L(y,1))];
    
    fprintf('\n\n(c) Satellites location at the final site.\n');
    
    rho_sez_bar_new(y,1:3) = inv(ROTn(zz,zz+2,1:3))*rho_ijk_new(y,1:3);
    rho_f(y,1) = norm(rho_sez_bar_new(y,1:3));
    El_f(y,1) = asind(rho_sez_bar_new(y,3)/rho_f(y,1));
    Az_f(y,1) = acosd(-rho_sez_bar_new(y,1)/(rho_f(y,1)*cosd(El_f(y,1))));

    if rho_sez_bar_new(y,2) < 0
        Az_f(y,1) = 360-Az_f(y,1);
    end

    fprintf('\nThe final range, elevation, and azimuth are:\n');
    fprintf('\n rho:\t\t range\t\t\t = %0.3f [km]', rho_f(y,1));
    fprintf('\n Az:\t\t azimuth\t\t\t = %0.3f [degrees]', Az_f(y,1));
    fprintf('\n El:\t\t elevation\t\t\t = %0.3f [degrees]', El_f(y,1));

    y = y + 1;
    zz = zz + 3;
end
