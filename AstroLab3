% Cole Stoumbaugh
clear;clc

%% Given Values
LST = [250.5; 262.5; 294.122]*(pi/180);
L = [30; 30; 38.83]*(pi/180);
He = [2; 2; 1.839];

P = [5000;600;430]; % 600; 430];
Az = [45; 270; 201.1]*(pi/180);
El = [90; 90; 87]*(pi/180);
dP = [0; 0; .5];
dAz = [0; 0; 2]*(pi/180);
dEl = [1; (0.015*(180/pi)); 1]*(pi/180);

Re = 6378.1;
we = [0; 0; (7.292082578*10^(-5))];

x = numel(P);
y = 1;
z = 1;

%% Calculations
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
