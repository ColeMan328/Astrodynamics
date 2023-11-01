% Cole Stoumbaugh
clear;clc


TLE = [1, 25544, 98067, 23107.54116911, .00020699, 00000+0, 37063-3, 0, 9998,...
        2, 25544, 51.6393, 269.0787, 0006070, 202.4487, 263.9445, 15.4991466039, 2381];



function startingVectors(UTC, TOF, a, E, FTA, i, RAAN, AOP, mu)
    GST = (UTC - (4*3600));
    R = (a*(1-E^2))/(1 + E*cosd(FTA));
    R_pqw = [R*cosd(FTA); R*sind(FTA); 0];
    V_pqw = sqrt(mu/(a*(1-E^2)))*[-sind(FTA); E+cosd(FTA); 0];
    %Rotation for Inclination
    ROT_i = [1 0 0;
              0 cosd(-i) -sind(-i);
              0 sind(-i) cosd(-i)];
    %Rotation for Raan
    ROT_RAAN = [cosd(-RAAN) -sind(-RAAN) 0;
                sind(-RAAN) cosd(-RAAN) 0;
                0 0 1];
    %Rotation for Argument of Perigee
    ROT_AOP = [cosd(-AOP) -sind(-AOP) 0;
                sind(-AOP) cosd(-AOP) 0;
                0 0 1];
    T = ROT_AOP*ROT_i*ROT_RAAN;
    %Updating R and V vectors
    R_Vector_2 = T\R_pqw;
    V_Vector_2 = T\V_pqw;
    fprintf("Final Range Tracking Site Location at: \n")
    fprintf("R = " + R_Vector_2(1) + " I + " + R_Vector_2(2) + " J + " + R_Vector_2(3) + " K \n");
    fprintf("V = " + V_Vector_2(1) + " I + " + V_Vector_2(2) + " J + " + V_Vector_2(3) + " K \n\n");  
    %finalTracking(L, longitude, H, GST, TOF, R_Vector_2)
    %findingSuitableTrackingStation(GST, TOF, R_Vector_2)
    %deltaV(sqrt(sum(V_Vector_2.^2)),sqrt(sum(R_Vector_2.^2)),20000)
end


function [Card_nummb1, Sat_numb1, Inter_desig, Epoch_time, Drag_rate, Drag_acceleration, Drag_term, Ephemeris_type, Element_numb,...
       Card_numb2, Sat_numb2, i, RAAN, e, w, mean_anomoly, mean_motion, revolution] = TLEread(TLE,TOF)
    
    Card_nummb1 = TLE(1);
    Sat_numb1 = TLE(2);
    Inter_desig = TLE(3);
    Epoch_time = TLE(4);
    Drag_rate = TLE(5);
    Drag_acceleration = TLE(6);
    Drag_term = TLE(7);
    Ephemeris_type = TLE(8);
    Element_numb = TLE(9);
    Card_numb2 = TLE(10);
    Sat_numb2 = TLE(11);
    i = TLE(12);
    RAAN = TLE(13);
    e = TLE(14);
    w = TLE(15);
    mean_anomoly = TLE(16);
    mean_motion = TLE(17);
    revolution = TLE(18);

    mu = 398600.5;
    a = (mu/((TLE(17)*((2*pi)/(3600*24)))^2))^(1/3);

    startingVectors(Epoch_time, TOF, a, e, mean_anomoly, i, RAAN, w, mu)
end
