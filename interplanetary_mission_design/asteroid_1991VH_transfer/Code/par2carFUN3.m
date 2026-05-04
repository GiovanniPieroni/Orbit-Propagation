function [rr, vv] = par2carFUN3(a ,e ,i ,OM, TH, mu)

% TRASFORMATION FROM KEPLERIAN PARAMETERS TO CARTESIAN COORDINATES AS A 
% FUNCTION OF INITIAL TRUE ANOMALY, FINAL TRUE ANOMALY AND THE 
% ARGUMENT OF PERIAPSIS OF THE TRANSFER ORBIT DURING A DIRECT TRANSFER
% BETWEEN TWO ORBITS THAT DO NOT INTERSECT
%
% [rr, vv] = par2carFUN3(a ,e ,i ,OM, TH, mu)
%
% Input arguments:
% a           [1x1] semi_major axis            [km]
% e           [1x1] eccentricity               [-]
% i           [1x1] inclination                [rad]
% OM          [1x1] RAAN                       [rad]
% TH          [1x1] true anomaly               [rad]
% mu          [1x1] gravitational parameter    [km^3/s^2]
% Output arguments:
% rr          [3x1] position vector            [km]
% vv          [3x1] velocity vector            [km/s]



% semi.latus rectum (as a function of initial true anomaly, final true
% anomaly and the argument of periapsis of the transfer orbit)
p =  @(th1,th2, om) a(th1,th2,om) * (1 - (e(th1,th2,om))^2);

% Position vector modulus (as a function of initial true anomaly, final 
% true anomaly and the argument of periapsis of the transfer orbit)
r =  @(th1,th2,om) p(th1,th2,om)/(1 + (e(th1,th2,om))*cos(TH(th1,th2,om)));

% Position and velocity vectors in the perifocal coordinate system (as a 
% function of initial true anomaly, final true anomaly and the argument of 
% periapsis of the transfer orbit)
r_pf = @(th1,th2,om)  r(th1,th2,om) .* [cos(TH(th1,th2,om)); sin(TH(th1,th2,om)); 0];
v_pf = @(th1,th2,om) sqrt(mu/p(th1,th2,om)) * [-sin(TH(th1,th2,om)); e(th1,th2,om) + cos(TH(th1,th2,om)); 0];

% Definition of rotation matrices
R_OM = @(th1,th2) [cos(OM(th1,th2)), sin(OM(th1,th2)), 0;  % Rotation through 
                  -sin(OM(th1,th2)), cos(OM(th1,th2)), 0;  % OM about the 
                         0,                 0,         1]; % inertial Z-axis
                                                           % (as a function 
                                                           % of initial true 
                                                           % anomaly and final 
                                                           % true anomaly) 
                                                          
R_i = @(th1,th2) [1,           0,                0;        % Rotation through 
                  0,    cos(i(th1,th2)), sin(i(th1,th2));  % i about the 
                  0,   -sin(i(th1,th2)), cos(i(th1,th2))]; % inertial X'-axis 
                                                           % (as a function 
                                                           % of initial true 
                                                           % anomaly and final 
                                                           % true anomaly)

R_om = @(om) [cos(om), sin(om), 0;  % Rotation through om about the inertial 
             -sin(om), cos(om), 0;  % Z''-axis (as a function of the argument 
                 0,       0,    1]; % of periapsis of the transfer orbit)

% Total rotation matrix (as a function of initial true anomaly, final true
% anomaly and the argument of periapsis of the transfer orbit)
T_ECItoPF = @(th1,th2,om) R_om(om) * R_i(th1,th2) * R_OM(th1,th2);
T_PFtoECI = @(th1,th2,om) (T_ECItoPF(th1,th2,om))';

% Position vector (as a function of initial true anomaly, final true anomaly 
% and the argument of periapsis of the transfer orbit)
rr = @(th1,th2,om) T_PFtoECI(th1,th2,om) * r_pf(th1,th2,om);
% Velocity vector (as a function of initial true anomaly, final true anomaly 
% and the argument of periapsis of the transfer orbit)
vv = @(th1,th2,om) T_PFtoECI(th1,th2,om) * v_pf(th1,th2,om);

end
