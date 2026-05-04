function [rr, vv] = par2carFUN1(a ,e ,i ,OM, th, mu)

% TRASFORMATION FROM KEPLERIAN PARAMETERS TO CARTESIAN COORDINATES EXPRESSED
% AS A FUNCTION OF THE ARGUMENT OF PERIAPSIS OF THE TRANSFER ORBIT DURING A 
% DIRECT TRANSFER BETWEEN TWO ORBITS THAT DO NOT INTERSECT
%
% [rr, vv] = par2carFUN1(a ,e ,i ,OM, th, mu)
%
% Input arguments:
% a       [1x1] semi_major axis            [km]
% e       [1x1] eccentricity               [-]
% i       [1x1] inclination                [rad]
% OM      [1x1] RAAN                       [rad]
% th      [1x1] true anomaly               [rad]
% mu      [1x1] gravitational parameter    [km^3/s^2]
%
% Output arguments:
% rr      [3x1] position vector             [km]
% vv      [3x1] velocity vector             [km/s]



% Semi-latus rectum
p =  @(om) a(om) * (1 - (e(om))^2);

% Position vector modulus (as a function of the argument of periapsis of
% the transfer orbit)
r =  @(om) p(om)/(1 + (e(om))*cos(th(om)));

% Position and velocity vectors in the perifocal coordinate system (as 
% functions of the argument of periapsis of the transfer orbit)
r_pf = @(om)  r(om) .* [cos(th(om)); sin(th(om)); 0];
v_pf = @(om) sqrt(mu/p(om)) * [-sin(th(om)); e(om) + cos(th(om)); 0];

% Definition of rotation matrices
R_OM = [cos(OM), sin(OM), 0; % Rotation through OM about the inertial Z-axis
       -sin(OM), cos(OM), 0;
          0,        0,    1];

R_i = [1,      0,       0;  % Rotation through i about the inertial X'-axis
       0,    cos(i),  sin(i);
       0,   -sin(i),  cos(i)];

R_om = @(om) [cos(om), sin(om), 0;  % Rotation through om about the inertial 
             -sin(om), cos(om), 0;  % Z''-axis (as a function of the argument 
                 0,       0,    1]; % of periapsis of the transfer orbit)
             
                 

% Total rotation matrix (as a function of the argument of periapsis of the 
% transfer orbit)
T_ECItoPF = @(om) R_om(om) * R_i * R_OM;
T_PFtoECI = @(om) (T_ECItoPF(om))';

% Position vector (as a function of the argument of periapsis of
% the transfer orbit)
rr = @(om) T_PFtoECI(om) * r_pf(om);
% Velocity vector (as a function of the argument of periapsis of
% the transfer orbit)
vv = @(om) T_PFtoECI(om) * v_pf(om);

end




