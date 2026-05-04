function [rr, vv] = par2car(a ,e ,i ,OM, om, th, mu)

% TRASFORMATION FROM KEPLERIAN PARAMETERS TO CARTESIAN COORDINATES
%
% [rr, vv] = par2car(a ,e ,i ,OM, om, th, mu)
%
% Input arguments:
% a       [1x1] semi_major axis            [km]
% e       [1x1] eccentricity               [-]
% i       [1x1] inclination                [rad]
% OM      [1x1] RAAN                       [rad]
% om      [1x1] argument of periapsis      [rad]
% th      [1x1] true anomaly               [rad]
% mu      [1x1] gravitational parameter    [km^3/s^2]
%
% Output arguments:
% rr      [3x1] position vector             [km]
% vv      [3x1] velocity vector             [km]



% Semi-latus rectum
p = a*(1-e^2);

% Position vector modulus
r = p/(1 + e*cos(th));

% Position and velocity vectors in the perifocal coordinate system
r_pf = r .* [cos(th); sin(th); 0];
v_pf = sqrt(mu/p) * [-sin(th); e + cos(th); 0];

% Define rotation matrices
R_OM = [cos(OM), sin(OM), 0; % Rotation through OM about the inertial Z-axis
       -sin(OM), cos(OM), 0;
           0,       0,    1];

R_i = [1,    0,      0;  % Rotation through i about the inertial X'-axis
       0,  cos(i), sin(i);
       0, -sin(i), cos(i)];

R_om = [cos(om), sin(om), 0; % Rotation through om about the inertial Z''-axis
       -sin(om), cos(om), 0;
           0,       0,    1];

% Total rotation matrix
T_ECItoPF = R_om * R_i * R_OM;
T_PFtoECI = T_ECItoPF';

% Position vector
rr = T_PFtoECI * r_pf;
% Velocity vector
vv = T_PFtoECI * v_pf;

end












