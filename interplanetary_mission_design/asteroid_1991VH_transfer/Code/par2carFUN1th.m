function [rr, vv] = par2carFUN1th(a ,e ,i ,OM, om, mu)

% TRASFORMATION FROM KEPLERIAN PARAMETERS TO CARTESIAN COORDINATES (AS A
% FUNCTION OF TRUE ANOMALY)
%
% [rr, vv] = par2carFUN1th(a ,e ,i ,OM, om, mu)
%
% Input arguments:
% a         [1x1] semi_major axis            [km]
% e         [1x1] eccentricity               [-]
% i         [1x1] inclination                [rad]
% OM        [1x1] RAAN                       [rad]
% om        [1x1] argument of the periapsis  [rad]
% mu        [1x1] gravitational parameter    [km^3/s^2]
%
% Output arguments:
% rr        [3x1] position vector            [km]
% vv        [3x1] velocity vector            [km/s]



% Semi-latus rectum
p =  a * (1 - (e)^2);

% Position vector modulus (as a function the the true anomaly)
r = @(th)  p/(1 + (e)*cos(th));

% Position and velocity vectors in the perifocal coordinate system (as a 
% function the the true anomaly)
r_pf =  @(th) r(th) .* [cos(th); sin(th); 0];
v_pf =  @(th) sqrt(mu/p) * [-sin(th); e + cos(th); 0];


% Definition of rotation matrices
R_OM = [cos(OM), sin(OM), 0; % Rotation through OM about the inertial Z-axis
       -sin(OM), cos(OM), 0;
          0,        0,    1];

R_i = [1,       0,       0;  % Rotation through i about the inertial X'-axis
       0,    cos(i),  sin(i);
       0,   -sin(i),  cos(i)];

R_om =  [cos(om), sin(om), 0; % Rotation through om about the inertial Z''-axis
        -sin(om), cos(om), 0;
            0,       0,    1];

% Total rotation matrix
T_ECItoPF =  R_om * R_i * R_OM;
T_PFtoECI =  (T_ECItoPF)';

% Position vector (as a function the the true anomaly)
rr = @(th) T_PFtoECI * r_pf(th);
% Velocity vector (as a function the the true anomaly)
vv = @(th) T_PFtoECI * v_pf(th);

end




