function [a, e, i, OM, om, th] = car2par(rr, vv, mu)

% TRASFORMATION FROM CARTESIAN COORDINATES TO KEPLERIAN PARAMETERS
%
% [a, e, i, OM, om, th] = car2par(rr, vv, mu)
% 
% 
% 
% Input arguments:
% rr          [3x1] position vector            [km]
% vv          [3x1] velocity vector            [km]
% mu          [1x1] gravitational parameter    [km^3/s^2]
%
% Output arguments:
% a           [1x1] semi_major axis            [km]
% e           [1x1] eccentricity               [-]
% i           [1x1] inclination                [rad]
% OM          [1x1] RAAN                       [rad]
% om          [1x1] argument of periapsis      [rad]
% th          [1x1] true anomaly               [rad]



% Compute modulus of the position and velocity vectors
r = norm(rr);
v = norm(vv);

% Specific orbital energy
E = 1/2 * v^2 - mu/r;

% Semi-major axis
a = - mu/(2*E);

% Specific relative angular momentum (vector and modulus)
hh = cross(rr,vv);
h = norm(hh);

% Eccentricity (vector and modulus)
ee = cross(vv, hh)/mu - rr/r;
e = norm(ee);

% Inclination
i = acos(hh(3)/h);

% Line of nodes
kk = [0, 0, 1];
N = cross(kk, hh)/norm(cross(kk, hh));

% RAAN
if N(2) >= 0
    OM = acos(N(1));
else
    OM = 2*pi - acos(N(1));
end

% Argument of periapsis
if ee(3) >= 0
    om = acos(dot(N, ee)/e);
else
    om = 2*pi - acos(dot(N,ee)/e);
end


v_r = dot(vv,rr)/r; % Radial velocity modulus

% True anomaly
if v_r >= 0
    th = acos(dot(rr,ee)/(r*e));
else
    th = 2*pi - acos(dot(rr,ee)/(r*e));
end
end





















