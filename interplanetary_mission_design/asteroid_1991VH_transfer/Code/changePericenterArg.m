function [DeltaV, th_i, th_f] = changePericenterArg(a, e, om_i, om_f, mu)

% Change of Pericenter Argument maneuver
%
% [DeltaV,th_i,th_f] = changePericenterArg(a, e, om_i, om_f, mu)
%
%
% Input arguments:
% a       [1x1]    semi-major axis                [km]
% e       [1x1]    eccentricity                   [-]
% om_i    [1x1]    initial argument of periapsis  [rad]
% om_f    [1x1]    final argument of periapsis    [rad]
% mu      [1x1]    gravitational parameter        [km^3/s^2]
%
% Output arguments:
% DeltaV  [1x1]    maneuver impulse             [km/s]
% th_i    [2x1]    initial true anomalies       [rad]
% th_f    [2x1]    final true anomalies         [rad]



% This is a radial maneuver: its impulse does not depend on the choice
% between the 2 anomalies because the result is the same

% Variation of the argument of periapsis
deltaom = (om_f - om_i);

th_i = wrapTo2Pi([deltaom / 2, pi + deltaom / 2]); 
th_f = wrapTo2Pi([ - deltaom / 2, pi - deltaom / 2]);

% Semi-latus rectum
p = a * (1 - e ^ 2);

% The impulse depends on the variaton of the argument of periapsis
DeltaV = abs(2 * sqrt(mu / p) * e * sin(deltaom / 2));

end