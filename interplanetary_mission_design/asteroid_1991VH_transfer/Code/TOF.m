function deltat = TOF(a, e, th1, th2, mu)

% Time of Flight
% 
% deltat = TOF(a, e, th1, th2)
% 
% 
% Input arguments:
% a           [1x1]    semi_major axis          [km]
% e           [1x1]    eccentricity             [-]
% th1         [1x1]    initial true anomaly     [rad]
% th2         [1x1]    final true anomaly       [rad]
% mu          [1x1]    gravitational parameter  [km^3/s^2] 
% 
%
% Output argument:
% deltat      [1x1]    time of flight           [s]


cosE1 = (e + cos(th1)) / (1 + e * cos(th1));
cosE2 = (e + cos(th2)) / (1 + e * cos(th2));
sinE1 = sqrt(1 - e^2) * sin(th1) / (1 + e * cos(th1));
sinE2 = sqrt(1 - e^2) * sin(th2) / (1 + e * cos(th2));

% Eccentric anomalies
E1 = wrapTo2Pi(atan2(sinE1, cosE1));
E2 = wrapTo2Pi(atan2(sinE2, cosE2));


if th2>th1
    deltat = sqrt(a^3 / mu) * (E2 - E1 - e * (sinE2 - sinE1));
else
    deltat = sqrt(a^3 / mu) * (E2 - E1 - e * (sinE2 - sinE1)) + 2 * pi * sqrt(a^3 / mu);
end
end