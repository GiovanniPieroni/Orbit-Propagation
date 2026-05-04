function [DeltaV, omf, th] = changeOrbitalPlane(a, e, i_i, OMi, omi, i_f, OMf, mu)

% Change of Plane maneuver
%
% [DeltaV, omf, th] = changeOrbitalPlane(a, e, i_i, OMi, omi, i_f, OMf, mu)
%
%
% Input arguments:
% a         [1x1]   semi-major axis                 [km]
% e         [1x1]   eccentricity                    [-]
% i_i       [1x1]   initial inclination             [rad]
% OMi       [1x1]   initial RAAN                    [rad]
% omi       [1x1]   initial argument of periapsis   [rad]
% i_f       [1x1]   final inclination               [rad] 
% OMf       [1x1]   final RAAN                      [rad]
% mu        [1x1]   gravitational parameter         [km^3/s^2]
%
% Output arguments:
% DeltaV    [1x1]   maneuver impulse                [km/s]
% omf       [1x1]   final argument of periapsis     [rad]
% th        [1x1]   true anomaly at maneuver        [rad]



DeltaOM = OMf - OMi; % Variation of RAAN
Deltai = i_f - i_i; % Variation of inclination
alpha = wrapTo2Pi(acos(cos(i_i) * cos(i_f) + sin(i_i) * sin(i_f) * cos(DeltaOM)));

if DeltaOM>0 && Deltai>0
   cos_ui = (- cos(i_f) + cos(alpha) * cos(i_i)) / (sin(alpha) * sin(i_i));
   cos_uf = (cos(i_i) - cos(alpha) * cos(i_f)) / (sin(alpha) * sin(i_f));
   sin_ui = sin(DeltaOM) / sin(alpha) * sin(i_f);
   sin_uf = sin(DeltaOM) / sin(alpha) * sin(i_i);
   ui = wrapTo2Pi(atan2(sin_ui, cos_ui));
   uf = wrapTo2Pi(atan2(sin_uf, cos_uf));
   th = wrapTo2Pi(ui - omi);
   omf = wrapTo2Pi(uf - th);
elseif DeltaOM>0 && Deltai<0
   cos_ui = (cos(i_f) - cos(alpha) * cos(i_i)) / (sin(alpha) * sin(i_i));
   cos_uf = (-cos(i_i) + cos(alpha) * cos(i_f)) / (sin(alpha) * sin(i_f));
   sin_ui = sin(DeltaOM) / sin(alpha) * sin(i_f);
   sin_uf = sin(DeltaOM) / sin(alpha) * sin(i_i);
   ui = wrapTo2Pi(atan2(sin_ui, cos_ui));
   uf = wrapTo2Pi(atan2(sin_uf, cos_uf));
   th = wrapTo2Pi(2 * pi - ui - omi);
   omf = wrapTo2Pi(2 * pi - uf - th);
elseif DeltaOM<0 && Deltai>0
   DeltaOM = abs(DeltaOM);
   cos_ui = (- cos(i_f) + cos(alpha) * cos(i_i)) / (sin(alpha) * sin(i_i));
   cos_uf = (cos(i_i) - cos(alpha) * cos(i_f)) / (sin(alpha) * sin(i_f));
   sin_ui = sin(DeltaOM) / sin(alpha) * sin(i_f);%Controllare abs(DeltaOM)
   sin_uf = sin(DeltaOM) / sin(alpha) * sin(i_i);%Controllare abs(DeltaOM)
   ui = wrapTo2Pi(atan2(sin_ui, cos_ui));
   uf = wrapTo2Pi(atan2(sin_uf, cos_uf));
   th = wrapTo2Pi(2 * pi - ui - omi);
   omf = wrapTo2Pi(2 * pi - uf - th);
elseif DeltaOM<0 && Deltai<0
   DeltaOM = abs(DeltaOM);
   cos_ui = (cos(i_f) - cos(alpha) * cos(i_i)) / (sin(alpha) * sin(i_i));
   cos_uf = (-cos(i_i) + cos(alpha) * cos(i_f)) / (sin(alpha) * sin(i_f));
   sin_ui = sin(DeltaOM) / sin(alpha) * sin(i_f);%Controllare abs(DeltaOM)
   sin_uf = sin(DeltaOM) / sin(alpha) * sin(i_i);%Controllare abs(DeltaOM)
   ui = wrapTo2Pi(atan2(sin_ui, cos_ui));
   uf = wrapTo2Pi(atan2(sin_uf, cos_uf));
   th = wrapTo2Pi(ui - omi);
   omf = wrapTo2Pi(uf - th);
end


if cos(th)>0
   th = wrapTo2Pi(th + pi); %To minimize DeltaV
end

p = a * (1 - e^2); % Semi-latus rectum
v_t = sqrt(mu / p) * (1 + e * cos(th)); % Transverse velocity
DeltaV = 2 * v_t * sin(alpha / 2);


end