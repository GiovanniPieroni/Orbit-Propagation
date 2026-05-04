function Deltat = customTOF(om, a, e, th1, th2, mu)

% Time of flight as a function of the argument of periapsis of the orbit
% 
% Deltat = customTOF(om, a, e, th1, th2, mu)
% 
% 
% Input arguments:
% om   [1x1]   argument of periapsis                        [rad]
% a    [1x1]   semi-major axis                              [km]
% e    [1x1]   eccentricity                                 [-]
% th1  [1x1]   initial true anomaly within transfer orbit   [rad]
% th2  [1x1]   final true anomaly within transfer orbit     [rad]
% mu   [1x1]   gravitational parameter                      [km^3/s]
% 
% 
% Output argument:
% Deltat  [1x1]  time of flight   [s]

try
    e_val = e(om); % Compute the value of eccentricity
    if ~isreal(e_val) || e_val <= 0 || e_val >= 1 % Consider only closed orbits
        Deltat = NaN; 
        return;
    end
    a_val = a(om); % Compute the value of the semi-major axis
    th1_val = th1(om); % Compute the value of the initial true anomaly within 
                       % the transfer orbit
    th2_val = th2(om); % Compute the value of the final true anomaly within 
                       % the transfer orbit
        if any(~isreal([a_val, th1_val, th2_val])) % Consider only real values
            Deltat = NaN;
            return;
        end
        Deltat = TOF(a_val, e_val, th1_val, th2_val, mu);
catch
    Deltat = NaN;
end
end