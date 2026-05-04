function Deltat = customTOF3(th_i, th_f, om, a_T, e_T, th1, th2, mu)
   
% Time of flight as a function of the initial true anomaly, the final true 
% anomaly and the argument of periapsis of the transfer orbit during a
% direct transfer between two orbits that do not intersect
% 
% Deltat = customTOF3(th_i, th_f, om, a, e, th1, th2, mu)
% 
% 
% Input arguments:
% th_i [1x1]   initial true anomaly                         [rad]
% th_f [1x1]   final true anomaly                           [rad]
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
    e_val = e_T(th_i, th_f, om); % Compute the value of eccentricity
        if ~isreal(e_val) || e_val <= 0 || e_val >= 1 % Consider only closed orbits
            Deltat = NaN;
            return;
        end
        a_val = a_T(th_i, th_f, om); % Compute the value of the semi-major axis
        th1_val = th1(th_i, th_f, om); % Compute the value of the initial 
                                       % true anomaly within transfer orbit
        th2_val = th2(th_i, th_f, om); % Compute the value of the final 
                                       % true anomaly within transfer orbit
        if any(~isreal([a_val, th1_val, th2_val])) % Consider only real values
            Deltat = NaN;
            return;
        end
        Deltat = TOF(a_val, e_val, th1_val, th2_val, mu); % Compute time of flight
    catch
        Deltat = NaN;
end
end