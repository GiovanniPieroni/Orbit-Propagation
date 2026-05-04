function [DeltaV1, DeltaV2, Deltat, a_T, e_T] = bitangentTransfer(a_i, e_i, a_f,e_f, type, mu)

% Bitangent transfer for elliptic orbits
% 
% [DeltaV1, DeltaV2, Deltat, a_T, e_T] = bitangentTransfer(a_i,e_i,a_f,e_f,type, mu)
% 
% 
% Input arguments:
% a_i    [1x1]   initial semi-major axis     [km]
% e_i    [1x1]   initial eccentricity        [-]
% a_f    [1x1]   final semi-major axis       [km]
% e_f    [1x1]   final eccentricity          [-]
% type   [char]  maneuver type    
% mu     [1x1]   gravitational parameter     [km^3/s]
% 
% 
% Output arguments:
% DeltaV1    [1x1]   1st maneuver impulse                    [km/s]
% DeltaV2    [1x1]   2nd maneuver impulse                    [km/s]
% Deltat     [1x1]   maneuver time                           [s]
% a_T        [1x1]   semi-major axis of the transfer orbit   [-]
% e_T        [1x1]   eccentricity of the transfer orbit      [-]


switch type
    case 'pa'
        r_pT = a_i * (1 - e_i); % Periapsis radius of the transfer orbit
        r_aT = a_f * (1 + e_f); % Apoapsis radius of the transfer orbit
        a_T = (r_pT + r_aT) / 2;
        DeltaV1 = sqrt(mu) * (sqrt(2 / r_pT - 1 / a_T) - sqrt(2 / r_pT - 1 / a_i));
        DeltaV2 = sqrt(mu) * (sqrt(2 / r_aT - 1 / a_f) - sqrt(2 / r_aT - 1 / a_T));
    case 'ap'
        r_pT = a_f * (1 - e_f);
        r_aT = a_i * (1 + e_i);
        a_T = (r_pT + r_aT) / 2;
        DeltaV1 = sqrt(mu) * (sqrt(2 / r_aT - 1 / a_T) - sqrt(2 / r_aT - 1 / a_i));
        DeltaV2 = sqrt(mu) * (sqrt(2 / r_pT - 1 / a_f) - sqrt(2 / r_pT - 1 / a_T));
    case 'pp'
        r_pT = a_i * (1 - e_i);
        r_aT = a_f * (1 - e_f);
        a_T = (r_pT + r_aT) / 2;
        DeltaV1 = sqrt(mu) * (sqrt(2 / r_pT - 1 / a_T) - sqrt(2 / r_pT - 1 / a_i));
        DeltaV2 = sqrt(mu) * (sqrt(2 / r_aT - 1 / a_f) - sqrt(2 / r_aT - 1 / a_T));
    case 'aa'
        r_pT = a_i * (1 + e_i);
        r_aT = a_f * (1 + e_f);
        a_T = (r_pT + r_aT) / 2;
        DeltaV1 = sqrt(mu) * (sqrt(2 / r_pT - 1 / a_T) - sqrt(2 / r_pT - 1 / a_i));
        DeltaV2 = sqrt(mu) * (sqrt(2 / r_aT - 1 / a_f) - sqrt(2 / r_aT - 1 / a_T));
    otherwise
         warning('Inserire una categoria corretta di trasferimento')
end

DeltaV1 = abs(DeltaV1);
DeltaV2 = abs(DeltaV2);
e_T = (r_aT - r_pT) / (r_aT + r_pT);
Deltat = pi * sqrt(a_T^3 / mu);


end