function [DeltaV1, DeltaV2, DeltaV3, Deltat1, Deltat2, a_T1, e_T1, a_T2, e_T2] = biellipticTransfer(ai, ei, af, ef, r_aT, mu)

% Bitangent transfer for elliptic orbits
%
% [DeltaV1,DeltaV2,DeltaV3,Deltat1,Deltat2] = biellipticTransfer(ai,ei,af,ef,r_aT, mu)
% 
%
% Input arguments:
% ai    [1x1]   initial semi-major axis              [km]
% ei    [1x1]   initial eccentricity                 [-]
% af    [1x1]   final semi-major axis                [km]
% ef    [1x1]   final eccentricity                   [-]
% r_aT  [1x1]   transfer orbits apocenter distance   [km]
% mu    [1x1]   gravitational parameter              [km^3/s^2]
%
%
% Output arguments:
% DeltaV1   [1x1]   1st maneuver impulse                         [km/s] 
% DeltaV2   [1x1]   2nd maneuver impulse                         [km/s] 
% DeltaV3   [1x1]   3rd maneuver impulse                         [km/s]
% Deltat1   [1x1]   maneuver time 1                              [s]
% Deltat2   [1x1]   maneuver time 2                              [s]
% a_T1      [1x1]   initial semi-major axis                      [km]
% e_T1      [1x1]   eccentricity of the first transfer orbit     [-]
% a_T2      [1x1]   initial semi-major axis                      [km]
% e_T2      [1x1]   eccentricity of the second transfer orbit    [-]

r_pT1 = ai * (1 - ei); % Periapsis radius of the first transfer orbit
r_pT2 = af * (1 - ef); % Periapsis radius of the second transfer orbit

a_T1 = (r_pT1 + r_aT) / 2;
a_T2 = (r_pT2 + r_aT) / 2;

e_T1 = r_aT / a_T1 - 1;
e_T2 = r_aT / a_T2 - 1;

DeltaV1 = abs(sqrt(mu) * (sqrt(2 / r_pT1 - 1 / a_T1) - sqrt(2 / r_pT1 - 1 / ai)));
DeltaV2 = abs(sqrt(mu) * (sqrt(2 / r_aT - 1/ a_T2) - sqrt(2 / r_aT - 1 / a_T1)));
DeltaV3 = abs(sqrt(mu) * (sqrt(2 / r_pT2 - 1 / af) - sqrt(2 / r_pT2 - 1 / a_T2)));

Deltat1 = pi * sqrt(a_T1^3 / mu);
Deltat2 = pi * sqrt(a_T2^3 / mu);



end