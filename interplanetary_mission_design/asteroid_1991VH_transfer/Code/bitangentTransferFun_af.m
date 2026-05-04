function [DeltaV1, DeltaV2, Deltat, a_T, e_T] = bitangentTransferFun_af(a_i, e_i, e_f, type, mu)
% Bitangent transfer for elliptic orbits - function of final semi-major
% axis a_f
% 
% [DeltaV1, DeltaV2, Deltat] = bitangentTransfer(a_i, e_i, e_f, type, mu)
% 
% 
% Input arguments:
% a_i    [1x1]   initial semi-major axis    [km]
% e_i    [1x1]   initial eccentricity       [-]
% e_f    [1x1]   final eccentricity         [-]
% type   [char]  maneuver type              [-]
% mu     [1x1]   gravitational parameter    [km^3./s]
% 
% 
% 
% Output arguments - function handle:
% DeltaV1    [1x1]   1st maneuver impulse       [km./s]
% DeltaV2    [1x1]   2nd maneuver impulse       [km./s]
% Deltat     [1x1]   maneuver time              [s]
% a_t        [1x1]   maneuver semi-major axis   [km]
% e_t        [1x1]   maneuver eccentricity      [-]

switch type
    case 'pa'
        r_pT = a_i .* (1 - e_i);
        r_aT = @(a_f) a_f .* (1 + e_f);
        a_T = @(a_f) (r_pT + r_aT(a_f)) ./ 2;

        DeltaV1 = @(a_f) sqrt(mu) .* ( sqrt(2 ./ r_pT - 1 ./ a_T(a_f)) - sqrt(2 ./ r_pT - 1 ./ a_i) );
        DeltaV2 = @(a_f) sqrt(mu) .* ( sqrt(2 ./ r_aT(a_f) - 1 ./ a_f) - sqrt(2 ./ r_aT(a_f) - 1 ./ a_T(a_f)));

    case 'ap'
        r_pT = @(a_f) a_f .* (1 - e_f);
        r_aT = a_i .* (1 + e_i);
        a_T = @(a_f) (r_pT(a_f) + r_aT) ./ 2;

        DeltaV1 = @(a_f) sqrt(mu) .* ( sqrt(2 ./ r_aT - 1 ./ a_T(a_f)) - sqrt(2 ./ r_aT - 1 ./ a_i) );
        DeltaV2 = @(a_f) sqrt(mu) .* ( sqrt(2 ./ r_pT(a_f) - 1 ./ a_f) - sqrt(2 ./ r_pT(a_f) - 1 ./ a_T(a_f)));

    case 'pp'
        r_pT = a_i .* (1 - e_i);
        r_aT = @(a_f) a_f .* (1 - e_f);
        a_T = @(a_f) (r_pT + r_aT(a_f)) ./ 2;

        DeltaV1 = @(a_f) sqrt(mu) .* ( sqrt(2 ./ r_pT - 1 ./ a_T(a_f)) - sqrt(2 ./ r_pT - 1 ./ a_i) );
        DeltaV2 = @(a_f) sqrt(mu) .* ( sqrt(2 ./ r_aT(a_f) - 1 ./ a_f) - sqrt(2 ./ r_aT(a_f) - 1 ./ a_T(a_f)));

    case 'aa'
        r_pT = a_i .* (1 + e_i);
        r_aT = @(a_f) a_f .* (1 + e_f);
        a_T = @(a_f) (r_pT + r_aT(a_f)) ./ 2;

        DeltaV1 = @(a_f) sqrt(mu) .* ( sqrt(2 ./ r_pT - 1 ./ a_T(a_f)) - sqrt(2 ./ r_pT - 1 ./ a_i) );
        DeltaV2 = @(a_f) sqrt(mu) .* ( sqrt(2 ./ r_aT(a_f) - 1 ./ a_f) - sqrt(2 ./ r_aT(a_f) - 1 ./ a_T(a_f)));

    otherwise
        warning('Inserire una categoria corretta di trasferimento')
end

DeltaV1 = @(a_f) abs(DeltaV1(a_f));
DeltaV2 = @(a_f) abs(DeltaV2(a_f));
e_T = @(a_f) (r_aT(a_f) - r_pT) ./ (r_aT(a_f) + r_pT);
Deltat = @(a_f) pi .* sqrt(a_T(a_f)^3 ./ mu);


end