function [DeltaV1, DeltaV2, Deltat, a_T, e_T] = bitangentTransferFun_ai(e_i, a_f, e_f, type, mu)
% Bitangent transfer for elliptic orbits - function of initial semi-major
% axis a_i
% 
% [DeltaV1, DeltaV2, Deltat] = bitangentTransfer(e_i, a_f, e_f, type, mu)
% 
% 
% Input arguments:
% e_i    [1x1]   initial eccentricity       [-]
% a_f    [1x1]   final semi-major axis      [km]
% e_f    [1x1]   final eccentricity         [-]
% type   [char]  maneuver type              [-]
% mu     [1x1]   gravitational parameter    [km^3./s]
% 
% 
% 
% Output arguments:
% DeltaV1    [1x1]   1st maneuver impulse       [km./s]
% DeltaV2    [1x1]   2nd maneuver impulse       [km./s]
% Deltat     [1x1]   maneuver time              [s]
% a_t        [1x1]   maneuver semi-major axis   [km]
% e_t        [1x1]   maneuver eccentricity      [-]


switch type
    case 'pa'
        r_pT = @(a_i) a_i .* (1 - e_i);
        r_aT =  a_f .* (1 + e_f);
        a_T = @(a_i) (r_pT(a_i) + r_aT) ./ 2;

        DeltaV1 = @(a_i) abs(sqrt(mu) .* (sqrt(2 ./ r_pT(a_i) - 1 ./ a_T(a_i)) - sqrt(2 ./ r_pT(a_i) - 1 ./ a_i)));
        DeltaV2 = @(a_i) abs(sqrt(mu) .* (sqrt(2 ./ r_aT - 1 ./ a_f) - sqrt(2 ./ r_aT - 1 ./ a_T(a_i))));

        e_T = @(a_i) (r_aT - r_pT(a_i)) ./ (r_aT + r_pT(a_i));
        Deltat = @(a_i) pi .* sqrt(a_T(a_i).^3 ./ mu);

    case 'ap'
        r_pT = a_f .* (1 - e_f);
        r_aT = @(a_i) a_i .* (1 + e_i);
        a_T = @(a_i) (r_pT + r_aT(a_i)) ./ 2;

        DeltaV1 = @(a_i) abs(sqrt(mu) .* (sqrt(2 ./ r_aT(a_i) - 1 ./ a_T(a_i)) - sqrt(2 ./ r_aT(a_i) - 1 ./ a_i)));
        DeltaV2 = @(a_i) abs(sqrt(mu) .* (sqrt(2 ./ r_pT - 1 ./ a_f) - sqrt(2 ./ r_pT - 1 ./ a_T(a_i))));

        e_T = @(a_i) (r_aT(a_i) - r_pT) ./ (r_aT(a_i) + r_pT);
        Deltat = @(a_i) pi .* sqrt(a_T(a_i).^3 ./ mu);

    case 'pp'
        r_pT = @(a_i) a_i .* (1 - e_i);
        r_aT = a_f .* (1 - e_f);
        a_T = @(a_i) (r_pT(a_i) + r_aT) ./ 2;

        DeltaV1 = @(a_i) abs(sqrt(mu) .* (sqrt(2 ./ r_pT(a_i) - 1 ./ a_T(a_i)) - sqrt(2 ./ r_pT(a_i) - 1 ./ a_i)));
        DeltaV2 = @(a_i) abs(sqrt(mu) .* (sqrt(2 ./ r_aT - 1 ./ a_f) - sqrt(2 ./ r_aT - 1 ./ a_T(a_i))));

        e_T = @(a_i) (r_aT - r_pT(a_i)) ./ (r_aT + r_pT(a_i));
        Deltat = @(a_i) pi .* sqrt(a_T(a_i).^3 ./ mu);

    case 'aa'
        r_pT = @(a_i) a_i .* (1 + e_i);
        r_aT = a_f .* (1 + e_f);
        a_T = @(a_i) (r_pT(a_i) + r_aT) ./ 2;

        DeltaV1 = @(a_i) abs(sqrt(mu) .* (sqrt(2 ./ r_pT(a_i) - 1 ./ a_T(a_i)) - sqrt(2 ./ r_pT(a_i) - 1 ./ a_i)));
        DeltaV2 = @(a_i) abs(sqrt(mu) .* (sqrt(2 ./ r_aT - 1 ./ a_f) - sqrt(2 ./ r_aT - 1 ./ a_T(a_i))));

        e_T = @(a_i) (r_aT - r_pT(a_i)) ./ (r_aT + r_pT(a_i));
        Deltat = @(a_i) pi .* sqrt(a_T(a_i).^3 ./ mu);

    otherwise
        warning('Inserire una categoria corretta di trasferimento')
end


end