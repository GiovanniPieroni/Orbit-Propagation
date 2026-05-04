%%%%%%%%%%%%%%%%%%%%%%%% General Setup - Run first %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                                                   %%%%%%%%%%%%
%%%%%%%%%%%%                    GROUP 22                       %%%%%%%%%%%%
%%%%%%%%%%%%                                                   %%%%%%%%%%%%
%%%%%%%%%%%%            Giovanni Pieroni - 10854347            %%%%%%%%%%%%
%%%%%%%%%%%%            Alessandro Ponti - 10901276            %%%%%%%%%%%%
%%%%%%%%%%%%            Jacopo Rotta     - 10696753            %%%%%%%%%%%%
%%%%%%%%%%%%    SC02 - Interplanetary trajectory to asteroid   %%%%%%%%%%%%
%%%%%%%%%%%%               35107 (1991 VH)                     %%%%%%%%%%%%
%%%%%%%%%%%%                                                   %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
clc

format short g

% Standard gravitational parameter of the sun
mu_s = 132712440018;   % [km^3/s^2]
AU = 149597870.7;      % [km]

% Initial and final parameters definition
a_i = 1.00000011 * AU; % [km]
e_i = 0.01671022;      % [-]
i_i = 9.1920e-5;       % [rad]
OM_i = 2.7847;         % [rad]
om_i = 5.2643;         % [rad]


a_f = 1.137342 * AU;   % [km]
e_f = 0.144236;        % [-]
i_f = deg2rad(13.91);  % [rad]
OM_f = deg2rad(139.33);% [rad]
om_f = deg2rad(206.94);% [rad]

% Unit vector I of the X-axis
I = [1 0 0];
% Unit vector J of the Y-axis
J = [0 1 0];
% Unit vector K of the Z-axis
K = [0 0 1];


I_s = 300 * 9.807;  % Specific gravimetric impulse
m_0 = 300;          % Empty satellite's mass

% Maximum possible DeltaV with 300 s of specific impulse and mass ratio 25
DeltaV_MAX = I_s * log(25) / 1000;  %Tsiolkovsky
% Maximum possible time of transfer of 240 days
DeltaT_MAX = 240 * 24 * 60 * 60;


%% Min Delta V - vectorized
clearvars -except AU mu_s a_i e_i i_i OM_i om_i th_i a_f e_f i_f OM_f om_f...
                  I_s m_0 DeltaV_MAX DeltaT_MAX I J K




% Initial guess
th_i_ott = 0;
th_f_ott = 0;
om_t_ott = 0;

% Indexes
jj = 1; % Discretization resolution
ii = 2; % While cicle index
N_max = 50; % Max number of iterations
tol = 1e-6; % Error tolerance between while cicle steps (DeltaV)
% Initial values for DeltaV - just to start the cicle
DV_cicle = [100, 90; 100, 90; 100, 90];
DV_cicle = [DV_cicle, zeros(3, N_max)];

flag = 1; % Exit flag at convergence

% Loop until tolerance is obtained
while  (abs(DV_cicle(1, ii-1) - DV_cicle(1, ii)) > tol || abs(DV_cicle(2, ii-1)...
        - DV_cicle(2, ii)) > tol || abs(DV_cicle(3, ii-1) - DV_cicle(3, ii)) > tol) ...
        && flag ~= 2 && ii < N_max

    % Discretizing th_i, th_f and om_t
    N = 50; % Discretization value
    th_i = wrapTo2Pi(th_i_ott - pi/jj: pi/(N*jj) : th_i_ott + pi/jj); % True anomaly at the initial orbit
    th_f = wrapTo2Pi(th_f_ott - pi/jj: pi/(N*jj) : th_f_ott + pi/jj); % NOTE: length(th_i) == length(th_f) == om_t
    om_t = wrapTo2Pi(om_t_ott - pi/jj: pi/(N*jj) : om_t_ott + pi/jj); % Argument of periapsis of transfer orbit
    dim = size(th_i, 2);    % Vectors dimension

    % Initializing vectors
    rr_i = zeros(3, dim);   % Initial ECI position
    vv_i = zeros(3, dim);   % Initial ECI velocity
    rr_f = zeros(3, dim);   % Final ECI position
    vv_f = zeros(3,dim);    % Final ECI velocity

    R_om_comb = zeros(3, 1);% Rotation matrix of om_t
    DeltaV_ott = 10;        % Initial optimal DeltaV

    % Combining initial and final vectors rr & vv; creating rotation matrix of om_t
    for idx = 1 : dim
        % Initial position and velocity
        [rr_i(:,idx), vv_i(:,idx)] = par2car(a_i, e_i, i_i, OM_i, om_i, th_i(idx), mu_s);
        rv_i = [rr_i; vv_i];
        % Final position and velocity
        [rr_f(:,idx),vv_f(:,idx)] = par2car(a_f, e_f, i_f, OM_f, om_f, th_f(idx), mu_s);
        rv_f = [rr_f; vv_f];

        % Construction of rotation matrix for every om_t
        R_om = [cos(om_t(idx)), sin(om_t(idx)), 0;
            -sin(om_t(idx)), cos(om_t(idx)), 0;
            0, 0, 1];
        R_om_comb(:, end+1:end+3) = R_om;
    end

    % Breaking down the rotation matrix from 2 to 3 dimensions
    R_om_comb(:,1) = [];
    R_om_blocks = reshape(R_om_comb, 3, 3, []);

    % Calculating every combination of rr_i and rr_f:
    comb = combinations(1:dim, 1:dim); % Every pair between 1 and dim
    comb = table2array(comb);
    comb_i = comb(:,1)';
    comb_f = comb(:,2)';

    rr_comb = [rr_i(:,comb_i); rr_f(:,comb_f)];  % [rr_i, rr_f]
    vv_comb = [vv_i(:,comb_i); vv_f(:,comb_f)];  % [vv_i, vv_f]

    % Reshaping vectors
    rr_i_blocks = reshape(rr_comb(1:3,:), 3, 1, []);
    rr_f_blocks = reshape(rr_comb(4:6,:), 3, 1, []);
    vv_i_blocks = reshape(vv_comb(1:3,:), 3, 1, []);
    vv_f_blocks = reshape(vv_comb(4:6,:), 3, 1, []);


    % Transfer orbit parameters and Delta V calculation:
    for j = 1 : size(comb,1)  % For every combination

        % Extracting one combination of initial and final ECI vectors
        rr_i_current = rr_i_blocks(:,:,j);
        rr_f_current = rr_f_blocks(:,:,j);
        vv_i_current = vv_i_blocks(:,:,j);
        vv_f_current = vv_f_blocks(:,:,j);

        rr_i_norm = norm(rr_i_current);
        rr_f_norm = norm(rr_f_current);

        % Transfer orbital plane direction
        h_t = cross(rr_i_current, rr_f_current)./norm(cross(rr_i_current, rr_f_current));

        % Transfer orbit inclination
        i_t = acos(h_t(3));

        % Transfer orbit node line
        N_T = cross(K, h_t)./norm(cross(K, h_t));

        % RAAN
        if N_T(2) >= 0
            OM_t = acos(N_T(1));
        else
            OM_t = 2*pi - acos(N_T(1));
        end

        % Rotation of OM around K
        R_OM = [cos(OM_t), sin(OM_t), 0;
            -sin(OM_t), cos(OM_t), 0;
            0, 0, 1];

        % Rotation of i around N
        R_i = [1,      0,       0;
            0, cos(i_t), sin(i_t);
            0, -sin(i_t), cos(i_t)];

        % Creating every possible rotation matrix
        R_comb = pagemtimes(R_om_blocks,R_i*R_OM);

        % Position vectors on PF frame, for every om_t
        rr_i_pf_current = pagemtimes(R_comb, rr_i_current);
        rr_f_pf_current = pagemtimes(R_comb, rr_f_current);

        % True anomaly of the initial transfer position
        cos_th_1t = rr_i_pf_current(1,:,:)./rr_i_norm;
        sin_th_1t = rr_i_pf_current(2,:,:)./rr_i_norm;
        th_1t = wrapTo2Pi(atan2(sin_th_1t(:), cos_th_1t(:)))';

        % True anomaly of the final transfer position
        cos_th_2t = rr_f_pf_current(1,:,:)./rr_f_norm;
        sin_th_2t = rr_f_pf_current(2,:,:)./rr_f_norm;
        th_2t = wrapTo2Pi(atan2(sin_th_2t(:), cos_th_2t(:)))';

        % Eccentricity modulus of current transfer orbit
        e_t_current =  (rr_i_norm - rr_f_norm)./(rr_f_norm .*  cos(th_2t) -rr_i_norm .* cos(th_1t));
        % Discarding non acceptable negative eccentricity
        e_t_current(e_t_current < 0) = NaN;

        % Semi-major axis of current transfer orbit
        a_t_current = (rr_i_norm.*(1+e_t_current.*cos(th_1t)))./(1-(e_t_current.^2));
        % Discarding non acceptable semi-major axis (minimum distance PSP 6.856.340 km)
        a_t_current(a_t_current < 1e7) = NaN;


        % Calculation of vv_1t and vv_2t
        p = a_t_current.*(1-(e_t_current.^2));   % Current semi-latus rectum

        vv_1t_PF = sqrt(mu_s./p) .* [-sin(th_1t);
            e_t_current + cos(th_1t);
            zeros(1, length(th_1t))];
        vv_1t_PF = reshape(vv_1t_PF,3,1,[]);
        vv_2t_PF = sqrt(mu_s./p) .* [-sin(th_2t);
            e_t_current + cos(th_2t);
            zeros(1, length(th_2t))];
        vv_2t_PF = reshape(vv_2t_PF,3,1,[]);

        % Rotation matrix from PF to ECI
        R_PFtoECI = pagetranspose(R_comb);

        vv_1t = pagemtimes(R_PFtoECI, vv_1t_PF);
        vv_2t_PF = pagemtimes(R_PFtoECI, vv_2t_PF);


        % Delta V calculation
        DeltaV1 = vecnorm(vv_1t - vv_i_current, 2, 1); % First Delta V
        DeltaV2 = vecnorm(vv_f_current - vv_2t_PF, 2, 1); % Second Delta V
        DeltaV = DeltaV1(:) + DeltaV2(:);
        % if all(isreal(DeltaV))==0
        %     error('Velocità immaginarie');
        % end
        [DeltaV_candidate, index] = min(DeltaV);

        % Extracting current optimum
        if DeltaV_candidate < DeltaV_ott

            % Optimal Delta V
            DeltaV_ott = DeltaV_candidate;
            DeltaVV_1 = vv_1t - vv_i_current;
            DeltaVV_1 = DeltaVV_1(:, :, index);
            DeltaVV_2 = vv_f_current - vv_2t_PF;
            DeltaVV_2 = DeltaVV_2(:, :, index);

            % Orbital parameters of optimal orbit
            i_t_ott = i_t;
            OM_t_ott = OM_t;
            e_t_ott = e_t_current(index);
            a_t_ott = a_t_current(index);
            th_1t_ott = th_1t(index);
            th_2t_ott = th_2t(index);
            om_t_ott = om_t(index);

            if (rem(j,length(th_1t))) == 0
                th_i_ott = th_i(ceil(j/length(th_1t)));
                th_f_ott = th_f(index);
                flag = 2;
                break
            else      
                th_i_ott = th_i(ceil(j/length(th_1t)));
                th_f_ott = th_f(rem(j,length(th_1t)));
            end
            
        end
    end

    DV_cicle(:, ii) = DeltaVV_1;
    ii = ii + 1;
    jj = jj * ii/2;
end

% Time of transfer
DeltaT = TOF(a_t_ott, e_t_ott, wrapTo2Pi(th_1t_ott), wrapTo2Pi(th_2t_ott), mu_s);

% Propellant mass calculation
m_p = m_0 * (exp(DeltaV_ott * 1e3/I_s) - 1);

% Displaing results
fprintf("================================================================\n"  + ...
    "Optimal  ΔV =  %.6f km/s (%.1f%% of  ΔV_MAX)\n" + ...
    "Respective  Δt = %.4f s = %.4f d \n" + ...
    "Propellant mass (with Is = 300 s & M_empty = 300 kg): %.2f kg\n" + ...
    "================================================================\n\n" ...
    , DeltaV_ott, 100*DeltaV_ott/DeltaV_MAX, DeltaT, DeltaT/(3600*24), m_p);



% Table A6_2
A6_2(1,:) = [0, a_i/AU, e_i, i_i, OM_i, om_i, th_i_ott, norm(DeltaVV_1)];
A6_2(2,:) = [0, a_t_ott/AU, e_t_ott, i_t_ott, OM_t_ott, om_t_ott, th_1t_ott, 0];
A6_2(3,:) = [DeltaT/3600, a_t_ott/AU, e_t_ott, i_t_ott, OM_t_ott, om_t_ott, th_2t_ott, norm(DeltaVV_2)];
A6_2(4,:) = [DeltaT/3600, a_f/AU, e_f, i_f, OM_f, om_f, th_f_ott, 0];

TAB6_2 = array2table(A6_2, 'VariableNames', ...
    {'Time [h]', 'a [AU]', 'e [-]', 'i [rad]', 'OM [rad]', 'om [rad]', 'theta [rad]', 'DeltaV [km/s]'});
disp(TAB6_2)


% Plotting the optimal transfer
figure('Name', 'A6 - Min Delta vectorized', 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

% Plot the Sun
sun_sphere
hold on
grid on
box on
axis equal

% Initial orbit and position(Earth)
plotOrbit(a_i, e_i, i_i, OM_i, om_i, 0, 2*pi, pi/10000, mu_s);
rr_start = par2car(a_i, e_i, i_i, OM_i, om_i, th_i_ott, mu_s);
scatter3(rr_start(1), rr_start(2), rr_start(3),'*', 'LineWidth', 4)

% Final orbit and position (Asteroid)
plotOrbit(a_f, e_f, i_f, OM_f, om_f, 0, 2*pi, pi/10000, mu_s);
rr_finish = par2car(a_f, e_f, i_f, OM_f, om_f, th_f_ott,mu_s);
scatter3(rr_finish(1), rr_finish(2), rr_finish(3),'*', 'LineWidth', 4)

% Transfer orbit
plotOrbit(a_t_ott, e_t_ott, i_t_ott, OM_t_ott, om_t_ott, 0, 2*pi, pi/1000, mu_s);
plotOrbit(a_t_ott, e_t_ott, i_t_ott, OM_t_ott, om_t_ott, th_1t_ott, th_2t_ott, pi/1000, mu_s, 3);

leg = legend('','Inital orbit','Inital point (1st impulse)','Final orbit','Final point (2nd impulse)','','Transfer','Interpreter','latex');
set(leg,'FontSize',22);
title('Direct transfer between the two orbits with $\Delta V_{MIN}$', 'Interpreter', 'latex')

%% fmincon setup - Run before the following sections
clearvars -except AU mu_s a_i e_i i_i OM_i om_i th_i a_f e_f i_f OM_f om_f I_s m_0 DeltaV_MAX DeltaT_MAX I J K



% % %  % Everything (if needed) is function of th_i, th_f and om_t % % % %
% Initial position and velocity
[rr_i, vv_i] = par2carFUN1th(a_i, e_i, i_i, OM_i, om_i, mu_s);

% Final position and velocity
[rr_f, vv_f] = par2carFUN1th(a_f, e_f, i_f, OM_f, om_f, mu_s);


% Direction of the angular momentum vector of the transfer orbit
h_T = @(th_i, th_f) cross(rr_i(th_i), rr_f(th_f))/norm(cross(rr_i(th_i), rr_f(th_f)));

% Inclination of the transfer orbit
i_T = @(th_i, th_f) acos(dot(h_T(th_i,th_f),K));

% Line of nodes
N_T = @(th_i, th_f) cross(K, h_T(th_i, th_f))/norm(cross(K, h_T(th_i, th_f)));

% Right ascension of the ascending node
OM_T = @(th_i, th_f) wrapTo2Pi(atan2(dot(N_T(th_i,th_f), J), dot(N_T(th_i,th_f), I)));

% Rotation matrices for transfer orbit frame transformation
R_OM_T = @(th_i, th_f) [cos(OM_T(th_i,th_f)), sin(OM_T(th_i,th_f)),   0;
    -sin(OM_T(th_i,th_f)), cos(OM_T(th_i,th_f)),   0;
    0,                 0,               1];

R_i_T = @(th_i, th_f) [1,            0,                   0;
    0,    cos(i_T(th_i,th_f)), sin(i_T(th_i,th_f));
    0,   -sin(i_T(th_i,th_f)), cos(i_T(th_i,th_f))];

R_om_T = @(om_T) [cos(om_T), sin(om_T), 0;
    -sin(om_T), cos(om_T), 0;
    0,         0,      1];

T_EliotoPF = @(th_i, th_f, om_T) R_om_T(om_T) * R_i_T(th_i, th_f) * R_OM_T(th_i, th_f);

rr_i_PF = @(th_i, th_f, om_T) T_EliotoPF(th_i, th_f, om_T) * rr_i(th_i);
rr_f_PF = @(th_i, th_f, om_T) T_EliotoPF(th_i, th_f, om_T) * rr_f(th_f);

% Cosine and sine of true anomaly at departure (in PF frame)
CosTh1_T = @(th_i,th_f,om_T) dot(rr_i_PF(th_i,th_f,om_T), I) / norm(rr_i(th_i));
SinTh1_T = @(th_i,th_f,om_T) dot(rr_i_PF(th_i,th_f,om_T), J) / norm(rr_i(th_i));

% Cosine and sine of true anomaly at arrival (in PF frame)
CosTh2_T = @(th_i,th_f,om_T) dot(rr_f_PF(th_i,th_f,om_T), I) / norm(rr_f(th_f));
SinTh2_T = @(th_i,th_f,om_T) dot(rr_f_PF(th_i,th_f,om_T), J) / norm(rr_f(th_f));

% True anomalies on the transfer orbit
Th1_T = @(th_i,th_f,om_T) wrapTo2Pi(atan2(SinTh1_T(th_i,th_f,om_T), CosTh1_T(th_i,th_f,om_T)));
Th2_T = @(th_i,th_f,om_T) wrapTo2Pi(atan2(SinTh2_T(th_i,th_f,om_T), CosTh2_T(th_i,th_f,om_T)));

% Eccentricity of the transfer orbit
e_T = @(th_i, th_f, om_T) (norm(rr_f(th_f)) - norm(rr_i(th_i))) ...
    / (CosTh1_T(th_i,th_f,om_T) * norm(rr_i(th_i)) -...
    CosTh2_T(th_i,th_f,om_T) * norm(rr_f(th_f)));

% Semi-major axis of the transfer orbit
a_T = @(th_i, th_f, om_T) norm(rr_i(th_i)) * (1 + e_T(th_i,th_f,om_T) * ...
    CosTh1_T(th_i,th_f,om_T)) / (1 - (e_T(th_i,th_f,om_T))^2);

% HCI positions at the initial and final true anomalies on the transfer
% orbit
[rr1_T, vv1_T] = par2carFUN3(a_T, e_T, i_T, OM_T, Th1_T, mu_s);
[rr2_T, vv2_T] = par2carFUN3(a_T, e_T, i_T, OM_T, Th2_T, mu_s);

% Maneuver cost as a function of the true anomalies of the departure and
% arrival positions, and the argument of periapsis of the transfer orbit
DeltaV_fun = @(x) norm(vv1_T(x(1), x(2), x(3)) - vv_i(x(1))) + ...
    norm(vv_f(x(2)) - vv2_T(x(1), x(2), x(3)));


% Time of flight for the direct transfer as a function of the true anomalies
% of the initial and final positions and the argument of periapsis
% of the transfer orbit
DeltaT_fun = @(x) customTOF3(x(1),x(2),x(3), a_T, e_T, Th1_T, Th2_T, mu_s);

% Constraints for fmincon (only inequalities)
cond = @(x) deal([e_T(x(1),x(2),x(3)) - 1 + 1e-6;       % e < 1
    -e_T(x(1),x(2),x(3));                  % e >= 0
    (DeltaV_fun(x)/DeltaV_MAX - 1) ;% DeltaV <= DeltaV_MAX
    (DeltaT_fun(x)/DeltaT_MAX - 1) ;% DeltaT <= DeltaT_MAX
    (-(a_T(x(1),x(2),x(3)))*(1-e_T(x(1),x(2),x(3))) + 696340) / 1e6  % periapsis > R_SUN
    ], []);

% Objective function for optimization:
% weighted combination of normalized DeltaV and time of flight
ObjFun = @(w, x) w * (DeltaV_fun(x) / DeltaV_MAX) + (1-w) * (DeltaT_fun(x) / DeltaT_MAX);
% w = weight between 0 and 1 (0 = only time, 1 = only DeltaV)

% Bounds
lb = [0, 0, 0];           % lower bounds for [th_i, th_f, om_T]
ub = [2*pi, 2*pi, 2*pi];  % upper bounds for [th_i, th_f, om_T]

% Optimization options
opts_fmincon = optimoptions('fmincon', ...
    'Display', 'off', ...
    'Algorithm', 'sqp', ...
    'ConstraintTolerance', 1e-6, ...
    'OptimalityTolerance', 1e-6,...
    'FunctionTolerance',1e-8);


%% Min Delta V - fmincon. Run fmincon setup if not already done
clearvars -except AU mu_s a_i e_i i_i OM_i om_i th_i a_f e_f i_f OM_f om_f ...
    I_s m_0 DeltaV_MAX DeltaT_MAX ObjFun CosTh1_T ...
    CosTh2_T DeltaT_fun DeltaV_fun I J K N_T OM_T R_OM_T R_i_T R_om_T ...
    SinTh1_T SinTh2_T T_EliotoPF Th1_T Th2_T a_T e_T i_T ...
    cond h_T lb ub opts_fmincon rr1_T rr2_T rr_f ...
    rr_f_PF rr_i rr_i_PF vv1_T vv2_T vv_f vv_i


if ~exist('ObjFun', 'var')
    error('Run fmincon setup before continuing')
end
% Discretizing th_i, th_f and om_t
N_disc = 50; % Discretization value
th_1i = 0:pi/N_disc:2*pi;
th_2f = th_1i;
om_t = th_1i;
dim = size(th_1i,2);   % Vectors dimension

% Setting optimum equals to maximum
DeltaV_ott = DeltaV_MAX;

% Initializing vectors
RR_i = zeros(3, dim);
VV_i = zeros(3, dim);
RR_f = zeros(3, dim);
VV_f = zeros(3, dim);
R_om_comb = zeros(3, 1);% Rotation matrix of om_t


% Combining initial and final vectors rr & vv; creating rotation matrix of om_t
for idx = 1 : dim
    % Initial position and velocity
    [RR_i(:,idx),VV_i(:,idx)] = par2car(a_i,e_i,i_i,OM_i,om_i,th_1i(idx),mu_s); % r iniziale
    rv_i = [RR_i;VV_i];
    % Final position and velocity
    [RR_f(:,idx),VV_f(:,idx)] = par2car(a_f,e_f,i_f,OM_f,om_f,th_2f(idx),mu_s); % r finale
    rv_f = [RR_f;VV_f];

    % Construction of rotation matrix for every om_t
    R_om = [cos(om_t(idx)), sin(om_t(idx)), 0;
        -sin(om_t(idx)), cos(om_t(idx)), 0;
        0, 0, 1];

    R_om_comb(:, end+1:end+3) = R_om;
end

% Breaking down the rotation matrix from 2 to 3 dimensions
R_om_comb(:,1) = [];
R_om_blocks = reshape(R_om_comb, 3, 3, []);

% Calculating every combination of rr_i and rr_f:
comb = combinations(1:dim, 1:dim); % Every pair between 1 and dim
comb = table2array(comb);
comb_i = comb(:,1)';
comb_f = comb(:,2)';

rr_comb = [RR_i(:,comb_i);RR_f(:,comb_f)];  % [RR_i;RR_f]
vv_comb = [VV_i(:,comb_i);VV_f(:,comb_f)];  % [VV_i;VV_f]

% Reshaping vectors
rr_i_blocks = reshape(rr_comb(1:3,:),3,1,[]);
rr_f_blocks = reshape(rr_comb(4:6,:),3,1,[]);
vv_i_blocks = reshape(vv_comb(1:3,:),3,1,[]);
vv_f_blocks = reshape(vv_comb(4:6,:),3,1,[]);


% Transfer orbit parameters and Delta V calculation:
for j = 1 : size(comb,1)  % For every combination

    % Extracting one combination of initial and final ECI vectors
    rr_i_current = rr_i_blocks(:,:,j);
    rr_f_current = rr_f_blocks(:,:,j);
    vv_i_current = vv_i_blocks(:,:,j);
    vv_f_current = vv_f_blocks(:,:,j);

    rr_i_norm = norm(rr_i_current);
    rr_f_norm = norm(rr_f_current);

    % Transfer orbital plane direction
    h_t = cross(rr_i_current, rr_f_current)./norm(cross(rr_i_current, rr_f_current));

    % Transfer orbit inclination
    i_t = acos(h_t(3));

    % Transfer orbit node line
    N = cross(K, h_t)./norm(cross(K, h_t));

    % RAAN
    if N(2) >= 0
        OM_t = acos(N(1));
    else
        OM_t = 2*pi - acos(N(1));
    end

    % Rotation of OM around K
    R_OM = [cos(OM_t), sin(OM_t), 0;
        -sin(OM_t), cos(OM_t), 0;
        0, 0, 1];

    % Rotation of i around N
    R_i = [1,      0,       0;
        0, cos(i_t), sin(i_t);
        0, -sin(i_t), cos(i_t)];

    % Creating every possible rotation matrix
    R_comb = pagemtimes(R_om_blocks, R_i*R_OM);

    % Position vectors on PF frame, for every om_t
    rr_i_pf_current = pagemtimes(R_comb, rr_i_current);
    rr_f_pf_current = pagemtimes(R_comb, rr_f_current);

    % True anomaly of the initial transfer position
    cos_th_1t = rr_i_pf_current(1,:,:)./rr_i_norm;
    sin_th_1t = rr_i_pf_current(2,:,:)./rr_i_norm;
    th_1t = wrapTo2Pi(atan2(sin_th_1t(:),cos_th_1t(:)))';

    % True anomaly of the final transfer position
    cos_th_2t = rr_f_pf_current(1,:,:)./rr_f_norm;
    sin_th_2t = rr_f_pf_current(2,:,:)./rr_f_norm;
    th_2t = wrapTo2Pi(atan2(sin_th_2t(:),cos_th_2t(:)))';

    % Eccentricity modulus of current transfer orbit
    e_t_current =  (rr_i_norm - rr_f_norm)./(rr_f_norm .*  cos(th_2t) - rr_i_norm .* cos(th_1t));

    % It takes only closed orbits
    e_t_current(e_t_current>=1) = NaN;
    e_t_current(e_t_current<0) = NaN;

    % Semi-major axis
    a_t_current = (rr_i_norm.*(1+e_t_current.*cos(th_1t)))./(1-(e_t_current.^2));
    r_peri_current = a_t_current.*(1-e_t_current);

    % It takes only the orbits with a periapsis radius longer than the
    % radius of the Sun
    a_t_current(r_peri_current<= 696340) = NaN;

    % Calculation of vv_1t and vv_2t
    p = a_t_current.*(1-(e_t_current.^2));   % Current semi-latus rectum

    vv_1t_PF = sqrt(mu_s./p) .* [-sin(th_1t);
        e_t_current + cos(th_1t);
        zeros(1, length(th_1t))];
    vv_1t_PF = reshape(vv_1t_PF, 3, 1, []);
    vv_2t_PF= sqrt(mu_s./p) .* [-sin(th_2t);
        e_t_current + cos(th_2t);
        zeros(1, length(th_2t))];
    vv_2t_PF = reshape(vv_2t_PF, 3, 1, []);

    % Rotation matrix from PF to ECI
    R_PFtoECI = pagetranspose(R_comb);

    vv_1t = pagemtimes(R_PFtoECI, vv_1t_PF);
    vv_2t = pagemtimes(R_PFtoECI, vv_2t_PF);

    % Delta V calculation
    DeltaV1 = vecnorm(vv_1t - vv_i_current, 2, 1);
    DeltaV2 = vecnorm(vv_f_current - vv_2t, 2, 1);
    DeltaVV = DeltaV1(:) + DeltaV2(:);

    cosE1 = (e_t_current + cos(th_1t)) ./ (1 + e_t_current .* cos(th_1t));
    cosE2 = (e_t_current + cos(th_2t)) ./ (1 + e_t_current .* cos(th_2t));
    sinE1 = sqrt(1 - e_t_current.^2) .* sin(th_1t) ./ (1 + e_t_current .* cos(th_1t));
    sinE2 = sqrt(1 - e_t_current.^2) .* sin(th_2t) ./ (1 + e_t_current .* cos(th_2t));

    E1 = wrapTo2Pi(atan2(sinE1(:), cosE1(:)))';
    E2 = wrapTo2Pi(atan2(sinE2(:), cosE2(:)))';


    DeltaT_vec = sqrt(a_t_current.^3 ./ mu_s) .* (E2 - E1 - e_t_current .* (sinE2 - sinE1)) + 2.*pi.*sqrt(a_t_current.^3 ./ mu_s).*(th_2t<th_1t);

    % It takes the times of transfer that meets DeltaV requirements
    DeltaT_vec(DeltaVV>DeltaV_MAX) = NaN;
    DeltaT_vec(DeltaT_vec<0) = NaN;

    [DeltaVcandidate,index] = min(DeltaVV);

    if DeltaVcandidate < DeltaV_ott
        DeltaV_ott = DeltaVcandidate;
        om_t_ott = om_t(index);
        th_1i_ott = th_1i(ceil(j/length(th_1t)));
        th_2f_ott = th_2f(rem(j,length(th_1t)));
    end

end

x0 = [th_1i_ott; th_2f_ott; om_t_ott];%Initial guess for fmincon.m

[OPT, minV_rel] = fmincon(@(x) ObjFun(1, x), x0, [], [], [], [], lb, ub, cond, opts_fmincon);

minV = minV_rel * DeltaV_MAX;
% Compute the time taken for the transfer
DeltaT_minV = DeltaT_fun(OPT);

% Propellant mass calculation
m_p = m_0 * (exp(minV * 1e3/I_s) - 1);

% Displaing results
fprintf("================================================================\n"  + ...
    "Optimal  ΔV =  %.6f km/s (%.1f%% of  ΔV_MAX)\n" + ...
    "Respective  Δt = %.4f s = %.4f d \n" + ...
    "Propellant mass (with Is = 300 s & M_empty = 300 kg): %.2f kg\n" + ...
    "================================================================\n\n" ...
    , minV, 100*minV/DeltaV_MAX, DeltaT_minV, DeltaT_minV/(3600*24), m_p);


% Compute the cost of the two impulses - modulus and vector
Deltavv1_minV = vv1_T(OPT(1), OPT(2), OPT(3)) - vv_i(OPT(1));
DeltaV1_minV = norm(vv1_T(OPT(1), OPT(2), OPT(3)) - vv_i(OPT(1)));
Deltavv2_minV = vv_f(OPT(2)) - vv2_T(OPT(1), OPT(2), OPT(3));
DeltaV2_minV = norm(vv_f(OPT(2)) - vv2_T(OPT(1), OPT(2), OPT(3)));

% Saving transfer orbit parameters
a_opt = a_T(OPT(1),OPT(2),OPT(3));
e_opt = e_T(OPT(1),OPT(2),OPT(3));
i_opt = i_T(OPT(1),OPT(2));
OM_opt = OM_T(OPT(1),OPT(2));
om_opt = OPT(3);
th_1t_opt = Th1_T(OPT(1),OPT(2),OPT(3));
th_2t_opt = Th2_T(OPT(1),OPT(2),OPT(3));

% Table A6
A6(1,:) = [0, a_i/AU, e_i, i_i, OM_i, om_i, OPT(1), DeltaV1_minV];
A6(2,:) = [0, a_opt/AU, e_opt, i_opt, OM_opt, om_opt, th_1t_opt, 0];
A6(3,:) = [DeltaT_minV/3600, a_opt/AU, e_opt, i_opt, OM_opt, om_opt, th_2t_opt, DeltaV2_minV];
A6(4,:) = [DeltaT_minV/3600, a_f/AU, e_f, i_f, OM_f, om_f, OPT(2), 0];

TAB6 = array2table(A6, 'VariableNames', ...
    {'Time [h]', 'a [AU]', 'e [-]', 'i [rad]', 'OM [rad]', 'om [rad]', 'theta [rad]', 'DeltaV [km/s]'});
disp(TAB6)


% Plotting the resulting optimal transfer
figure('Name', 'A6 - Min Delta V fmincon','Units', 'normalized', 'OuterPosition', [0 0 1 1]);

% Plot the Sun
sun_sphere
grid on
hold on
box on
axis equal

% Initial orbit and position
plotOrbit(a_i, e_i, i_i, OM_i, om_i, 0, 2*pi, pi/10000, mu_s);
rr_start = rr_i(OPT(1));
scatter3(rr_start(1),rr_start(2),rr_start(3),'b*');

% Final orbit and position
plotOrbit(a_f, e_f, i_f, OM_f, om_f, 0, 2*pi, pi/10000, mu_s);
rr_finish = rr_f(OPT(2));
scatter3(rr_finish(1),rr_finish(2),rr_finish(3),'*');

% Transfer orbit
plotOrbit(a_opt, e_opt, i_opt, OM_opt, om_opt, 0, 2*pi, pi/10000, mu_s)
plotOrbit(a_opt, e_opt, i_opt, OM_opt, om_opt,th_1t_opt, th_2t_opt, pi/10000, mu_s, 3)

leg = legend('', 'Initial orbit', 'Initial position','Final orbit', 'Final position', 'Transfer orbit', 'Direct transfer', 'Interpreter', 'latex');
set(leg,'FontSize',22);
title('Direct transfer between the two orbits with $\Delta V_{MIN}$', 'Interpreter', 'latex')






%% Min Delta T - fmincon. Run fmincon setup if not already done
clearvars -except AU mu_s a_i e_i i_i OM_i om_i th_i a_f e_f i_f OM_f om_f ...
    I_s m_0 DeltaV_MAX DeltaT_MAX ObjFun CosTh1_T ...
    CosTh2_T DeltaT_fun DeltaV_fun I J K N_T OM_T R_OM_T R_i_T R_om_T ...
    SinTh1_T SinTh2_T T_EliotoPF Th1_T Th2_T a_T e_T i_T ...
    cond h_T lb ub opts_fmincon rr1_T rr2_T rr_f ...
    rr_f_PF rr_i rr_i_PF vv1_T vv2_T vv_f vv_i


if ~exist('ObjFun', 'var')
    error('Run fmincon setup before continuing')
end

% Discretizing th_i, th_f and om_t
N_disc = 50; % Discretization value
th_1i = 0:pi/N_disc:2*pi;
th_2f = th_1i;
om_t = th_1i;
dim = size(th_1i,2);   % Vectors dimension

% Setting optimum equals to maximum
DeltaT_ott = DeltaT_MAX;
% Initializing vectors
RR_i = zeros(3, dim);
VV_i = zeros(3, dim);
RR_f = zeros(3, dim);
VV_f = zeros(3, dim);
R_om_comb = zeros(3, 1);% Rotation matrix of om_t


% Combining initial and final vectors rr & vv; creating rotation matrix of om_t
for idx = 1 : dim
    % Initial position and velocity
    [RR_i(:,idx),VV_i(:,idx)] = par2car(a_i,e_i,i_i,OM_i,om_i,th_1i(idx),mu_s); % r iniziale
    rv_i = [RR_i;VV_i];
    % Final position and velocity
    [RR_f(:,idx),VV_f(:,idx)] = par2car(a_f,e_f,i_f,OM_f,om_f,th_2f(idx),mu_s); % r finale
    rv_f = [RR_f;VV_f];

    % Construction of rotation matrix for every om_t
    R_om = [cos(om_t(idx)), sin(om_t(idx)), 0;
        -sin(om_t(idx)), cos(om_t(idx)), 0;
        0, 0, 1];

    R_om_comb(:, end+1:end+3) = R_om;
end

% Breaking down the rotation matrix from 2 to 3 dimensions
R_om_comb(:,1) = [];
R_om_blocks = reshape(R_om_comb, 3, 3, []);

% Calculating every combination of rr_i and rr_f:
comb = combinations(1:dim, 1:dim); % Every pair between 1 and dim
comb = table2array(comb);
comb_i = comb(:,1)';
comb_f = comb(:,2)';

rr_comb = [RR_i(:,comb_i);RR_f(:,comb_f)];  % [RR_i;RR_f]
vv_comb = [VV_i(:,comb_i);VV_f(:,comb_f)];  % [VV_i;VV_f]

% Reshaping vectors
rr_i_blocks = reshape(rr_comb(1:3,:),3,1,[]);
rr_f_blocks = reshape(rr_comb(4:6,:),3,1,[]);
vv_i_blocks = reshape(vv_comb(1:3,:),3,1,[]);
vv_f_blocks = reshape(vv_comb(4:6,:),3,1,[]);


% Transfer orbit parameters and Delta V calculation:
for j = 1 : size(comb,1)  % For every combination

    % Extracting one combination of initial and final ECI vectors
    rr_i_current = rr_i_blocks(:,:,j);
    rr_f_current = rr_f_blocks(:,:,j);
    vv_i_current = vv_i_blocks(:,:,j);
    vv_f_current = vv_f_blocks(:,:,j);

    rr_i_norm = norm(rr_i_current);
    rr_f_norm = norm(rr_f_current);

    % Transfer orbital plane direction
    h_t = cross(rr_i_current, rr_f_current)./norm(cross(rr_i_current, rr_f_current));

    % Transfer orbit inclination
    i_t = acos(h_t(3));

    % Transfer orbit node line
    N = cross(K, h_t)./norm(cross(K, h_t));

    % RAAN
    if N(2) >= 0
        OM_t = acos(N(1));
    else
        OM_t = 2*pi - acos(N(1));
    end

    % Rotation of OM around K
    R_OM = [cos(OM_t), sin(OM_t), 0;
        -sin(OM_t), cos(OM_t), 0;
        0, 0, 1];

    % Rotation of i around N
    R_i = [1,      0,       0;
        0, cos(i_t), sin(i_t);
        0, -sin(i_t), cos(i_t)];

    % Creating every possible rotation matrix
    R_comb = pagemtimes(R_om_blocks, R_i*R_OM);

    % Position vectors on PF frame, for every om_t
    rr_i_pf_current = pagemtimes(R_comb, rr_i_current);
    rr_f_pf_current = pagemtimes(R_comb, rr_f_current);

    % True anomaly of the initial transfer position
    cos_th_1t = rr_i_pf_current(1,:,:)./rr_i_norm;
    sin_th_1t = rr_i_pf_current(2,:,:)./rr_i_norm;
    th_1t = atan2(sin_th_1t(:),cos_th_1t(:))';
    th_1t = wrapTo2Pi(th_1t);

    % True anomaly of the final transfer position
    cos_th_2t = rr_f_pf_current(1,:,:)./rr_f_norm;
    sin_th_2t = rr_f_pf_current(2,:,:)./rr_f_norm;
    th_2t = atan2(sin_th_2t(:),cos_th_2t(:))';
    th_2t = wrapTo2Pi(th_2t);

    % Eccentricity modulus of current transfer orbit
    e_t_current =  (rr_i_norm - rr_f_norm)./(rr_f_norm .*  cos(th_2t) - rr_i_norm .* cos(th_1t));

    % It takes only closed orbits
    e_t_current(e_t_current>=1) = NaN;
    e_t_current(e_t_current<0) = NaN;

    % Semi-major axis
    a_t_current = (rr_i_norm.*(1+e_t_current.*cos(th_1t)))./(1-(e_t_current.^2));
    r_peri_current = a_t_current.*(1-e_t_current);

    % It takes only the orbits with a periapsis radius longer than the
    % radius of the Sun
    a_t_current(r_peri_current<= 696340) = NaN;

    % Calculation of vv_1t and vv_2t
    p = a_t_current.*(1-(e_t_current.^2));   % Current semi-latus rectum

    vv_1t_PF = sqrt(mu_s./p) .* [-sin(th_1t);
        e_t_current + cos(th_1t);
        zeros(1, length(th_1t))];
    vv_1t_PF = reshape(vv_1t_PF, 3, 1, []);
    vv_2t_PF= sqrt(mu_s./p) .* [-sin(th_2t);
        e_t_current + cos(th_2t);
        zeros(1, length(th_2t))];
    vv_2t_PF = reshape(vv_2t_PF, 3, 1, []);

    % Rotation matrix from PF to ECI
    R_PFtoECI = pagetranspose(R_comb);

    vv_1t = pagemtimes(R_PFtoECI, vv_1t_PF);
    vv_2t = pagemtimes(R_PFtoECI, vv_2t_PF);

    % Delta V calculation
    DeltaV1 = vecnorm(vv_1t - vv_i_current, 2, 1);
    DeltaV2 = vecnorm(vv_f_current - vv_2t, 2, 1);
    DeltaVV = DeltaV1(:) + DeltaV2(:);

    cosE1 = (e_t_current + cos(th_1t)) ./ (1 + e_t_current .* cos(th_1t));
    cosE2 = (e_t_current + cos(th_2t)) ./ (1 + e_t_current .* cos(th_2t));
    sinE1 = sqrt(1 - e_t_current.^2) .* sin(th_1t) ./ (1 + e_t_current .* cos(th_1t));
    sinE2 = sqrt(1 - e_t_current.^2) .* sin(th_2t) ./ (1 + e_t_current .* cos(th_2t));

    E1 = wrapTo2Pi(atan2(sinE1(:), cosE1(:)))';
    E2 = wrapTo2Pi(atan2(sinE2(:), cosE2(:)))';


    DeltaT_vec = sqrt(a_t_current.^3 ./ mu_s) .* (E2 - E1 - e_t_current .* (sinE2 - sinE1)) + 2.*pi.*sqrt(a_t_current.^3 ./ mu_s).*(th_2t<th_1t);

    % It takes the times of transfer that meets DeltaV requirements
    DeltaT_vec(DeltaVV>DeltaV_MAX) = NaN;
    DeltaT_vec(DeltaT_vec<0) = NaN;

    [DeltaTcandidate, index] = min(DeltaT_vec);

    if DeltaTcandidate < DeltaT_ott
        DeltaT_ott = DeltaTcandidate;
        om_t_ott = om_t(index);
        th_1i_ott = th_1i(ceil(j/length(th_1t)));
        th_2f_ott = th_2f(rem(j,length(th_1t)));
    end

end



x0 = [th_1i_ott; th_2f_ott; om_t_ott];  %Initial guess for fmincon.m
[OPT, minT_perc] = fmincon(@(x) ObjFun(0, x) , x0, [], [], [], [], lb, ub, cond, opts_fmincon);
minT = minT_perc * DeltaT_MAX;

% Compute the cost, in terms of DeltaV, for the transfer
DeltaV_minT = DeltaV_fun(OPT);

% Propellant mass calculation
m_p = m_0 * (exp(DeltaV_minT * 1e3/I_s) - 1);

% Displaing results
fprintf("================================================================\n"  + ...
    "Optimal  ΔV =  %.6f km/s (%.1f%% of  ΔV_MAX)\n" + ...
    "Respective  Δt = %.4f s = %.4f d \n" + ...
    "Propellant mass (with Is = 300 s & M_empty = 300 kg): %.2f kg\n" + ...
    "================================================================\n\n" ...
    , DeltaV_minT, 100*DeltaV_minT/DeltaV_MAX, minT, minT/(3600*24), m_p);


% Compute the cost of the two impulses - modulus and vector
Deltavv1_minT = vv1_T(OPT(1), OPT(2), OPT(3)) - vv_i(OPT(1));
DeltaV1_minT = norm(vv1_T(OPT(1), OPT(2), OPT(3)) - vv_i(OPT(1)));
Deltavv2_minT = vv_f(OPT(2)) - vv2_T(OPT(1), OPT(2), OPT(3));
DeltaV2_minT = norm(vv_f(OPT(2)) - vv2_T(OPT(1), OPT(2), OPT(3)));

% Saving transfer orbit parameters
a_opt = a_T(OPT(1),OPT(2),OPT(3));
e_opt = e_T(OPT(1),OPT(2),OPT(3));
i_opt = i_T(OPT(1),OPT(2));
OM_opt = OM_T(OPT(1),OPT(2));
om_opt = OPT(3);
th_1t_opt = Th1_T(OPT(1),OPT(2),OPT(3));
th_2t_opt = Th2_T(OPT(1),OPT(2),OPT(3));



% Table A7
A7(1,:) = [0, a_i/AU, e_i, i_i, OM_i, om_i, OPT(1), DeltaV1_minT];
A7(2,:) = [0, a_opt/AU, e_opt, i_opt, OM_opt, om_opt, th_1t_opt, 0];
A7(3,:) = [minT/3600, a_opt/AU, e_opt, i_opt, OM_opt, om_opt, th_2t_opt, DeltaV2_minT];
A7(4,:) = [minT/3600, a_f/AU, e_f, i_f, OM_f, om_f, OPT(2), 0];

TAB7 = array2table(A7, 'VariableNames', ...
    {'Time [h]', 'a [AU]', 'e [-]', 'i [rad]', 'OM [rad]', 'om [rad]', 'theta [rad]', 'DeltaV [km/s]'});
disp(TAB7)


% Plotting the resulting optimal transfer
figure('Name', 'A7 - Min Delta t fmincon','Units', 'normalized', 'OuterPosition', [0 0 1 1]);

% Plot the Sun
sun_sphere
grid on
hold on
box on
axis equal

% Initial orbit and position
plotOrbit(a_i, e_i, i_i, OM_i, om_i, 0, 2*pi, pi/10000, mu_s);
rr_start = rr_i(OPT(1));
scatter3(rr_start(1), rr_start(2), rr_start(3),'b*');

% Final orbit and position
plotOrbit(a_f, e_f, i_f, OM_f, om_f, 0, 2*pi, pi/10000, mu_s);
rr_finish = rr_f(OPT(2));
scatter3(rr_finish(1), rr_finish(2), rr_finish(3),'*');

% Transfer orbit
plotOrbit(a_opt, e_opt, i_opt, OM_opt, om_opt, 0, 2*pi, pi/10000, mu_s)
plotOrbit(a_opt, e_opt, i_opt, OM_opt, om_opt,th_1t_opt, th_2t_opt, pi/10000, mu_s, 3)

leg = legend('', 'Initial orbit', 'Initial position','Final orbit', 'Final position', 'Transfer orbit', 'Direct transfer', 'Interpreter', 'latex');
set(leg,'FontSize',22);
title('Direct transfer between the two orbits with $\Delta t_{MIN}$', 'Interpreter', 'latex')


%% Optimal transfer - fmincon. Run fmincon setup if not already done

clearvars -except AU mu_s a_i e_i i_i OM_i om_i th_i a_f e_f i_f OM_f om_f ...
    I_s m_0 DeltaV_MAX DeltaT_MAX ObjFun CosTh1_T ...
    CosTh2_T DeltaT_fun DeltaV_fun I J K N_T OM_T R_OM_T R_i_T R_om_T ...
    SinTh1_T SinTh2_T T_EliotoPF Th1_T Th2_T a_T e_T i_T ...
    cond h_T lb ub opts_fmincon rr1_T rr2_T rr_f ...
    rr_f_PF rr_i rr_i_PF vv1_T vv2_T vv_f vv_i


if ~exist('ObjFun', 'var')
    error('Run fmincon setup before continuing')
end

% Discretizing th_i, th_f and om_t
N_disc = 50; % Discretization value
th_1i = 0:pi/N_disc:2*pi;
th_2f = th_1i;
om_t = th_1i;
dim = size(th_1i,2);   % Vectors dimension

% Setting optimum equals to maximum
DeltaT_ott = DeltaT_MAX;
DeltaV_ott = DeltaV_MAX;

% Initializing vectors
RR_i = zeros(3, dim);
VV_i = zeros(3, dim);
RR_f = zeros(3, dim);
VV_f = zeros(3, dim);
R_om_comb = zeros(3, 1);% Rotation matrix of om_t

w = linspace(0, 1, 33); %Vector of weights in ObjFun(w,x)
% w = linspace(0.8258,0.8264, 14); De-comment for focused analysis
% w = 0.82; % De-comment for mixed optimum


Ts = zeros(1,length(w));
Vs = zeros(1,length(w));


for m = 1:length(w)
    
    fprintf("Iterazione %d di %d \nPeso = %.4f\n", m, length(w), w(m))
    % Setting optimum equals to maximum
    DeltaV_ott = DeltaV_MAX;
    DeltaT_ott = DeltaT_MAX;
    % Initializing vectors

    RR_i = zeros(3, dim);
    VV_i = zeros(3, dim);
    RR_f = zeros(3, dim);
    VV_f = zeros(3, dim);
    R_om_comb = zeros(3, 1);% Rotation matrix of om_t

    OTT = w(m) * DeltaV_ott / DeltaV_MAX + (1 - w(m)) * DeltaT_ott / DeltaT_MAX;

    for idx = 1 : dim
        % Initial position and velocity
        [RR_i(:,idx),VV_i(:,idx)] = par2car(a_i,e_i,i_i,OM_i,om_i,th_1i(idx),mu_s); % r iniziale
        rv_i = [RR_i;VV_i];
        % Final position and velocity
        [RR_f(:,idx),VV_f(:,idx)] = par2car(a_f,e_f,i_f,OM_f,om_f,th_2f(idx),mu_s); % r finale
        rv_f = [RR_f;VV_f];

        % Construction of rotation matrix for every om_t
        R_om = [cos(om_t(idx)), sin(om_t(idx)), 0;
            -sin(om_t(idx)), cos(om_t(idx)), 0;
            0, 0, 1];

        R_om_comb(:, end+1:end+3) = R_om;
    end

    % Breaking down the rotation matrix from 2 to 3 dimensions
    R_om_comb(:,1) = [];
    R_om_blocks = reshape(R_om_comb, 3, 3, []);

    % Calculating every combination of rr_i and rr_f:
    comb = combinations(1:dim, 1:dim); % Every pair between 1 and dim
    comb = table2array(comb);
    comb_i = comb(:,1)';
    comb_f = comb(:,2)';

    rr_comb = [RR_i(:,comb_i);RR_f(:,comb_f)];  % [RR_i;RR_f]
    vv_comb = [VV_i(:,comb_i);VV_f(:,comb_f)];  % [VV_i;VV_f]

    % Reshaping vectors
    rr_i_blocks = reshape(rr_comb(1:3,:),3,1,[]);
    rr_f_blocks = reshape(rr_comb(4:6,:),3,1,[]);
    vv_i_blocks = reshape(vv_comb(1:3,:),3,1,[]);
    vv_f_blocks = reshape(vv_comb(4:6,:),3,1,[]);


    % Transfer orbit parameters and Delta V calculation:
    for j = 1 : size(comb,1)  % For every combination

        % Extracting one combination of initial and final ECI vectors
        rr_i_current = rr_i_blocks(:,:,j);
        rr_f_current = rr_f_blocks(:,:,j);
        vv_i_current = vv_i_blocks(:,:,j);
        vv_f_current = vv_f_blocks(:,:,j);

        rr_i_norm = norm(rr_i_current);
        rr_f_norm = norm(rr_f_current);

        % Transfer orbital plane direction
        h_t = cross(rr_i_current, rr_f_current)./norm(cross(rr_i_current, rr_f_current));

        % Transfer orbit inclination
        i_t = acos(h_t(3));

        % Transfer orbit node line
        N = cross(K, h_t)./norm(cross(K, h_t));

        % RAAN
        if N(2) >= 0
            OM_t = acos(N(1));
        else
            OM_t = 2*pi - acos(N(1));
        end

        % Rotation of OM around K
        R_OM = [cos(OM_t), sin(OM_t), 0;
            -sin(OM_t), cos(OM_t), 0;
            0, 0, 1];

        % Rotation of i around N
        R_i = [1,      0,       0;
            0, cos(i_t), sin(i_t);
            0, -sin(i_t), cos(i_t)];

        % Creating every possible rotation matrix
        R_comb = pagemtimes(R_om_blocks, R_i*R_OM);

        % Position vectors on PF frame, for every om_t
        rr_i_pf_current = pagemtimes(R_comb, rr_i_current);
        rr_f_pf_current = pagemtimes(R_comb, rr_f_current);

        % True anomaly of the initial transfer position
        cos_th_1t = rr_i_pf_current(1,:,:)./rr_i_norm;
        sin_th_1t = rr_i_pf_current(2,:,:)./rr_i_norm;
        th_1t = atan2(sin_th_1t(:),cos_th_1t(:))';
        th_1t = wrapTo2Pi(th_1t);

        % True anomaly of the final transfer position
        cos_th_2t = rr_f_pf_current(1,:,:)./rr_f_norm;
        sin_th_2t = rr_f_pf_current(2,:,:)./rr_f_norm;
        th_2t = atan2(sin_th_2t(:),cos_th_2t(:))';
        th_2t = wrapTo2Pi(th_2t);

        % Eccentricity modulus of current transfer orbit
        e_t_current =  (rr_i_norm - rr_f_norm)./(rr_f_norm .*  cos(th_2t) - rr_i_norm .* cos(th_1t));

        % It takes only closed orbits
        e_t_current(e_t_current>=1) = NaN;
        e_t_current(e_t_current<0) = NaN;

        % Semi-major axis
        a_t_current = (rr_i_norm.*(1+e_t_current.*cos(th_1t)))./(1-(e_t_current.^2));
        r_peri_current = a_t_current.*(1-e_t_current);

        % It takes only the orbits with a periapsis radius longer than the
        % radius of the Sun
        a_t_current(r_peri_current<= 696340) = NaN;

        % Calculation of vv_1t and vv_2t
        p = a_t_current.*(1-(e_t_current.^2));   % Current semi-latus rectum

        vv_1t_PF = sqrt(mu_s./p) .* [-sin(th_1t);
            e_t_current + cos(th_1t);
            zeros(1, length(th_1t))];
        vv_1t_PF = reshape(vv_1t_PF, 3, 1, []);
        vv_2t_PF= sqrt(mu_s./p) .* [-sin(th_2t);
            e_t_current + cos(th_2t);
            zeros(1, length(th_2t))];
        vv_2t_PF = reshape(vv_2t_PF, 3, 1, []);

        % Rotation matrix from PF to ECI
        R_PFtoECI = pagetranspose(R_comb);

        vv_1t = pagemtimes(R_PFtoECI, vv_1t_PF);
        vv_2t = pagemtimes(R_PFtoECI, vv_2t_PF);

        % Delta V calculation
        DeltaV1 = vecnorm(vv_1t - vv_i_current, 2, 1);
        DeltaV2 = vecnorm(vv_f_current - vv_2t, 2, 1);
        DeltaVV = DeltaV1(:) + DeltaV2(:);

        cosE1 = (e_t_current + cos(th_1t)) ./ (1 + e_t_current .* cos(th_1t));
        cosE2 = (e_t_current + cos(th_2t)) ./ (1 + e_t_current .* cos(th_2t));
        sinE1 = sqrt(1 - e_t_current.^2) .* sin(th_1t) ./ (1 + e_t_current .* cos(th_1t));
        sinE2 = sqrt(1 - e_t_current.^2) .* sin(th_2t) ./ (1 + e_t_current .* cos(th_2t));

        E1 = wrapTo2Pi(atan2(sinE1(:), cosE1(:)))';
        E2 = wrapTo2Pi(atan2(sinE2(:), cosE2(:)))';


        DeltaT_vec = sqrt(a_t_current.^3 ./ mu_s) .* (E2 - E1 - e_t_current .* (sinE2 - sinE1)) + 2.*pi.*sqrt(a_t_current.^3 ./ mu_s).*(th_2t<th_1t);

        % It takes the times of transfer that meets DeltaV requirements
        DeltaT_vec(DeltaVV>DeltaV_MAX) = NaN;
        DeltaT_vec(DeltaT_vec<0) = NaN;

        Obj_vec = w(m) * DeltaVV / DeltaV_MAX + (1-w(m)) * DeltaT_vec' / DeltaT_MAX;

        Obj_vec(Obj_vec<0) = NaN;
        Obj_vec(DeltaVV>DeltaV_MAX) = NaN;
        Obj_vec(DeltaT_vec<0) = NaN;
        Obj_vec(DeltaT_vec>DeltaT_MAX) = NaN;

        [Opt_candidate,index] = min(Obj_vec);

        if Opt_candidate < OTT
            OTT = Opt_candidate;
            om_t_ott = om_t(index);
            th_1i_ott = th_1i(ceil(j/length(th_1t)));
            th_2f_ott = th_2f(rem(j,length(th_1t)));

        end

    end

    th_1i_ott = wrapTo2Pi(th_1i_ott);
    th_2f_ott = wrapTo2Pi(th_2f_ott);

    x0(m,:) = [th_1i_ott, th_2f_ott, om_t_ott]; % Initial guess for fmincon.m

    [OPT, ~] = fmincon(@(x) ObjFun(w(m),x), x0(m,:), [], [], [], [], lb, ub, cond, opts_fmincon);

    % Compute the time taken for the transfer
    DeltaT_MIX = DeltaT_fun(OPT);
    % Compute the cost, in terms of DeltaV, for the transfer
    DeltaV_MIX = DeltaV_fun(OPT);

    Ts(m) = DeltaT_MIX/3600;
    Vs(m) = DeltaV_MIX;

end

figure('Name', 'Optimization','Units', 'normalized', 'OuterPosition', [0 0 1 1]);

grid on
yyaxis left
scatter(w, Vs, 'filled');
set(ylabel('$\Delta V$ [km/s]', 'Interpreter', 'latex'), 'Rotation', 0);

yyaxis right
scatter(w, Ts, 'filled');
set(ylabel('$\Delta t$ [h]', 'Interpreter', 'latex'), 'Rotation', 0);


title('Variation of $\Delta V$ and $\Delta t$ with weight $w$ ', 'Interpreter', 'latex')
leg = legend('$\Delta V$', '$\Delta t$','interpreter','latex','location','east');
set(leg,'FontSize',22);
box on
axis padded
ax = gca;        
ax.FontSize = 18; 
xlabel('$w$','interpreter','latex','FontSize',22);







