%%%%%%%%%%%%%%%%%%%%%%%% General Setup - Run first %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                                                   %%%%%%%%%%%%
%%%%%%%%%%%%                    GROUP 22                       %%%%%%%%%%%%
%%%%%%%%%%%%                                                   %%%%%%%%%%%%
%%%%%%%%%%%%            Giovanni Pieroni - 10854347            %%%%%%%%%%%%
%%%%%%%%%%%%            Alessandro Ponti - 10901276            %%%%%%%%%%%%
%%%%%%%%%%%%            Jacopo Rotta     - 10696753            %%%%%%%%%%%%
%%%%%%%%%%%%    SC03 - Escape and arrival trajectories with    %%%%%%%%%%%%
%%%%%%%%%%%%               patched conics method               %%%%%%%%%%%%
%%%%%%%%%%%%                                                   %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
clc

%% SETUP - This section (SC2 - DeltaV min) needs to be executed before calculating SC3 (following sections) 

format short g

% Gravitational constant [km^3/(kg*s^2)]
G = 6.67259e-20; 

% Masses [kg]
M_T = 5.974e24;
M_S = 1.989e30;
M_A = 8.3961e11;

% Gravitational parameters [km^3/s^2]
mu_s = 132712440018;   
mu_t = 398600; 
mu_A = G * (M_A);

% Astronoical unit [km]
AU = 149597870.7;

% 1991 VH orbit
a_A = 1.137342 * AU;    % Semi-major axis [km]
e_A = 0.144236;         % Eccentricity [-]
i_A = deg2rad(13.91);   % Inclination [rad]
OM_A = deg2rad(139.33); % RAAN [rad]
om_A = deg2rad(206.94); % pericenter argument [rad]
D = 0.93;               % Diameter [km]


% SOI [km]
rSOI_T = 149597870.7 * (M_T/M_S)^(2/5);
rSOI_A = a_A * (M_A/M_S)^(2/5);

Dati = importdata('DatiSC1-2025.mat','T');
SC_1 = table2array(Dati(22,:));


% Inital orbit data (end of SC1)
rr_SC1 = (SC_1(8:10))';
vv_SC1 = (SC_1(11:13))';
[a_SC1, e_SC1, i_SC1, OM_SC1, om_SC1, ~] = car2par(rr_SC1, vv_SC1, mu_t);
th_SC1 = 0;
vp_SC1 = sqrt(mu_t/(a_SC1*(1 - e_SC1^2))) * (1 + e_SC1);


% Earth orbit
a_i = 1.00000011 * AU; % Semi-major axix [km]
e_i = 0.01671022;      % Eccentricity [-]
i_i = 9.1920e-5;       % Inclination [rad]
OM_i = 2.7847;         % RAAN [rad]
om_i = 5.2643;         % pericenter argument [rad]


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
% % %  % Everything (if needed) is function of th_i, th_f and om_t % % % %

% Initial position and velocity
[rr_i, vv_i] = par2carFUN1th(a_i, e_i, i_i, OM_i, om_i, mu_s);

% Final position and velocity
[rr_f, vv_f] = par2carFUN1th(a_A, e_A, i_A, OM_A, om_A, mu_s);


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
cond_multi = @(x) deal([e_T(x(1),x(2),x(3)) - 1 + 1e-6;       % e < 1
    -e_T(x(1),x(2),x(3));                  % e >= 0
    (DeltaV_fun(x)/DeltaV_MAX - 1) ;% DeltaV <= DeltaV_MAX
    (DeltaT_fun(x)/DeltaT_MAX - 1) ;% DeltaT <= DeltaT_MAX
    (-(a_T(x(1),x(2),x(3)))*(1-e_T(x(1),x(2),x(3))) + 696340) / 1e6  % periapsis > R_SUN
    ], []);

% Objective function for optimization:
% weighted combination of normalized DeltaV and time of flight
ObjectiveFunction = @(w, x) w * (DeltaV_fun(x) / DeltaV_MAX) + (1-w) * (DeltaT_fun(x) / DeltaT_MAX);
% w = weight between 0 and 1 (0 = only time, 1 = only DeltaV)

% Bounds
lb = [0, 0, 0];           % lower bounds for [th_i, th_f, om_T]
ub = [2*pi, 2*pi, 2*pi];  % upper bounds for [th_i, th_f, om_T]

% Optimization options
opts_fmincon_multi = optimoptions('fmincon', ...
    'Display', 'off', ...
    'Algorithm', 'sqp', ...
    'ConstraintTolerance', 1e-16, ...
    'FunctionTolerance',1e-16);

if ~exist('ObjectiveFunction', 'var')
    error('Run fmincon setup before continuing')
end
% Discretizing th_i, th_f and om_t
N_disc = 100; % Discretization value
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
    [RR_f(:,idx),VV_f(:,idx)] = par2car(a_A,e_A,i_A,OM_A,om_A,th_2f(idx),mu_s); % r finale
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


    DeltaT_vec = sqrt(a_t_current.^3 ./ mu_s) .* (E2 - E1 - e_t_current .* (sinE2 - sinE1)) + 2.*pi.*sqrt(a_t_current.^3 ./ mu_s).*(th_2t>th_1t);

    % It takes the times of transfer that meets DeltaV requirements
    DeltaT_vec(DeltaVV>DeltaV_MAX) = NaN;


    [DeltaVcandidate,index] = min(DeltaVV);

    if DeltaVcandidate < DeltaV_ott
        DeltaV_ott = DeltaVcandidate;
        om_t_ott = om_t(index);
        th_1i_ott = th_1i(ceil(j/length(th_1t)));
        th_2f_ott = th_2f(rem(j,length(th_1t)));
    end

end

x0 = [th_1i_ott; th_2f_ott; om_t_ott];% Initial guess for fmincon.m

[OPT, minV_rel] = fmincon(@(x) ObjectiveFunction(1, x), x0, [], [], [], [], lb, ub, cond_multi, opts_fmincon_multi);


% Compute the time taken for the transfer
DeltaT_minV = DeltaT_fun(OPT);


% Compute the cost of the two impulses - modulus and vector
Deltavv1_minV = vv1_T(OPT(1), OPT(2), OPT(3)) - vv_i(OPT(1));
DeltaV1_minV = norm(vv1_T(OPT(1), OPT(2), OPT(3)) - vv_i(OPT(1)));
Deltavv2_minV = vv_f(OPT(2)) - vv2_T(OPT(1), OPT(2), OPT(3));
DeltaV2_minV = norm(vv_f(OPT(2)) - vv2_T(OPT(1), OPT(2), OPT(3)));


%% SC3 with coplanarity hypotesis





% Transfer orbit form SC2
a_t = a_T(OPT(1),OPT(2),OPT(3));
e_t = e_T(OPT(1),OPT(2),OPT(3));
i_t = i_T(OPT(1),OPT(2));
OM_t = OM_T(OPT(1),OPT(2));
om_t = OPT(3);
Th1_t = Th1_T(OPT(1),OPT(2),OPT(3));
Th2_t = Th2_T(OPT(1),OPT(2),OPT(3));


% Excess hyperbolic speed
vinf_1 = norm(Deltavv1_minV);
vinf_2 = norm(Deltavv2_minV);

% Escape hyperbola
a_H1 = -mu_t/vinf_1^2;

rp_H1 = a_SC1*(1 - e_SC1);

e_H1 = 1 - rp_H1/a_H1;

v_fuga1 = sqrt(2*mu_t/rp_H1);

vp_H1 = sqrt(v_fuga1^2 + vinf_1^2);

DeltaV_H1 = abs(vp_H1 - vp_SC1);


% Arrival hyperbola to a circular orbit
a_H2 = - mu_A/vinf_2^2;

% Testing varius r_circ from the surface to r_SOI
rp_H2 = D/2 : 1 : rSOI_A;

e_H2 = 1 - rp_H2./a_H2;

v_fuga2 = sqrt(2*mu_A./rp_H2);

vp_H2 = sqrt(v_fuga2.^2 + vinf_2^2);

v_circ = sqrt(mu_A./rp_H2);

DeltaV_H2 = abs(v_circ - vp_H2);

Delta_impact = - a_H2 .* sqrt(e_H2.^2 - 1);

% Total Delta V for the patched conics method (circular arrival orbit)
DeltaV_tot = DeltaV_H1 + DeltaV_H2;



% Arrival hyperbola to an elliptical orbit (optimizing DeltaV)

e_H2_ell = @(rp_H2) 1 - rp_H2./a_H2;

v_fuga2_ell = @(rp_H2) sqrt(2*mu_A./rp_H2);

vp_H2_ell = @(rp_H2) sqrt(v_fuga2_ell(rp_H2).^2 + vinf_2.^2);

rp_ellisse = @(rp_H2) rp_H2;

margin = 0.1; % km

ra_ellisse = rSOI_A ;

e_ellisse =  @(rp_H2) (ra_ellisse - rp_ellisse(rp_H2))./(ra_ellisse + rp_ellisse(rp_H2));

a_ellisse = @(rp_H2) rp_ellisse(rp_H2)  ./ (1 - e_ellisse(rp_H2));

p_ellisse = @(rp_H2) a_ellisse(rp_H2)  .* (1 - e_ellisse(rp_H2).^2);

vp_ellisse = @(rp_H2) sqrt(mu_A./p_ellisse(rp_H2)) .* (1 + e_ellisse(rp_H2));

DeltaV_H2_ell = @(rp_H2) abs(vp_ellisse(rp_H2)  - vp_H2_ell(rp_H2));

opts_ga = optimoptions('ga', 'Display', 'off', 'PlotFcn', []);

% Upper and lower bound for rp_H2
lb = D/2 + margin;           
ub = rSOI_A;  

rng(42); % rng seed for consistent results

% Initial guess for fmincon
[x0, ~] = ga(DeltaV_H2_ell, 1,[],[],[],[], lb, ub, [], opts_ga);
opts_fmincon = optimoptions('fmincon', 'Display', 'off', 'PlotFcn', []);

% Minimixing DeltaV
[OPT_ell, ~] = fmincon(DeltaV_H2_ell, x0, [], [], [], [], lb, ub, [], opts_fmincon);

rp_H2_ell = OPT_ell;

% Total Delta V for the patched conics method (elliptic arrival orbit)
DeltaV_ell = DeltaV_H1 + DeltaV_H2_ell(OPT_ell);

Delta_impact_ell = -a_H2 * sqrt(e_H2_ell(rp_H2_ell)^2 -1);




% Plotting results

% Delta V plot 
figure('Name','Delta V comparison','Units', 'normalized', 'OuterPosition', [0 0 1 1])
subplot(2,1,1)
hold on
grid on
box on
scatter(rp_H2-D/2, DeltaV_tot,'b','filled')
yline(DeltaV_ell,'k--','LineWidth',1.5)
xline(rSOI_A,'--r','LineWidth',1.5)
xlabel("Km from the surface [Km]",Interpreter="latex")
ylabel("$\Delta V$ [Km/s]",Interpreter="latex")
leg = legend('$\Delta V_{tot}$  circular orbit','$\Delta V_{tot}$  best elliptical orbit','SOI','Interpreter','latex');
set(leg,'FontSize',22,'Location','best');

axis padded

% Delta impact plot
subplot(2,1,2)
hold on
grid on
box on
scatter(rp_H2-D/2, Delta_impact,'b','filled')
yline(Delta_impact_ell,'k--','LineWidth',1.5)
xline(rSOI_A,'--r','LineWidth',1.5)
xlabel("Km from the surface [Km]",Interpreter="latex")
ylabel('$\Delta$ [Km]','Interpreter','latex')
leg = legend('$\Delta$  circular orbit','$\Delta$  best elliptical orbit','SOI','Interpreter','latex');
set(leg,'FontSize',22,'Location','best');

axis padded







% Escape trajectory
figure('Name','Escape trajectory','Units', 'normalized', 'OuterPosition', [0 0 1 1])
hold on
earth_sphere('km')
grid on
box on
axis equal

% Initial orbit
plotOrbit(a_SC1, e_SC1, i_SC1, OM_SC1, om_SC1, 0, 2*pi, pi/1000, mu_t)

% Escape hyperbola
plotOrbitC(a_H1, e_H1, i_SC1, OM_SC1, om_SC1, 0, pi/3, pi/1000, mu_t, 43,[0.1 0.7 0.1],'-',0);

% Perigee
[rr_pSC1, ~] = par2car(a_SC1, e_SC1, i_SC1, OM_SC1, om_SC1, th_SC1, mu_t);
scatter3(rr_pSC1(1),rr_pSC1(2),rr_pSC1(3),'sk')

% Earth's SOI
plotOrbitC(rSOI_T, 0, i_SC1, OM_SC1, om_SC1, 0, 2*pi, pi/1000, mu_t,1,'r','-')

% Intersection with SOI
th_H1 = acos(a_H1 / e_H1 / rSOI_T * (1 - e_H1^2) - 1 / e_H1);
[rr_H1, ~] = par2car(a_H1, e_H1, i_SC1, OM_SC1, om_SC1, th_H1, mu_t);
scatter3(rr_H1(1),rr_H1(2),rr_H1(3),'ok')


leg = legend('',"Initial orbit", "Hyperbolic orbit","Perigee", "Earth SOI", "Intersection with Earth SOI",'Interpreter','latex');
set(leg,'FontSize', 22)



% Arrival trajectory
figure('Name','Arrival trajectory','Units', 'normalized', 'OuterPosition', [0 0 1 1])
hold on
grid on
box on
axis equal
view(2)

% Asteroid
[X_A,Y_A,Z_A] = sphere(50);
surf((D/2)*X_A, (D/2)*Y_A,(D/2)*Z_A,'FaceColor', [0.7, 0.7, 0.7], 'EdgeColor', 'none');

% SOI and intersections
plotOrbitC(rSOI_A, 0, i_t, OM_t, om_t, 0, 2*pi, pi/1000, mu_A, 2,'r','--')
th_H2 = acos(a_H2 ./ e_H2 / rSOI_A .* (1 - e_H2.^2) - 1 ./ e_H2);

% circular orbits
for j =1:length(e_H2)
plotOrbitC(a_H2, e_H2(j), i_t, OM_t, om_t, 0, 2*pi, pi/1000, mu_A, 40/(1.4^(j-1)),'g','-',0)
plotOrbitC(rp_H2(j), 0, i_t , OM_t, om_t, 0, 2*pi, pi/1000, mu_A,1 ,[0.1 0.7 0.1])
[rr_H2, ~] = par2car(a_H2, e_H2(j), i_t, OM_t, om_t, th_H2(j), mu_A);
scatter3(rr_H2(1), rr_H2(2), rr_H2(3),'Marker','o','MarkerEdgeColor','k','LineWidth',1.5)
end

% elliptic orbit
plotOrbitC(a_H2, e_H2_ell(rp_H2_ell), i_t, OM_t, om_t, 0, 2*pi, pi/1000, mu_A, 40,'b','-',0)
plotOrbitC(a_ellisse(rp_H2_ell), e_ellisse(rp_H2_ell), i_t, OM_t, om_t, 0, 2*pi, pi/1000, mu_A, 1,'c')


leg = legend('1991 VH','SOI','Arrival trajectory (circular)','Circular orbits','SOI intersections', ...
    '','','','','','','','','','','','','','','','','','','','','', ...
  'Arrival trajectory (elliptic)','Elliptic orbit','Interpreter','latex');
set(leg, 'FontSize',22,'Location','bestoutside');

xlim([-15 15])
ylim([-15 15])




figure('Name','Elliptical DeltaV');
hold on
grid on
box on
rr = linspace(D/2, rSOI_A, 1000);
plot(rr, DeltaV_H2_ell(rr), 'LineWidth',1.5)
xline(D/2,'r--')
xline(rSOI_A,'k--')

legend('$\Delta V_H$', '$R_A$','SOI', 'interpreter', 'Latex', 'FontSize',30,'Location','northwest')
xlabel('Pericenter distance $r_{p,f}$', 'Interpreter','latex', 'FontSize', 24)
ylabel('Required $\Delta V$','Interpreter','latex', 'FontSize', 24)


% Escape trajectory table
A8(1,:) = [0,a_SC1,e_SC1,i_SC1,OM_SC1,om_SC1,0,0];
A8(2,:) = [0,a_SC1,e_SC1,i_SC1,OM_SC1,om_SC1,0,DeltaV_H1];
A8(3,:) = [0,a_H1,e_H1,i_SC1,OM_SC1,om_SC1,0,0];

TAB8 = array2table(A8, 'VariableNames', ...
    {'Time [h]', 'a [AU]', 'e [-]', 'i [rad]', 'OM [rad]', 'om [rad]', 'theta [rad]', 'DeltaV [km/s]'});
disp(TAB8)

% Elliptical arrival table
A10(1,:) = [0,a_H2,e_H2_ell(OPT_ell),i_T(OPT(1),OPT(2)),OM_T(OPT(1),OPT(2)),OPT(3),0,0];
A10(2,:) = [0,a_H2,e_H2_ell(OPT_ell),i_T(OPT(1),OPT(2)),OM_T(OPT(1),OPT(2)),OPT(3),0,DeltaV_H2_ell(OPT_ell)];
A10(3,:) = [0,a_ellisse(OPT_ell),e_ellisse(OPT_ell),i_T(OPT(1),OPT(2)),OM_T(OPT(1),OPT(2)),OPT(3),0,0];

TAB10 = array2table(A10, 'VariableNames', ...
    {'Time [h]', 'a [AU]', 'e [-]', 'i [rad]', 'OM [rad]', 'om [rad]', 'theta [rad]', 'DeltaV [km/s]'});
disp(TAB10)

% Circular arrival table (1 km above the surface)
A11(1,:) = [0,a_H2,e_H2(2),i_T(OPT(1),OPT(2)),OM_T(OPT(1),OPT(2)),OPT(3),0,0];
A11(2,:) = [0,a_H2,e_H2(2),i_T(OPT(1),OPT(2)),OM_T(OPT(1),OPT(2)),OPT(3),0,DeltaV_H2(2)];
A11(3,:) = [0,rp_H2(2),0,i_T(OPT(1),OPT(2)),OM_T(OPT(1),OPT(2)),OPT(3),0,0];

TAB11 = array2table(A11, 'VariableNames', ...
    {'Time [h]', 'a [AU]', 'e [-]', 'i [rad]', 'OM [rad]', 'om [rad]', 'theta [rad]', 'DeltaV [km/s]'});
disp(TAB11)

%% SC3 with plane change to match v_inf1 direction

clear non_conv_iter flags


opts = optimoptions('fsolve','FunctionTolerance',1e-16,'Display','off');


Dati = importdata('DatiSC1-2025.mat','T');
SC_1 = table2array(Dati(22,:));

% Earth's axis inclination
epsilon = deg2rad(23.45);
T = [1 0 0; 
     0 cos(epsilon) sin(epsilon); 
     0 -sin(epsilon) cos(epsilon)];

% v_inf in ECI reference frame
v_inf1 = T' * Deltavv1_minV;
r_inf = v_inf1 / norm(v_inf1);

a_H1 =-mu_t / dot(v_inf1,v_inf1);

% initialization
K = [0 0 1];



TH = linspace(0, 2*pi, 1000);
DeltaV1_ottSC3 = 10;
alphas = zeros(1, length(TH));
M = zeros(6, length(TH));
flags = zeros(1, length(TH));
id_acc = zeros(1, length(TH));



% First impulse
for j = 1:length(TH)

% Iterating through all initial positions
[rr_g,~] = par2car(a_SC1, e_SC1, i_SC1, OM_SC1, om_SC1, TH(j), mu_t);
rr_versore = rr_g/norm(rr_g);
cos_alfa = dot(r_inf, rr_g) / norm(rr_g);
alfa = acos(cos_alfa);
alphas(j) = alfa;

% Calculating e_H, th_inf, th_H
[SOL, F, Err] = fsolve(@(var) sist(var, alfa, norm(rr_g), a_H1), [1; pi/2; 0],opts);
flags(j) = Err;
e_H = SOL(1);
th_inf = wrapTo2Pi(SOL(2));
th_H = SOL(3);

    % Calculating orbital parameters
    h_H = cross(rr_versore, r_inf)/norm(cross(rr_versore, r_inf));
    N_H = cross(K, h_H) / norm(cross(K, h_H));
    i_H = acos(dot(h_H, K));
    if N_H(2) >= 0
        OM_H = acos(N_H(1));
    else
        OM_H = 2*pi - acos(N_H(1));
    end
    
    cosB = dot(N_H, rr_versore);
    sinB = dot(cross(N_H, rr_versore), h_H);
    
    Beta = atan2(sinB, cosB);
    om_H = Beta - th_H;

 if (th_inf > pi || th_inf < pi/2) || (a_H1*(1-e_H) < 6378)
      continue
 end

 % check fsolve convergence
  if Err ~= 1 
    M(:,j) = [NaN;NaN;NaN;NaN;NaN;NaN];
    id_acc(j) = NaN;
 else
    M(:,j) = [SOL; i_H; OM_H; om_H];
    id_acc(j) = j;
  end

end


% Extract acceptable points
id_acc = id_acc(~isnan(id_acc));
e_H_acc = M(1,id_acc);
th_inf_acc = M(2,id_acc);
th_H_acc = M(3,id_acc);
i_H_acc = M(4,id_acc);
OM_H_acc = M(5,id_acc);
om_H_acc = M(6,id_acc);
TH_acc = TH(id_acc);
alpha_acc = alphas(id_acc);

% Calculating all acceptable DeltaV
for j = 1:length(id_acc)
    p = a_H1*(1-e_H_acc(j)^2);

    % v on the hyperbola
    v_r = -sqrt(mu_t/p) * sin(th_H_acc(j));
    v_th = sqrt(mu_t/p) * ( e_H_acc(j) + cos(th_H_acc(j)));
   
    v_perif = [ v_r; 
                v_th; 
                0];


    R_perif_ECI= [ cos(OM_H_acc(j))*cos(om_H_acc(j)) - sin(OM_H_acc(j))*sin(om_H_acc(j))*cos(i_H_acc(j)),  -cos(OM_H_acc(j))*sin(om_H_acc(j)) - sin(OM_H_acc(j))*cos(om_H_acc(j))*cos(i_H_acc(j)),  sin(OM_H_acc(j))*sin(i_H_acc(j));
                   sin(OM_H_acc(j))*cos(om_H_acc(j)) + cos(OM_H_acc(j))*sin(om_H_acc(j))*cos(i_H_acc(j)),  -sin(OM_H_acc(j))*sin(om_H_acc(j)) + cos(OM_H_acc(j))*cos(om_H_acc(j))*cos(i_H_acc(j)), -cos(OM_H_acc(j))*sin(i_H_acc(j));
                   sin(om_H_acc(j))*sin(i_H_acc(j)),                                     cos(om_H_acc(j))*sin(i_H_acc(j)),                                      cos(i_H_acc(j)) ];
    
    v_ECI = R_perif_ECI * v_perif;
    
    
    % v on the initial orbit
    [~,v_ell] = par2car(a_SC1,e_SC1,i_SC1,OM_SC1,om_SC1,TH_acc(j),mu_t); 
    
    
    
    DeltaV1_SC3 = norm(v_ECI - v_ell);

    if DeltaV1_SC3 <= DeltaV1_ottSC3 % Saving optimum
        DeltaV1_ottSC3 = DeltaV1_SC3;
        DeltaV1_ottSC3_dir = v_ECI - v_ell;
        id_ott = j;

    end

end

% Second impulse
% Ecces hyperbolic speed
vinf_2 = norm(Deltavv2_minV);

% Arrival hyperbola to a circular orbit
a_H2 = - mu_A/vinf_2^2;

% Testing varius r_circ from the surface to r_SOI
rp_H2 = D/2 : 1 : rSOI_A;

e_H2 = 1 - rp_H2./a_H2;

v_fuga2 = sqrt(2*mu_A./rp_H2);

vp_H2 = sqrt(v_fuga2.^2 + vinf_2^2);

v_circ = sqrt(mu_A./rp_H2);

DeltaV_H2 = abs(v_circ - vp_H2);

Delta_impact_CP = - a_H2 .* sqrt(e_H2.^2 - 1);

% Total Delta V for the patched conics method (circular arrival orbit)
DeltaV_tot_CP = DeltaV1_ottSC3 + DeltaV_H2;



% Arrival hyperbola to an elliptical orbit

e_H2_ell = @(rp_H2) 1 - rp_H2./a_H2;

v_fuga2_ell = @(rp_H2) sqrt(2*mu_A./rp_H2);

vp_H2_ell = @(rp_H2) sqrt(v_fuga2_ell(rp_H2).^2 + vinf_2.^2);

rp_ellisse = @(rp_H2) rp_H2;

margin = 0.05; % km
ra_ellisse = rSOI_A ;

e_ellisse =  @(rp_H2) (ra_ellisse - rp_ellisse(rp_H2))./(ra_ellisse + rp_ellisse(rp_H2));

a_ellisse = @(rp_H2) rp_ellisse(rp_H2)  ./ (1 - e_ellisse(rp_H2));

p_ellisse = @(rp_H2) a_ellisse(rp_H2)  .* (1 - e_ellisse(rp_H2).^2);

vp_ellisse = @(rp_H2) sqrt(mu_A./p_ellisse(rp_H2)) .* (1 + e_ellisse(rp_H2));

DeltaV_H2_ell = @(rp_H2) abs(vp_ellisse(rp_H2)  - vp_H2_ell(rp_H2));

opts_ga = optimoptions('ga', 'Display', 'off', 'PlotFcn', []);
lb = D/2 + margin; % lower bounds for rp_H2
ub = rSOI_A; % upper bounds for rp_H2

rng(42); % rng seed for consistent results
[x0, ~] = ga(DeltaV_H2_ell, 1,[],[],[],[], lb, ub, [], opts_ga);
opts_fmincon = optimoptions('fmincon', 'Display', 'off', 'PlotFcn', []);

% Calculation of the minimum cost for the direct transfer between the two assigned orbits,
% the argument of periapsis of the corresponding transfer orbit,
% and the true anomalies of the departure and arrival positions
[OPT_ell, minV] = fmincon(DeltaV_H2_ell, x0, [], [], [], [], lb, ub, [], opts_fmincon);

rp_H2_ell = OPT_ell;
DeltaV_ell_CP = DeltaV1_ottSC3 + DeltaV_H2_ell(OPT_ell);

Delta_impact_ell_CP = -a_H2 * sqrt(e_H2_ell(rp_H2_ell)^2 -1);

% Plotting results

% Delta V plot 
figure('Name','Delta V comparison with plane change')
subplot(2,1,1)
hold on
grid on
box on
scatter(rp_H2-D/2, DeltaV_tot_CP,'b','filled')
yline(DeltaV_ell_CP,'k--','LineWidth',1.5)
xline(rSOI_A,'--r','LineWidth',1.5)
xlabel("Km from the surface [Km]",Interpreter="latex")
ylabel("$\Delta V$ [Km/s]",Interpreter="latex")
leg = legend('$\Delta V_{tot}$  circular orbit','$\Delta V_{tot}$ elliptical orbit','SOI','Interpreter','latex');
set(leg,'FontSize',22,'Location','best');

axis padded


subplot(2,1,2)
hold on
grid on
box on
scatter(rp_H2-D/2, Delta_impact_CP,'b','filled')
yline(Delta_impact_ell_CP,'k--','LineWidth',1.5)
xline(rSOI_A,'--r','LineWidth',1.5)
xlabel("Km from the surface [Km]",Interpreter="latex")
ylabel('$\Delta$ [Km]','Interpreter','latex')
leg = legend('$\Delta$  circular orbit','$\Delta$ elliptical orbit','SOI','Interpreter','latex');
set(leg,'FontSize',22,'Location','best');

axis padded




figure('Name','Non coplanar escape trajectory','Units', 'normalized', 'OuterPosition', [0 0 1 1])
hold on
earth_sphere('km');
grid on
box on
axis equal


% Initial orbit
plotOrbit(a_SC1, e_SC1, i_SC1, OM_SC1, om_SC1, 0, 2*pi, pi/1000, mu_t)

% Perigee
[rr_pSC1, ~] = par2car(a_SC1, e_SC1, i_SC1, OM_SC1, om_SC1, th_SC1, mu_t);
scatter3(rr_pSC1(1),rr_pSC1(2),rr_pSC1(3),'xk')

% Hyperbola
plotOrbitC(a_H1,e_H_acc(id_ott),i_H_acc(id_ott),OM_H_acc(id_ott),om_H_acc(id_ott),0, 2*pi, pi/10000,mu_t,6,[0.1 0.7 0.1],'-',th_H_acc(id_ott));

% Intersection
[rr_ott,~]= par2car(a_H1,e_H_acc(id_ott),i_H_acc(id_ott),OM_H_acc(id_ott),om_H_acc(id_ott),th_H_acc(id_ott),mu_t);
scatter3(rr_ott(1),rr_ott(2),rr_ott(3),'sk');


leg = legend('',"Initial orbit",'Perigee', "Hyperbolic orbit","Intersection with initial orbit", 'Interpreter','latex');
set(leg,'FontSize', 22)



Deltat_CP = TOF(a_SC1, e_SC1, 0, TH(id_ott), mu_t);

A9(1,:) = [0,a_SC1,e_SC1,i_SC1,OM_SC1,om_SC1,0,0];
A9(2,:) = [Deltat_CP,a_SC1,e_SC1,i_SC1,OM_SC1,om_SC1,TH(id_ott),DeltaV1_ottSC3];
A9(3,:) = [0,a_H1,e_H_acc(id_ott),i_H_acc(id_ott),OM_H_acc(id_ott),wrapTo2Pi(om_H_acc(id_ott)),th_H_acc(id_ott),0];

TAB9 = array2table(A9, 'VariableNames', ...
    {'Time [h]', 'a [AU]', 'e [-]', 'i [rad]', 'OM [rad]', 'om [rad]', 'theta [rad]', 'DeltaV [km/s]'});
disp(TAB9)


