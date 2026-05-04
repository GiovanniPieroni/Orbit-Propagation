%% 
% This part isn't required and is just a guess on how to further improve
% the transfer by modifing some of the orbit's fixed data.
%
% In particular, scenarios 2 and 3 are calculated as before. Then, the final 
% orbit for scenario 1 is modified to have the given a and e, while Ω, ω, and i 
% are computed from the h_H of the escape hyperbola. By doing so, the direct 
% transfer with minimum ΔV for scenario 1 must be recomputed, along with the 
% new escape hyperbola, since the previous one no longer intersects the updated 
% initial orbit.
%
% This approach ensures that the initial orbit and the escape hyperbola lie in 
% more similar planes. If repeated iteratively, it is possible to obtain a 
% coplanar escape hyperbola; however, the plane change in scenario 1 becomes 
% extremely costly. Performing this adjustment only once helps reduce the 
% escape ΔV.

%% SETUP
close
clear
clc

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

minV = minV_rel * DeltaV_MAX;
% Compute the time taken for the transfer
DeltaT_minV = DeltaT_fun(OPT);


% Compute the cost of the two impulses - modulus and vector
Deltavv1_minV = vv1_T(OPT(1), OPT(2), OPT(3)) - vv_i(OPT(1));
DeltaV1_minV = norm(vv1_T(OPT(1), OPT(2), OPT(3)) - vv_i(OPT(1)));
Deltavv2_minV = vv_f(OPT(2)) - vv2_T(OPT(1), OPT(2), OPT(3));
DeltaV2_minV = norm(vv_f(OPT(2)) - vv2_T(OPT(1), OPT(2), OPT(3)));

%% Arrival and escape trajectories (Do not run twice, always run after setup)

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
M = [];
alphas = [];
id_acc = [];
TH = linspace(0, 2*pi, 1000);
DeltaV1_ottSC3 =10;
H_H = [];


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
    H_H(:,j) = h_H;
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
 else
    M(:,j) = [SOL; i_H; OM_H; om_H];
    id_acc = [id_acc,j];
 end

end


% Extract acceptable points

e_H_acc = M(1,id_acc);
th_inf_acc = M(2,id_acc);
th_H_acc = M(3,id_acc);
i_H_acc = M(4,id_acc);
OM_H_acc = M(5,id_acc);
om_H_acc = M(6,id_acc);
TH_acc = TH(id_acc);
alpha_acc = alphas(id_acc);
H_H_acc = H_H(:,id_acc);

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
lb = D/2 + margin;  % upper and lower bounds
ub = rSOI_A;  

rng(42); % rng seed for consistent results
[x0, ~] = ga(DeltaV_H2_ell, 1,[],[],[],[], lb, ub, [], opts_ga);
opts_fmincon = optimoptions('fmincon', 'Display', 'off', 'PlotFcn', []);

% Minimizing DeltaV of arrival ellpiptical trajectory
[OPT_ell, minV] = fmincon(DeltaV_H2_ell, x0, [], [], [], [], lb, ub, [], opts_fmincon);

rp_H2_ell = OPT_ell;
DeltaV_ell_CP = DeltaV1_ottSC3 + DeltaV_H2_ell(OPT_ell);

Delta_impact_ell_CP = -a_H2 * sqrt(e_H2_ell(rp_H2_ell)^2 -1);




% recalculating SC1 final orbit to match h_H fo the hyperbola
h = H_H_acc(:,id_ott);
a_SC1 = 47899;
e_SC1 = 0.53109;
k = [0 0 1];

i_SC1 = acos(h(3)/norm(h));
N = cross(k,h);

if N(2) >= 0
    OM_SC1 = acos(N(1));
else
    OM_SC1 = 2*pi - acos(N(1));
end

om_SC1 = wrapTo2Pi(om_H_acc(id_ott));

% running again SC3 to recalculate the escape hyperbola with the new
% improved SC1 final orbit


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
M = [];
alphas = [];
id_acc = [];
TH = linspace(0, 2*pi, 1000);
DeltaV1_ottSC3 =10;
H_H = [];


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
    H_H(:,j) = h_H;
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
 else
    M(:,j) = [SOL; i_H; OM_H; om_H];
    id_acc = [id_acc,j];
 end

end


% Extract acceptable points

e_H_acc = M(1,id_acc);
th_inf_acc = M(2,id_acc);
th_H_acc = M(3,id_acc);
i_H_acc = M(4,id_acc);
OM_H_acc = M(5,id_acc);
om_H_acc = M(6,id_acc);
TH_acc = TH(id_acc);
alpha_acc = alphas(id_acc);
H_H_acc = H_H(:,id_acc);

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
lb = D/2 + margin;           % upper and lower bounds
ub = rSOI_A;  

rng(42); % rng seed for consistent results
[x0, ~] = ga(DeltaV_H2_ell, 1,[],[],[],[], lb, ub, [], opts_ga);
opts_fmincon = optimoptions('fmincon', 'Display', 'off', 'PlotFcn', []);

[OPT_ell, minV] = fmincon(DeltaV_H2_ell, x0, [], [], [], [], lb, ub, [], opts_fmincon);

rp_H2_ell = OPT_ell;
DeltaV_ell_CP = DeltaV1_ottSC3 + DeltaV_H2_ell(OPT_ell);

Delta_impact_ell_CP = -a_H2 * sqrt(e_H2_ell(rp_H2_ell)^2 -1);


% Recalculating direct transfer for SC1


format short g
mu_t = 398600;
mu_s = 132712440018;
Dati = importdata('DatiSC1-2025.mat','T');
SC_1 = table2array(Dati(22,:));

% Initial orbit data
a_i = SC_1(2);
e_i = SC_1(3);
i_i = SC_1(4);
OM_i = SC_1(5);
om_i = SC_1(6);
th_i = SC_1(7);
[rrr_i,vvv_i] = par2car(a_i,e_i,i_i,OM_i,om_i,th_i,mu_t);

% Final orbit data
a_f = a_SC1;
e_f = e_SC1;
i_f = i_SC1;
OM_f = OM_SC1;
om_f = om_SC1;
th_f = 0;
[rrr_f,vvv_f] = par2car(a_f,e_f,i_f,OM_f,om_f,th_f,mu_t);

% Arbitrary choices (discretization of th_1i,th_2f, om_t)
th_1i = 0:pi/100:2*pi;
th_2f = th_1i;          % NOTE: length(th_1i) == length(th_2f) == om_t
om_t = th_1i; 
dim = size(th_1i,2);

% Initialization
rr_i = zeros(3,length(th_1i));
vv_i = zeros(3,length(th_1i));
rr_f = zeros(3,length(th_2f));
vv_f = zeros(3,length(th_2f));
k = [0,0,1];        % for N computation
R_om_comb = zeros(3,1);
DeltaV_ott = 100;

% Compute rr and vv vectors at start and end
for h = 1:length(th_1i)
    [rr_i(:,h),vv_i(:,h)] = par2car(a_i,e_i,i_i,OM_i,om_i,th_1i(h),mu_t); % initial r
    rv_i = [rr_i;vv_i];
    [rr_f(:,h),vv_f(:,h)] = par2car(a_f,e_f,i_f,OM_f,om_f,th_2f(h),mu_t); % final r
    rv_f = [rr_f;vv_f];
    
    % Build all rotation matrices around om
    R_om = [cos(om_t(h)), sin(om_t(h)), 0;
            -sin(om_t(h)), cos(om_t(h)), 0;
            0, 0, 1];

    R_om_comb = [R_om_comb,R_om];
end

% Remove first column (zeros) from R_comb matrix
R_om_comb(:,1) = [];
% Decompose R_om_comb into 3D blocks
R_om_blocks = reshape(R_om_comb, 3, 3, []);

% Compute all combinations of rr_i and rr_f
comb = combinations(1:dim,1:dim); % all pairs between 1 and dim
comb = table2array(comb);
comb_i = comb(:,1)';
comb_f = comb(:,2)';

rr_comb = [rr_i(:,comb_i);rr_f(:,comb_f)];  % [rr_i,rr_f]
vv_comb = [vv_i(:,comb_i);vv_f(:,comb_f)];  % [vv_i,vv_f]

% Reshape vectors
rr_i_blocks = reshape(rr_comb(1:3,:),3,1,[]);
rr_f_blocks = reshape(rr_comb(4:6,:),3,1,[]);
vv_i_blocks = reshape(vv_comb(1:3,:),3,1,[]);
vv_f_blocks = reshape(vv_comb(4:6,:),3,1,[]);

for j = 1:size(comb,1)

    % Extract a combination
    rr_i_current = rr_i_blocks(:,:,j);
    rr_f_current = rr_f_blocks(:,:,j);
    vv_i_current = vv_i_blocks(:,:,j);
    vv_f_current = vv_f_blocks(:,:,j);

    rr_i_norm = norm(rr_i_current);
    rr_f_norm = norm(rr_f_current);

    h_t = cross(rr_i_current,rr_f_current)./norm(cross(rr_i_current,rr_f_current)); % transfer orbital plane
    i_t = acos(h_t(3)); % transfer inclination

    N = cross(k,h_t)./norm(cross(k,h_t)); % Line of nodes

    if N(2) >= 0    % RAAN
        OM_t = acos(N(1));
    else
        OM_t = 2*pi - acos(N(1));
    end

    R_OM = [cos(OM_t), sin(OM_t), 0;            % Omega rotation around k
            -sin(OM_t), cos(OM_t), 0;
            0, 0, 1];
    
    R_i = [1,      0,       0;              % Inclination rotation around i'
           0, cos(i_t), sin(i_t);
           0, -sin(i_t), cos(i_t)];

    R_comb = pagemtimes(R_om_blocks,R_i*R_OM); % Build all possible rotation matrices

    % Position vectors in perifocal system (for each om)
    rr_i_pf_current = pagemtimes(R_comb,rr_i_current);
    rr_f_pf_current = pagemtimes(R_comb,rr_f_current);

    cos_th_1t = rr_i_pf_current(1,:,:)./rr_i_norm;  % true anomaly at start of transfer
    sin_th_1t = rr_i_pf_current(2,:,:)./rr_i_norm;
    th_1t = atan2(sin_th_1t(:),cos_th_1t(:))';

    cos_th_2t = rr_f_pf_current(1,:,:)./rr_f_norm;  % true anomaly at end of transfer
    sin_th_2t = rr_f_pf_current(2,:,:)./rr_f_norm;
    th_2t = atan2(sin_th_2t(:),cos_th_2t(:))';

    e_t_current =  (rr_i_norm - rr_f_norm)./(rr_f_norm .*  cos(th_2t) -rr_i_norm .* cos(th_1t)); % eccentricity 
    % discard unacceptable eccentricities
    e_t_current(e_t_current<0) = NaN;

    a_t_current = (rr_i_norm.*(1+e_t_current.*cos(th_1t)))./(1-(e_t_current.^2)); % semi-major axis
    % discard unacceptable semi-major axes
    a_t_current(a_t_current< 100 + 6378) = NaN;

    % Compute vv_1t and vv_2t without par2car to avoid for-loop
    p = a_t_current.*(1-(e_t_current.^2));
    
    vv_calcolo_vv_1t = sqrt(mu_t./p) .* [-sin(th_1t); e_t_current + cos(th_1t); zeros(1,length(th_1t))];
    vv_calcolo_vv_1t = reshape(vv_calcolo_vv_1t,3,1,[]);
    vv_calcolo_vv_2t = sqrt(mu_t./p) .* [-sin(th_2t); e_t_current + cos(th_2t); zeros(1,length(th_2t))];
    vv_calcolo_vv_2t = reshape(vv_calcolo_vv_2t,3,1,[]);

    R_PFtoECI = pagetranspose(R_comb);

    vv_1t = pagemtimes(R_PFtoECI,vv_calcolo_vv_1t);
    vv_2t = pagemtimes(R_PFtoECI,vv_calcolo_vv_2t);

    % Compute DeltaV
    
    DeltaV1 = vecnorm(vv_1t-vv_i_current,2,1); % first deltaV
    DeltaV2 = vecnorm(vv_f_current-vv_2t,2,1); % second deltaV
    DeltaV = DeltaV1(:) + DeltaV2(:);
    [DeltaV_candidate,index] = min(DeltaV);

    if DeltaV_candidate < DeltaV_ott
        DeltaV_ott = DeltaV_candidate;
        i_t_ott = i_t; 
        OM_t_ott = OM_t;
        e_t_ott = e_t_current(index);
        a_t_ott = a_t_current(index);
        th_1t_ott = th_1t(index);
        th_2t_ott = th_2t(index);
        om_t_ott = om_t(index);
        th_1i_ott = th_1i(ceil(j/length(th_1t)));
        th_2f_ott = th_2f(rem(j,length(th_1t)));

        % direction of the two deltaVs
        DeltaV1_dir = (vv_1t-vv_i_current);
        DeltaV1_dir_ott = DeltaV1_dir(:,:,index);

        DeltaV2_dir = (vv_f_current-vv_2t);
        DeltaV2_dir_ott = DeltaV2_dir(:,:,index);
    end
end


TotalDeltaV = norm(DeltaV1_dir_ott)+norm(DeltaV2_dir_ott)+DeltaV_tot_CP(2)

th_1t_ott = wrapTo2Pi(th_1t_ott);
th_2t_ott = wrapTo2Pi(th_2t_ott);

DeltaT1 = TOF(a_i,e_i,th_i,th_1i_ott,mu_t);
DeltaT2 = TOF(a_t_ott,e_t_ott,th_1t_ott,th_2t_ott,mu_t);
DeltaT3 = TOF(a_f,e_f,th_2f_ott,TH_acc(id_ott),mu_t);
DeltaTtot = DeltaT1+DeltaT2+DeltaT3;

A12(1,:) = [0,a_i,e_i,i_i,OM_i,om_i,th_i,0];
A12(2,:) = [DeltaT1,a_i,e_i,i_i,OM_i,om_i,th_1i_ott,norm(DeltaV1_dir_ott)];
A12(3,:) = [0,a_t_ott,e_t_ott,i_t_ott,OM_t_ott,om_t_ott,th_1t_ott,0];
A12(4,:) = [DeltaT1+DeltaT2,a_t_ott,e_t_ott,i_t_ott,OM_t_ott,om_t_ott,th_2t_ott,norm(DeltaV2_dir_ott)];
A12(5,:) = [0,a_f,e_f,i_f,OM_f,om_f,th_2f_ott,0];
A12(6,:) = [DeltaTtot,a_f,e_f,i_f,OM_f,om_f,TH_acc(id_ott),DeltaV1_ottSC3];
A12(7,:) = [0,a_H1,e_H_acc(id_ott),i_H_acc(id_ott),OM_H_acc(id_ott),wrapTo2Pi(om_H_acc(id_ott)),th_H_acc(id_ott),0];