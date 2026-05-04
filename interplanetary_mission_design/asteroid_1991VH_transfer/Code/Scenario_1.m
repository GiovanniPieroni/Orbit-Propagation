%%%%%%%%%%%%%%%%%%%%%%%% General Setup - Run first %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                                                   %%%%%%%%%%%%
%%%%%%%%%%%%                    GROUP 22                       %%%%%%%%%%%%
%%%%%%%%%%%%                                                   %%%%%%%%%%%%
%%%%%%%%%%%%            Giovanni Pieroni - 10854347            %%%%%%%%%%%%
%%%%%%%%%%%%            Alessandro Ponti - 10901276            %%%%%%%%%%%%
%%%%%%%%%%%%            Jacopo Rotta     - 10696753            %%%%%%%%%%%%
%%%%%%%%%%%%    SC01 - Transfer from GTO to parking orbit      %%%%%%%%%%%%
%%%%%%%%%%%%               around the Earth                    %%%%%%%%%%%%
%%%%%%%%%%%%                                                   %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
clc


%% Initial and final parameters

format short g
mu_t = 398600; % Earth standard gravitational parameter [km^3/s^2]

Dati = importdata('DatiSC1-2025.mat','T');
SC_1 = table2array(Dati(22,:));

% Initial orbit parameters
a_i = SC_1(2);      % semi-major axis [km]
e_i = SC_1(3);      % eccentricity [-]
i_i = SC_1(4);      % inclination [rad]
OM_i = SC_1(5);     % RAAN [rad]
om_i = SC_1(6);     % pericenter argument [rad]
th_i = SC_1(7);     % true anomaly [rad]

% Inital orbit cartesian components
[rr_i, vv_i] = par2car(a_i, e_i, i_i, OM_i, om_i, th_i, mu_t); 


% Final orbit cartesian components
rr_f = (SC_1(8:10))';
vv_f = (SC_1(11:13))';

% Final orbit parameters
[a_f, e_f, i_f, OM_f, om_f, ~] = car2par(rr_f, vv_f, mu_t);
th_f = 0; 

DeltaVV = [];
DeltaTT = [];

I_s = 300 * 9.81;  % Specific gravimetric impulse [m/s]
m_0 = 300;         % Empty satellite's mass [kg]

%% Table A1 - Standard
clearvars -except mu_t a_i e_i i_i OM_i om_i th_i rr_i vv_i a_f e_f i_f + ...
OM_f om_f th_f rr_f vv_f DeltaVV DeltaTT A1 A2 A3 + ...
A4_grid A4_ga A5_grid A5_ga m_0 I_s
close all


% Plane change 
[DeltaVCP, omCP, thCP] = changeOrbitalPlane(a_i, e_i, i_i, OM_i, om_i, i_f, OM_f, mu_t);
DeltatCP = TOF(a_i, e_i, th_i, thCP, mu_t);

% Pericenter argument change 
[DeltaVom, th_iom, th_fom] = changePericenterArg(a_i, e_i, omCP, om_f+pi, mu_t);
Deltatom = TOF(a_i, e_i, thCP, th_iom(2), mu_t);

% Wait Bitangent
Deltat3 = TOF(a_i, e_i, th_fom(2), 0, mu_t);
% Bitangent transfer
[DeltaV1pp, DeltaV2pp, Deltat4, a_T, e_T] = bitangentTransfer(a_i, e_i, a_f, e_f, 'pp',  mu_t);

% Delta t
Deltat_A1 = DeltatCP + Deltatom + Deltat3 + Deltat4;
% Delta V
DeltaV_A1 = DeltaVCP + DeltaVom + DeltaV1pp + DeltaV2pp;


% Saving results
DeltaVV = [DeltaVV, DeltaV_A1]; 
DeltaTT = [DeltaTT, Deltat_A1];


% plotting results
figure('Name','A1 - Standard','Units', 'normalized', 'OuterPosition', [0 0 1 1])
hold on
grid on
box on
axis equal
earth_sphere('km');

% Initial orbit and position
plotOrbit(a_i, e_i, i_i, OM_i, om_i, 0, 2*pi, pi/10000, mu_t);
scatter3(rr_i(1),rr_i(2),rr_i(3),'*')

% Final orbit and position
plotOrbit(a_f, e_f, i_f, OM_f, om_f, 0, 2*pi, pi/10000, mu_t);
scatter3(rr_f(1),rr_f(2),rr_f(3),'*')

% Change of plane
plotOrbit(a_i, e_i, i_i, OM_i, om_i, th_i, thCP, pi/10000, mu_t, 3);
plotOrbit(a_i, e_i, i_f, OM_f, omCP, 0, 2*pi, pi/10000, mu_t);

% Change pericenter argument
plotOrbit(a_i, e_i, i_f, OM_f, om_f+pi, 0, 2*pi, pi/10000, mu_t);
plotOrbit(a_i, e_i, i_f, OM_f, omCP, thCP, th_iom(2), pi/10000, mu_t, 3);

% Bitangent transfer
plotOrbit(a_i, e_i, i_f, OM_f, om_f+pi, th_fom(2), 2*pi, pi/10000, mu_t, 3);
plotOrbit(a_T, e_T, i_f, OM_f, om_f+pi, 0, pi, pi/10000, mu_t, 3);


leg = legend('','Inital orbit','Inital point','Final orbit','Final point', 'Wait change of plane', '', '', 'Wait change of $\omega$', 'Wait 1st bitangent impulse', 'Wait 2nd bitangent impulse','Interpreter','latex');
set(leg, 'FontSize', 22);

m_p = m_0 *(exp(DeltaV_A1 * 1e3/I_s)-1); % Propellant mass

fprintf("================================================================\n"  + ...
        "Table: A1 \n" + ...
        "Optimal Delta V =  %.6f km/s \n" + ...
        "Respective Delta T = %.4f s = %.4f gg \n" + ...
        "Propellant mass (with Is = 300 s & M_empty = 300 kg): %.2f kg\n" + ...
        "================================================================\n\n", DeltaV_A1, Deltat_A1, Deltat_A1/(3600*24), m_p);


% Table A1
A1(1,:) = [0, a_i, e_i, i_i, OM_i, om_i, th_i, 0];
A1(2,:) = [DeltatCP, a_i, e_i, i_i, OM_i, om_i, thCP, DeltaVCP];
A1(3,:) = [DeltatCP, a_i, e_i, i_f, OM_f, omCP, thCP, 0];
A1(4,:) = [DeltatCP+Deltatom, a_i, e_i, i_f, OM_f, omCP, th_iom(2), DeltaVom];
A1(5,:) = [DeltatCP+Deltatom, a_i, e_i, i_f, OM_f, om_f+pi, th_fom(2), 0];
A1(6,:) = [DeltatCP+Deltatom+Deltat3, a_i, e_i, i_f, OM_f, om_f+pi, 0, DeltaV1pp];
A1(7,:) = [DeltatCP+Deltatom+Deltat3, a_T, e_T, i_f, OM_f, om_f+pi, 0, 0];
A1(8,:) = [DeltatCP+Deltatom+Deltat3+Deltat4, a_T, e_T, i_f, OM_f, om_f+pi, pi, DeltaV2pp];
A1(9,:) = [DeltatCP+Deltatom+Deltat3+Deltat4, a_f, e_f, i_f, OM_f, om_f, 0, 0];

TAB1 = array2table(A1, 'VariableNames', ...
    {'Time [s]', 'a [km]', 'e [-]', 'i [rad]', 'OM [rad]', 'om [rad]', 'theta [rad]', 'DeltaV [km/s]'});
disp(TAB1)


%% Table A2 - Circular orbit
clearvars -except mu_t a_i e_i i_i OM_i om_i th_i rr_i vv_i a_f e_f i_f + ...
OM_f om_f th_f rr_f vv_f DeltaVV DeltaTT A1 A2 A3 + ...
A4_grid A4_ga A5_grid A5_ga m_0 I_s
close all


% Wait first bitangent
Deltat1 = TOF(a_i,e_i,th_i,2*pi,mu_t);

% a2 = final orbit apoapsis radius
a_2 = a_f*(1+e_f);

% Bitangent transfer
[DeltaV1pa, DeltaV2pa, Deltat2, a_T, e_T] = bitangentTransfer(a_i, e_i, a_2,0, 'pa', mu_t);

% Change of plane
[DeltaVCP, omCP, thCP] = changeOrbitalPlane(a_2, 0, i_i, OM_i, om_i, i_f, OM_f, mu_t);

% Wait Change of plane
DeltatCP = TOF(a_2,0,pi,thCP+pi,mu_t);

% Wait last tangent impulse
Deltat4 = TOF(a_2,0,thCP-pi+omCP-om_f,pi,mu_t);

% Single tangent impulse
[DeltaV3pa, DeltaV4pa, ~, a_T2, e_T2] = bitangentTransfer(a_2, 0, a_f, e_f, 'pa', mu_t);

% Tangent Delta t
Deltat5 = TOF(a_f,e_f,pi,th_f,mu_t);

% Delta V 
DeltaV_A2 = DeltaVCP + DeltaV1pa+ DeltaV2pa + DeltaV3pa + DeltaV4pa;
% Delta t
Deltat_A2 = Deltat1 + Deltat2 + DeltatCP + Deltat4 + Deltat5;

% Saving results
DeltaVV = [DeltaVV, DeltaV_A2];
DeltaTT = [DeltaTT, Deltat_A2];


% plotting results
figure('Name','A2 - Circular orbit','Units', 'normalized', 'OuterPosition', [0 0 1 1])
hold on
grid on
box on
axis equal
earth_sphere('km');


% Initial orbit and position
plotOrbit(a_i, e_i, i_i, OM_i, om_i, 0, 2*pi, pi/10000, mu_t);
scatter3(rr_i(1),rr_i(2),rr_i(3),'b*');

% Final orbit and position
plotOrbit(a_f, e_f, i_f, OM_f, om_f, 0, 2*pi, pi/10000, mu_t);
scatter3(rr_f(1),rr_f(2),rr_f(3),'*');

% Wait first bitangent impulse 
plotOrbit(a_i,e_i,i_i,OM_i,om_i,th_i,2*pi,pi/10000,mu_t, 3);

% Bitangent transfer orbit
plotOrbit(a_T, e_T, i_i, OM_i, om_i, 0, 2*pi, pi/10000, mu_t);

% Wait second bitangent impulse
plotOrbit(a_T, e_T, i_i, OM_i, om_i, 0, pi, pi/10000, mu_t, 3)

% Orbit after bitangent
plotOrbit(a_2, 0, i_i, OM_i, om_i, 0, 2*pi, pi/10000, mu_t);

% Orbit after plane change
plotOrbit(a_2, 0, i_f, OM_f, omCP, 0, 2*pi, pi/10000, mu_t);

% Wait plane change
plotOrbit(a_2, 0, i_i, OM_i, om_i, pi, thCP+pi, pi/10000, mu_t, 3);

% Wait last tangent impulse
plotOrbit(a_2, 0, i_f, OM_f, om_f, thCP+pi+omCP-om_f, pi, pi/10000, mu_t, 3);

% Wait final position
plotOrbitC(a_f,e_f,i_f,OM_f, om_f, pi, 2*pi,pi/10000, mu_t, 3, 'b');


leg = legend('', 'Initial orbit', 'Initial position', 'Final orbit', 'Final position', 'Wait 1st bitangent 1st impulse', '', 'Wait 1st bitangent 2nd impulse', '', '', 'Wait change of plane', 'Wait tangent transfer', 'Wait final destination', 'Interpreter','latex');
set(leg,'FontSize',22);

m_p = m_0 *(exp(DeltaV_A2 * 1e3/I_s)-1); % Propellant mass

fprintf("================================================================\n"  + ...
        "Table: A2 \n" + ...
        "Optimal Delta V =  %.6f km/s \n" + ...
        "Respective Delta T = %.4f s = %.4f gg \n" + ...
        "Propellant mass (with Is = 300 s & M_empty = 300 kg): %.2f kg\n" + ...
        "================================================================\n\n", DeltaV_A2, Deltat_A2, Deltat_A2/(3600*24), m_p);


% Table A2
A2(1,:) = [0, a_i, e_i, i_i, OM_i, om_i, th_i, 0];
A2(2,:) = [Deltat1, a_i, e_i, i_i, OM_i, om_i, 0, DeltaV1pa];
A2(3,:) = [Deltat1, a_T, e_T, i_i, OM_i, om_i, 0, 0];
A2(4,:) = [Deltat1+Deltat2, a_T, e_T, i_i, OM_i, om_i, pi, DeltaV2pa];
A2(5,:) = [Deltat1+Deltat2, a_2, 0, i_i, OM_i, om_i, pi, 0];
A2(6,:) = [Deltat1+Deltat2+DeltatCP, a_2, 0, i_i, OM_i, om_i, thCP+pi, DeltaVCP];
A2(7,:) = [Deltat1+Deltat2+DeltatCP, a_2, 0, i_f, OM_f, omCP, thCP+pi, 0];
A2(8,:) = [Deltat1+Deltat2+DeltatCP, a_2, 0, i_f, OM_f, omCP, thCP+pi, DeltaV3pa];
A2(9,:) = [Deltat1+Deltat2+DeltatCP, a_T2, e_T2, i_f, OM_f, om_f, thCP-pi+omCP-om_f, 0];
A2(10,:) = [Deltat1+Deltat2+DeltatCP+Deltat4, a_T2, e_T2, i_f, OM_f, om_f, pi, DeltaV4pa];
A2(11,:) = [Deltat1+Deltat2+DeltatCP+Deltat4, a_f, e_f, i_f, OM_f, om_f, pi, 0];
A2(12,:) = [Deltat1+Deltat2+DeltatCP+Deltat4+Deltat5, a_f, e_f, i_f, OM_f, om_f, th_f, 0];

TAB2 = array2table(A2, 'VariableNames', ...
    {'Time [s]', 'a [km]', 'e [-]', 'i [rad]', 'OM [rad]', 'om [rad]', 'theta [rad]', 'DeltaV [km/s]'});
disp(TAB2)


%% Circular orbit Delta-V trend
clearvars -except  mu_t a_i e_i i_i OM_i om_i th_i rr_i vv_i a_f e_f i_f + ...
OM_f om_f th_f rr_f vv_f DeltaVV DeltaTT A1 A2 A3 + ...
A4_grid A4_ga A5_grid A5_ga m_0 I_s
close all



% 1 - Bitangent transfer
[DeltaV1, DeltaV2, DeltatBi1, ~, ~] = bitangentTransferFun_af(a_i, e_i, 0, 'pa', mu_t);

% 2 - Plane change (Obs: the choosen point isn't important since v is constant)
[DeltaVCP, ~, ~] = changeOrbitalPlaneFun_a(0, i_i, OM_i, om_i, i_f, OM_f, mu_t);

% 3 - Second bitangent transfer
[DeltaV3, DeltaV4, DeltatB2, a_T2, e_T2] = bitangentTransferFun_ai(0, a_f, e_f, 'ap', mu_t);

% Delta V calculation - function
DeltaV= @(a_circ) DeltaVCP(a_circ)  + DeltaV1(a_circ) + DeltaV2(a_circ)...
                    + DeltaV3(a_circ)  + DeltaV4(a_circ) ;

% Lower and upper boundaries
a_inf = a_i;
a_sup = a_f * 2;

figure
hold on
grid on
a_plot = linspace(a_inf, a_sup, 1000);
plot(a_plot/a_i, DeltaV(a_plot), 'LineWidth',1.25,'Color','b')
xline(a_i*(1 + e_i)/a_i,'r--')
xline(a_f*(1 + e_f)/a_i,'g--','LineWidth',1.5)
xlabel("$a_{\mathrm{circ}}/a_i$ [-]", 'Interpreter','latex', 'LineWidth',30)
ylabel("$\Delta V$ [km/s]", 'Interpreter','latex', 'LineWidth',30)
leg = legend('$\Delta V_{tot}$', '$r_{a,i}/a_i$', '$r_{a,f}/a_i$','interpreter','latex');
set(leg,'FontSize',22);
box on
axis padded
ax = gca;        
ax.FontSize = 18;  

%% Table A3 - Deltat_MIN (Inefficient)
clearvars -except mu_t a_i e_i i_i OM_i om_i th_i rr_i vv_i a_f e_f i_f + ...
OM_f om_f th_f rr_f vv_f DeltaVV DeltaTT A1 A2 A3 + ...
A4_grid A4_ga A5_grid A5_ga m_0 I_s
close all


% Plane change   
[DeltaVCP, omCP, thCP] = changeOrbitalPlane(a_i, e_i, i_i, OM_i, om_i, i_f, OM_f, mu_t);

thCP = thCP + pi; % Choosing closest point of plane change

% Ricalculating DeltaVCP
DeltaVCP = DeltaVCP / (1 + e_i * cos(thCP-pi)) * (1 + e_i * cos(thCP));

% Wait for plane change
DeltatCP = TOF(a_i, e_i, th_i, thCP, mu_t);

% Change of pericenter argument
[DeltaVom, th_iom, th_fom] = changePericenterArg(a_i, e_i, omCP, om_f+pi, mu_t);

% Wait pericenter argument change
Deltatom = TOF(a_i, e_i, thCP, th_iom(1), mu_t);

% Wait first bitangent impulse
Deltat3 = TOF(a_i, e_i, th_fom(1), 0, mu_t);
 
% Bitangent transfer 
[DeltaV1pp, DeltaV2pp, Deltat4, a_T, e_T] = bitangentTransfer(a_i, e_i, a_f, e_f, 'pp',  mu_t); 

% Deltat 
Deltat_A3 = DeltatCP + Deltatom + Deltat3 + Deltat4;
% DeltaV
DeltaV_A3 = DeltaVCP + DeltaVom + DeltaV1pp + DeltaV2pp;
 
% Saving results 
DeltaVV = [DeltaVV, DeltaV_A3]; 
DeltaTT = [DeltaTT, Deltat_A3];


% plotting results
figure('Name','A3 - Deltat_MIN (Inefficient)','Units', 'normalized', 'OuterPosition', [0 0 1 1])
hold on
grid on
box on
axis equal
earth_sphere('km');


% Initial orbit and position
plotOrbit(a_i, e_i, i_i, OM_i, om_i, 0, 2*pi, pi/10000, mu_t);
scatter3(rr_i(1),rr_i(2),rr_i(3),'b*');

% Final orbit and position
plotOrbit(a_f, e_f, i_f, OM_f, om_f, 0, 2*pi, pi/10000, mu_t);
scatter3(rr_f(1),rr_f(2),rr_f(3),'*');

% Wait change of plane
plotOrbit(a_i, e_i, i_i, OM_i, om_i, th_i, thCP, pi/10000, mu_t, 3);

% Orbit after plane change
plotOrbit(a_i, e_i, i_f, OM_f, omCP, 0, 2*pi, pi/10000, mu_t);

% Orbit after change of pericenter argument
plotOrbit(a_i, e_i, i_f, OM_f, om_f+pi, 0, 2*pi, pi/10000, mu_t);

% Wait for pericenter argument change
plotOrbit(a_i, e_i, i_f, OM_f, omCP, thCP, th_iom(1), pi/10000, mu_t, 3);

% Wait first bitangent
plotOrbit(a_i, e_i, i_f, OM_f, om_f+pi, th_fom(1), 0, pi/10000, mu_t, 3);

% Wait second bitangent
plotOrbit(a_T, e_T, i_f, OM_f, om_f+pi, 0, pi, pi/10000, mu_t, 3);



leg = legend('','Initial orbit', 'Initial position', 'Final orbit', 'Final position', 'Wait plane change', '', '', 'Wait change of $\omega$', 'Wait first bitangent impulse', 'Wait second bitangent impulse','Interpreter','latex');
set(leg, 'FontSize', 22)

m_p = m_0 *(exp(DeltaV_A3 * 1e3/I_s)-1); % Propellant mass

fprintf("================================================================\n"  + ...
        "Table: A3 \n" + ...
        "Optimal Delta T =  %.4f s = %.4f gg \n" + ...
        "Respective Delta V = %.6f km/s \n" + ...
        "Propellant mass (with Is = 300 s & M_empty = 300 kg): %.2f kg\n" + ...
        "================================================================\n\n", Deltat_A3, Deltat_A3/(3600*24) ,DeltaV_A3, m_p);

% Table A3
A3(1,:) = [0, a_i, e_i, i_i, OM_i, om_i, th_i, 0];
A3(2,:) = [DeltatCP, a_i, e_i, i_i, OM_i, om_i, thCP, DeltaVCP];
A3(3,:) = [DeltatCP, a_i, e_i, i_f, OM_f, omCP, thCP, 0];
A3(4,:) = [DeltatCP+Deltatom, a_i, e_i, i_f, OM_f, omCP, th_iom(1), DeltaVom];
A3(5,:) = [DeltatCP+Deltatom, a_i, e_i, i_f, OM_f, om_f+pi, th_fom(1), 0];
A3(6,:) = [DeltatCP+Deltatom+Deltat3, a_i, e_i, i_f, OM_f, om_f+pi, 0, DeltaV1pp];
A3(7,:) = [DeltatCP+Deltatom+Deltat3, a_T, e_T, i_f, OM_f, om_f+pi, 0, 0];
A3(8,:) = [DeltatCP+Deltatom+Deltat3+Deltat4, a_T, e_T, i_f, OM_f, om_f+pi, pi, DeltaV2pp];
A3(9,:) = [DeltatCP+Deltatom+Deltat3+Deltat4, a_f, e_f, i_f, OM_f, om_f, 0, 0];

TAB3 = array2table(A3, 'VariableNames', ...
    {'Time [s]', 'a [km]', 'e [-]', 'i [rad]', 'OM [rad]', 'om [rad]', 'theta [rad]', 'DeltaV [km/s]'});
disp(TAB3)

%% Table A4 - Direct transfer Delta t MIN (grid-search)
clearvars -except mu_t a_i e_i i_i OM_i om_i th_i rr_i vv_i a_f e_f i_f + ...
OM_f om_f th_f rr_f vv_f DeltaVV DeltaTT A1 A2 A3 + ...
A4_grid A4_ga A5_grid A5_ga m_0 I_s
close all


% fixed parameters (initial and final position)
th_1i = th_i;
th_2f = 0; 

% free parameter (pericenter argument)
om_t = 0:pi/1000:2*pi; 

dim = size(om_t,2);

DeltaV_margin = 10; % max Delta V



k = [0,0,1];
R_om_comb = zeros(3, 1); % Rotation matrix of om_t

% initial and final rr and vv
rv_i = [rr_i; vv_i]; 
rv_f = [rr_f; vv_f];



for h = 1 : dim
    
    % All possible rotation of om
    R_om = [cos(om_t(h)), sin(om_t(h)), 0;
            -sin(om_t(h)), cos(om_t(h)), 0;
            0, 0, 1];

    R_om_comb(:, end+1:end+3) = R_om;
end

R_om_comb(:,1) = [];
% matrices along the third dimention
R_om_blocks = reshape(R_om_comb, 3, 3, []);


    rr_i_norm = norm(rr_i);
    rr_f_norm = norm(rr_f);

    % transfer orbital plane
    h_t = cross(rr_i,rr_f)./norm(cross(rr_i,rr_f));
    % inclination
    i_t = acos(h_t(3));

    % line of nodes
    N = cross(k,h_t)./norm(cross(k,h_t));
    
    if N(2) >= 0    % RAAN
    OM_t = acos(N(1));
    else
        OM_t = 2*pi - acos(N(1));
    end

    % OM rotation
    R_OM = [cos(OM_t), sin(OM_t), 0;   
            -sin(OM_t), cos(OM_t), 0;
            0, 0, 1];
    % i rotation
    R_i = [1,      0,       0;             
           0, cos(i_t), sin(i_t);
           0, -sin(i_t), cos(i_t)];


    % all possible rotations
    R_comb = pagemtimes(R_om_blocks,R_i*R_OM); 
    
    % ECI --> PF
    rr_i_pf = pagemtimes(R_comb,rr_i);
    rr_f_pf = pagemtimes(R_comb,rr_f);

    % Initial true anomaly on transfer orbit
    cos_th_1t = rr_i_pf(1,:,:)./rr_i_norm; 
    sin_th_1t = rr_i_pf(2,:,:)./rr_i_norm;
    th_1t = atan2(sin_th_1t(:),cos_th_1t(:))';
    th_1t = wrapTo2Pi(th_1t);

    % Final true anomaly on transfer orbit
    cos_th_2t = rr_f_pf(1,:,:)./rr_f_norm;  
    sin_th_2t = rr_f_pf(2,:,:)./rr_f_norm;
    th_2t = atan2(sin_th_2t(:),cos_th_2t(:))';
    th_2t = wrapTo2Pi(th_2t);

    % Eccentricity
    e_t =  (rr_i_norm - rr_f_norm)./(rr_f_norm .*  cos(th_2t) -rr_i_norm .* cos(th_1t)); % eccentricità 

    % Considering only closed orbits 0<=e<1
    e_t(e_t>=1) = NaN;
    e_t(e_t<0) = NaN;
    
    % Semi-major axis
    a_t = (rr_i_norm.*(1+e_t.*cos(th_1t)))./(1-(e_t.^2));
    % Pericenter radius
    r_pericenter = a_t.*(1-e_t);

    % Pericenter above Karman line
    a_t(r_pericenter<= 6378 + 100) = NaN;

    % Calculating vv_1t, vv_2t without par2car to avoid for loop

    p = a_t.*(1-(e_t.^2));
    
    vv_calc_vv_1t = sqrt(mu_t./p) .* [-sin(th_1t); e_t + cos(th_1t); zeros(1,length(th_1t))];
    vv_calc_vv_1t = reshape(vv_calc_vv_1t,3,1,[]);
    vv_calc_vv_2t = sqrt(mu_t./p) .* [-sin(th_2t); e_t + cos(th_2t); zeros(1,length(th_2t))];
    vv_calc_vv_2t = reshape(vv_calc_vv_2t,3,1,[]);

    R_PFtoECI = pagetranspose(R_comb);

    vv_1t = pagemtimes(R_PFtoECI,vv_calc_vv_1t);
    vv_2t = pagemtimes(R_PFtoECI,vv_calc_vv_2t);

    % DeltaV
    
    DeltaV1 = vecnorm(vv_1t-vv_i,2,1);
    DeltaV2 = vecnorm(vv_f-vv_2t,2,1);
    DeltaV = DeltaV1(:) + DeltaV2(:);
    
    % TOF vectorialization
    cosE1 = (e_t + cos(th_1t)) ./ (1 + e_t .* cos(th_1t));
    cosE2 = (e_t + cos(th_2t)) ./ (1 + e_t .* cos(th_2t));
    sinE1 = sqrt(1 - e_t.^2) .* sin(th_1t) ./ (1 + e_t .* cos(th_1t));
    sinE2 = sqrt(1 - e_t.^2) .* sin(th_2t) ./ (1 + e_t .* cos(th_2t));
    
    E1 = wrapTo2Pi(atan2(sinE1(:), cosE1(:)))';
    E2 = wrapTo2Pi(atan2(sinE2(:), cosE2(:)))';
    
    DeltaT = sqrt(a_t.^3 ./ mu_t) .* (E2 - E1 - e_t .* (sinE2 - sinE1)) + 2.*pi.*sqrt(a_t.^3 ./ mu_t).*(th_2t<th_1t);

    % Removing transfers with DeltaV > DeltaV_margin
    DeltaT(DeltaV>DeltaV_margin) = NaN;
    DeltaT(DeltaT<0) = NaN;
    
    % Find Delta V MIN
    [DeltaV_ott,index] = min(DeltaV);

    % Saving optimum transfer
    i_t_ott = i_t; 
    OM_t_ott = OM_t;
    e_t_ott = e_t(index);
    a_t_ott = a_t(index);
    th_1t_ott = th_1t(index);
    th_2t_ott = th_2t(index);
    om_t_ott = om_t(index);
    th_1i_ott = th_1i(ceil(1/length(th_1t)));
    th_2f_ott = th_2f(rem(1,length(th_1t)));
    DeltaT_ott = DeltaT(index);
    DeltaV1_ott = vv_1t(:,:,index)-vv_i;
    DeltaV2_ott = vv_f-vv_2t(:,:,index);

    % DeltaV direction
    DeltaV1_dir = (vv_1t-vv_i);
    DeltaV1_dir_ott = DeltaV1_dir(:,:,index);

    DeltaV2_dir = (vv_f-vv_2t);
    DeltaV2_dir_ott = DeltaV2_dir(:,:,index);


% Rounding angles
th_1t_ott = wrapTo2Pi(th_1t_ott);
th_2t_ott = wrapTo2Pi(th_2t_ott);
th_1i_ott = wrapTo2Pi(th_1i_ott);
th_2f_ott = wrapTo2Pi(th_2f_ott);

% Saving results
DeltaVV = [DeltaVV,DeltaV_ott];
DeltaTT = [DeltaTT,DeltaT_ott];


% Plotting results

figure('Name','A4 - Direct transfer Delta t MIN (grid-search)','Units', 'normalized', 'OuterPosition', [0 0 1 1])
hold on
grid on
box on
axis equal
earth_sphere('km');

% Initial orbit and position
plotOrbit(a_i, e_i, i_i, OM_i, om_i, 0, 2*pi, pi/10000, mu_t);
scatter3(rr_i(1),rr_i(2),rr_i(3),'b*');

% Final orbit and position
plotOrbit(a_f, e_f, i_f, OM_f, om_f, 0, 2*pi, pi/10000, mu_t);
scatter3(rr_f(1),rr_f(2),rr_f(3),'*');

% Transfer orbit
plotOrbit(a_t_ott,e_t_ott,i_t_ott,OM_t_ott,om_t_ott,0,2*pi,pi/1000,mu_t);

% Transfer trajectory
plotOrbit(a_t_ott,e_t_ott,i_t_ott,OM_t_ott,om_t_ott,th_1t_ott,th_2t_ott,pi/1000,mu_t, 3);

leg = legend('','Inital orbit','Inital position','Final orbit','Final position','','Transfer trajectory','Interpreter','latex');
set(leg,'FontSize',22)


m_p = m_0 *(exp(DeltaV_ott * 1e3/I_s)-1); % Propellant mass

fprintf("================================================================\n"  + ...
        "Table: A4 with grid-search \n" + ...
        "Optimal Delta T = %.4f s = %.4f gg  \n" + ...
        "Respective Delta V = %.6f km/s \n" + ...
        "Propellant mass (with Is = 300 s & M_empty = 300 kg): %.2f kg\n" + ...
        "================================================================\n\n", DeltaT_ott, DeltaT_ott/(3600*24),DeltaV_ott, m_p);

% Table A4 grid-search
A4_grid(1,:) = [0,a_i,e_i,i_i,OM_i,om_i,th_i,norm(DeltaV1_ott)];
A4_grid(2,:) = [0,a_t_ott,e_t_ott,i_t_ott,OM_t_ott,om_t_ott,th_1t_ott,0];
A4_grid(3,:) = [DeltaT_ott,a_t_ott,e_t_ott,i_t_ott,OM_t_ott,om_t_ott,th_2t_ott,norm(DeltaV2_ott)];
A4_grid(4,:) = [DeltaT_ott,a_f,e_f,i_i,OM_i,om_i,0,0];

TAB4_grid = array2table(A4_grid, 'VariableNames', ...
    {'Time [s]', 'a [km]', 'e [-]', 'i [rad]', 'OM [rad]', 'om [rad]', 'theta [rad]', 'DeltaV [km/s]'});
disp(TAB4_grid)

%% Table A4 - Direct transfer Delta t MIN (ga)
clearvars -except mu_t a_i e_i i_i OM_i om_i th_i rr_i vv_i a_f e_f i_f + ...
OM_f om_f th_f rr_f vv_f DeltaVV DeltaTT A1 A2 A3 + ...
A4_grid A4_ga A5_grid A5_ga m_0 I_s
close all



% Transfer orbital plane
h_T = cross(rr_i, rr_f)/norm(cross(rr_i, rr_f));
% Inclination
i_T = acos(h_T(3));

% X axis
I = [1 0 0];
% Y axis
J = [0 1 0];
% Z axis
K = [0 0 1];

% Line of nodes
N = cross(K, h_T) / norm(cross(K, h_T));

% RAAN
if N(2)>0
OM_T = acos(N(1));
else
OM_T = 2 * pi - acos(N(1));
end

% Rotation matrices
R_OM_T = [cos(OM_T), sin(OM_T), 0;
         -sin(OM_T), cos(OM_T), 0;
             0,          0,     1];

R_i_T = [1,       0,        0;     
         0,    cos(i_T), sin(i_T);
         0,   -sin(i_T), cos(i_T)];

R_om_T = @(om_T) [cos(om_T), sin(om_T), 0;
                 -sin(om_T), cos(om_T), 0;
                     0,         0,      1];

% ECI --> PF
T_ECItoPF = @(om_T) R_om_T(om_T) * R_i_T * R_OM_T;

rr_i_PF = @(om_T) T_ECItoPF(om_T) * rr_i;
rr_f_PF = @(om_T) T_ECItoPF(om_T) * rr_f;

% Inital and final true anomalies on transfer orbit
CosTh1_T = @(om_T) dot(rr_i_PF(om_T), I) / norm(rr_i);
SinTh1_T = @(om_T) dot(rr_i_PF(om_T), J) / norm(rr_i);
    
CosTh2_T = @(om_T) dot(rr_f_PF(om_T), I) / norm(rr_f);
SinTh2_T = @(om_T) dot(rr_f_PF(om_T), J) / norm(rr_f);

Th1_T = @(om_T) wrapTo2Pi(atan2(SinTh1_T(om_T), CosTh1_T(om_T)));

Th2_T = @(om_T) wrapTo2Pi(atan2(SinTh2_T(om_T), CosTh2_T(om_T)));


% Eccentricity
e_T = @(om_T) (norm(rr_f) - norm(rr_i)) / (CosTh1_T(om_T) * norm(rr_i) - CosTh2_T(om_T) * norm(rr_f));

% Semi-major axis
a_T = @(om_T) norm(rr_i) * (1 + e_T(om_T) * CosTh1_T(om_T)) / (1 - (e_T(om_T))^2);
    
[~, vv1_T] = par2carFUN1(a_T, e_T, i_T, OM_T, Th1_T, mu_t);
[~, vv2_T] = par2carFUN1(a_T, e_T, i_T, OM_T, Th2_T, mu_t);

% Delta V
DeltaV = @(om_T) norm(vv1_T(om_T) - vv_i) + norm(vv_f - vv2_T(om_T));
    
opts = optimoptions('ga', 'Display', 'off', 'PlotFcn', []);
lb = 0;        % om_T lower bound
ub = 2*pi;     % om_T upper bound

% transfer orbit must be closed --> 0<=e_T<1
cond = @(x) deal([e_T(x) - 1;-e_T(x)], []);

% Minimizing DeltaV(om_T)
[om_opt_T, minV] = ga(DeltaV, 1, [], [], [], [], lb, ub, cond, opts);

% DeltaV1, DeltaV2
DeltaV1_minV = norm(vv1_T(om_opt_T) - vv_i);
DeltaV2_minV = norm(vv_f - vv2_T(om_opt_T));

% Deltat(om_T)
DeltaT = @(om) customTOF(om, a_T, e_T, Th1_T, Th2_T, mu_t);

% Calculate Delta t
DeltaTga = DeltaT(om_opt_T);

% Saving results
DeltaVV = [DeltaVV,minV];
DeltaTT = [DeltaTT,DeltaTga];


figure('Name','A4 - Direct transfer Delta t MIN (ga)','Units', 'normalized', 'OuterPosition', [0 0 1 1])
hold on
grid on
box on
axis equal
earth_sphere('km');

% Initial orbit and position
plotOrbit(a_i, e_i, i_i, OM_i, om_i, 0, 2*pi, pi/10000, mu_t);
scatter3(rr_i(1),rr_i(2),rr_i(3),'b*');

% Final orbit and position
plotOrbit(a_f, e_f, i_f, OM_f, om_f, 0, 2*pi, pi/10000, mu_t);
scatter3(rr_f(1),rr_f(2),rr_f(3),'*');

% Transfer orbit
plotOrbit(a_T(om_opt_T), e_T(om_opt_T), i_T, OM_T, om_opt_T, 0, 2*pi, pi/10000, mu_t)
% Transfer trajectory
plotOrbit(a_T(om_opt_T), e_T(om_opt_T), i_T, OM_T, om_opt_T, Th1_T(om_opt_T), Th2_T(om_opt_T), pi/10000, mu_t, 3)

leg = legend('', 'Inital orbit', 'Initial position', 'Final orbit', 'Final position','', 'Transfer trajectory','Interpreter','latex');
set(leg, 'FontSize',22)


m_p = m_0 *(exp(minV * 1e3/I_s)-1); % Propellant mass

fprintf("================================================================\n"  + ...
        "Table: A4 with ga \n" + ...
        "Optimal Delta T =  %.4f s = %.4f gg \n" + ...
        "Respective Delta V = %.6f km/s \n" + ...
        "Propellant mass (with Is = 300 s & M_empty = 300 kg): %.2f kg\n" + ...
        "================================================================\n\n", DeltaTga, DeltaTga/(3600*24), minV, m_p);


% Table A4 - Direct transfer Delta t MIN (ga)
A4_ga(1,:) = [0, a_i, e_i, i_i, OM_i, om_i, th_i, DeltaV1_minV];
A4_ga(2,:) = [0, a_T(om_opt_T), e_T(om_opt_T), i_T, OM_T, om_opt_T, Th1_T(om_opt_T), 0];
A4_ga(3,:) = [DeltaTga, a_T(om_opt_T), e_T(om_opt_T), i_T, OM_T, om_opt_T, Th2_T(om_opt_T), DeltaV2_minV];
A4_ga(4,:) = [DeltaTga, a_f, e_f, i_f, OM_f, om_f, th_f, 0];

TAB4_ga = array2table(A4_ga, 'VariableNames', ...
    {'Time [s]', 'a [km]', 'e [-]', 'i [rad]', 'OM [rad]', 'om [rad]', 'theta [rad]', 'DeltaV [km/s]'});
disp(TAB4_ga)


%% Table A5 - Direct transfer Delta V MIN (grid-search)
clearvars -except mu_t a_i e_i i_i OM_i om_i th_i rr_i vv_i a_f e_f i_f + ...
OM_f om_f th_f rr_f vv_f DeltaVV DeltaTT A1 A2 A3 + ...
A4_grid A4_ga A5_grid A5_ga m_0 I_s
close all


% Free parameters
th_1i = 0:pi/100:2*pi;
th_2f = th_1i;        
om_t = th_1i; 
dim = size(th_1i,2);

% Initial optimum
DeltaV_ott = 100;


% Initialization of all possible initial and final positions
rrr_i = zeros(3,length(th_1i));
vvv_i = zeros(3,length(th_1i));
rrr_f = zeros(3,length(th_2f));
vvv_f = zeros(3,length(th_2f));
k = [0,0,1];        
R_om_comb = zeros(3, 1); % Rotation matrix of om_t


% Calculating all possibile position and velocities
for h = 1:length(th_1i)
    [rrr_i(:,h),vvv_i(:,h)] = par2car(a_i,e_i,i_i,OM_i,om_i,th_1i(h),mu_t); % r iniziale
    rv_i = [rrr_i;vvv_i];
    [rrr_f(:,h),vvv_f(:,h)] = par2car(a_f,e_f,i_f,OM_f,om_f,th_2f(h),mu_t); % r finale
    rv_f = [rrr_f;vvv_f];
    
    % All possible om rotations
    R_om = [cos(om_t(h)), sin(om_t(h)), 0;
            -sin(om_t(h)), cos(om_t(h)), 0;
            0, 0, 1];

    R_om_comb(:, end+1:end+3) = R_om;
end

R_om_comb(:,1) = [];
% matrices along the third dimention
R_om_blocks = reshape(R_om_comb, 3, 3, []);


% Calculating all possibile couples of rrr_i and rrr_f

comb = combinations(1:dim,1:dim); 
comb = table2array(comb);
comb_i = comb(:,1)';
comb_f = comb(:,2)';

% Each column is a combination
rr_comb = [rrr_i(:,comb_i);rrr_f(:,comb_f)];  
vv_comb = [vvv_i(:,comb_i);vvv_f(:,comb_f)];  

% Vectors along the third dimention
rr_i_blocks = reshape(rr_comb(1:3,:),3,1,[]);
rr_f_blocks = reshape(rr_comb(4:6,:),3,1,[]);
vv_i_blocks = reshape(vv_comb(1:3,:),3,1,[]);
vv_f_blocks = reshape(vv_comb(4:6,:),3,1,[]);



for j = 1:size(comb,1)

    % Extracting a combination
    rr_i_current = rr_i_blocks(:,:,j);
    rr_f_current = rr_f_blocks(:,:,j);
    vv_i_current = vv_i_blocks(:,:,j);
    vv_f_current = vv_f_blocks(:,:,j);

    rr_i_norm = norm(rr_i_current);
    rr_f_norm = norm(rr_f_current);

    % Transfer orbital plane
    h_t = cross(rr_i_current,rr_f_current)./norm(cross(rr_i_current,rr_f_current));
    % inclination
    i_t = acos(h_t(3)); 

    % Line of nodes
    N = cross(k,h_t)./norm(cross(k,h_t)); 
    
    if N(2) >= 0    % RAAN
    OM_t = acos(N(1));
    else
        OM_t = 2*pi - acos(N(1));
    end

    % OM rotation
    R_OM = [cos(OM_t), sin(OM_t), 0;            
            -sin(OM_t), cos(OM_t), 0;
            0, 0, 1];
    
    % i rotation
    R_i = [1,      0,       0;              
           0, cos(i_t), sin(i_t);
           0, -sin(i_t), cos(i_t)];

    % All possibile rotations once OM and i are fixed 
    R_comb = pagemtimes(R_om_blocks,R_i*R_OM);
    
    % ECI --> PF for each om_t
    rr_i_pf_current = pagemtimes(R_comb,rr_i_current);
    rr_f_pf_current = pagemtimes(R_comb,rr_f_current);

    % Initial and final true anomaly on the transfer orbit (for each om_t)
    cos_th_1t = rr_i_pf_current(1,:,:)./rr_i_norm;  
    sin_th_1t = rr_i_pf_current(2,:,:)./rr_i_norm;
    th_1t = atan2(sin_th_1t(:),cos_th_1t(:))';

    cos_th_2t = rr_f_pf_current(1,:,:)./rr_f_norm; 
    sin_th_2t = rr_f_pf_current(2,:,:)./rr_f_norm;
    th_2t = atan2(sin_th_2t(:),cos_th_2t(:))';

    % Eccentricity 
    e_t_current =  (rr_i_norm - rr_f_norm)./(rr_f_norm .*  cos(th_2t) -rr_i_norm .* cos(th_1t)); % eccentricità 
    % Transfer orbit must be closed --> 0<=e<1
    e_t_current(e_t_current<0) = NaN;
    e_t_current(e_t_current>=1) = NaN;

    % Semi-major axis
    a_t_current = (rr_i_norm.*(1+e_t_current.*cos(th_1t)))./(1-(e_t_current.^2)); % semiasse maggiore
    r_pericenter = a_t_current.*(1-e_t_current);
    % Pericenter must be above Karman line
    a_t_current(r_pericenter< 6378 + 100) = NaN;

    % Calculating vv_1t, vv_2t without par2car to avoid for loop

    p = a_t_current.*(1-(e_t_current.^2));
    
    vv_calc_vv_1t = sqrt(mu_t./p) .* [-sin(th_1t); e_t_current + cos(th_1t); zeros(1,length(th_1t))];
    vv_calc_vv_1t = reshape(vv_calc_vv_1t,3,1,[]);
    vv_calc_vv_2t = sqrt(mu_t./p) .* [-sin(th_2t); e_t_current + cos(th_2t); zeros(1,length(th_2t))];
    vv_calc_vv_2t = reshape(vv_calc_vv_2t,3,1,[]);

    R_PFtoECI = pagetranspose(R_comb);

    vv_1t = pagemtimes(R_PFtoECI,vv_calc_vv_1t);
    vv_2t = pagemtimes(R_PFtoECI,vv_calc_vv_2t);

    %DeltaV
    DeltaV1 = vecnorm(vv_1t-vv_i_current,2,1);
    DeltaV2 = vecnorm(vv_f_current-vv_2t,2,1);
    DeltaV = DeltaV1(:) + DeltaV2(:);

    % If statement is applied only to the best of each rrr_i, rrr_f couple
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

        % DeltaV direction
        DeltaV1_dir = (vv_1t-vv_i_current);
        DeltaV1_dir_ott = DeltaV1_dir(:,:,index);

        DeltaV2_dir = (vv_f_current-vv_2t);
        DeltaV2_dir_ott = DeltaV2_dir(:,:,index);
        

    end

end

% Rounding angles
th_1i_ott = wrapTo2Pi(th_1i_ott);
th_1t_ott = wrapTo2Pi(th_1t_ott);
th_2t_ott = wrapTo2Pi(th_2t_ott);
th_2f_ott = wrapTo2Pi(th_2f_ott);

% Delta t
DeltaT1 = TOF(a_i,e_i,th_i,th_1i_ott,mu_t);
DeltaT2 = TOF(a_t_ott,e_t_ott,th_1t_ott,th_2t_ott,mu_t);
DeltaT3 = TOF(a_f,e_f,th_2f_ott,th_f,mu_t);
DeltaTtot = DeltaT1+DeltaT2+DeltaT3;

% Saving results
DeltaVV = [DeltaVV, DeltaV_ott];
DeltaTT = [DeltaTT, DeltaTtot];


% Plotting resuts
figure('Name','A5 - Direct transfer Delta V MIN (grid-search)','Units', 'normalized', 'OuterPosition', [0 0 1 1]);
hold on
grid on
box on
axis equal
earth_sphere('km');

% Initial orbit and position
plotOrbit(a_i, e_i, i_i, OM_i, om_i, 0, 2*pi, pi/10000, mu_t);
scatter3(rr_i(1),rr_i(2),rr_i(3),'b*');

% Final orbit and position
plotOrbit(a_f, e_f, i_f, OM_f, om_f, 0, 2*pi, pi/10000, mu_t);
scatter3(rr_f(1),rr_f(2),rr_f(3),'*');

%Transfer orbit
plotOrbit(a_t_ott,e_t_ott,i_t_ott,OM_t_ott,om_t_ott,0,2*pi,pi/1000,mu_t);

% Wait first impulse
plotOrbit(a_i,e_i,i_i,OM_i,om_i,th_i,th_1i_ott,pi/1000,mu_t, 3);
% Wait second impulse
plotOrbit(a_t_ott,e_t_ott,i_t_ott,OM_t_ott,om_t_ott,th_1t_ott,th_2t_ott,pi/1000,mu_t, 3);
%Wait final position
plotOrbit(a_f, e_f,i_f,OM_f,om_f,th_2f_ott,th_f+2*pi,pi/1000,mu_t, 3);


leg = legend('','Inital orbit','Inital position','Final orbit','Final position','','Wait 1st impulse','Transfer trajectory','Wait final destination','Interpreter','latex');
set(leg, 'FontSize', 22);


m_p = m_0 *(exp(DeltaV_ott * 1e3/I_s)-1); % Propellant mass

fprintf("================================================================\n"  + ...
        "Table: A5 with grid-search \n" + ...
        "Optimal Delta V =  %.6f km/s \n" + ...
        "Respective Delta T = %.4f s = %.4f gg \n" + ...
        "Propellant mass (with Is = 300 s & M_empty = 300 kg): %.2f kg\n" + ...
        "================================================================\n\n", DeltaV_ott, DeltaTtot, DeltaTtot/(3600*24), m_p);

% Table A5 - Direct transfer Delta V MIN (grid-search)
A5_grid(1,:) = [0,a_i,e_i,i_i,OM_i,om_i,th_i,0];
A5_grid(2,:) = [DeltaT1,a_i,e_i,i_i,OM_i,om_i,th_1i_ott,norm(DeltaV1_dir_ott)];
A5_grid(3,:) = [DeltaT1,a_t_ott,e_t_ott,i_t_ott,OM_t_ott,om_t_ott,th_1t_ott,0];
A5_grid(4,:) = [DeltaT1+DeltaT2,a_t_ott,e_t_ott,i_t_ott,OM_t_ott,om_t_ott,th_2t_ott,norm(DeltaV2_dir_ott)];
A5_grid(5,:) = [DeltaT1+DeltaT2,a_f,e_f,i_f,OM_f,om_f,th_2f_ott,0];
A5_grid(6,:) = [DeltaTtot,a_f,e_f,i_f,OM_f,om_f,th_f,0];

TAB5_grid = array2table(A5_grid, 'VariableNames', ...
    {'Time [s]', 'a [km]', 'e [-]', 'i [rad]', 'OM [rad]', 'om [rad]', 'theta [rad]', 'DeltaV [km/s]'});
disp(TAB5_grid)

%% Table A5 - Direct transfer Delta V MIN (ga)
clearvars -except mu_t a_i e_i i_i OM_i om_i th_i rr_i vv_i a_f e_f i_f + ...
OM_f om_f th_f rr_f vv_f DeltaVV DeltaTT A1 A2 A3 + ...
A4_grid A4_ga A5_grid A5_ga m_0 I_s
close all


% position and velocity @(true anomaly)
[RR_i, VV_i] = par2carFUN1th(a_i, e_i, i_i, OM_i, om_i, mu_t);
[RR_f, VV_f] = par2carFUN1th(a_f, e_f, i_f, OM_f, om_f, mu_t);

% X axis
I = [1 0 0];
% Y axis
J = [0 1 0];
% Z axis
K = [0 0 1];
   

% Transfer orbital plane @(initial and final true anomaly)
h_T = @(TH_i, TH_f) cross(RR_i(TH_i), RR_f(TH_f))/norm(cross(RR_i(TH_i), RR_f(TH_f)));
% Inclination @(initial and final true anomaly)
i_T = @(TH_i, TH_f) wrapTo2Pi(acos(dot(h_T(TH_i,TH_f),K)));


% Line of nodes @(initial and final true anomaly)
N = @(TH_i,TH_f) cross(K, h_T(TH_i,TH_f)) / norm(cross(K, h_T(TH_i,TH_f)));

% RAAN @(initial and final true anomaly)
OM_T = @(TH_i,TH_f) wrapTo2Pi(atan2(dot(N(TH_i,TH_f), J), dot(N(TH_i,TH_f), I)));
        

R_OM_T = @(TH_i,TH_f) [cos(OM_T(TH_i,TH_f)), sin(OM_T(TH_i,TH_f)), 0;
                      -sin(OM_T(TH_i,TH_f)), cos(OM_T(TH_i,TH_f)), 0;
                                 0,               0,               1];

R_i_T = @(TH_i,TH_f) [1,            0,                   0;     
                      0,    cos(i_T(TH_i,TH_f)), sin(i_T(TH_i,TH_f));
                      0,   -sin(i_T(TH_i,TH_f)), cos(i_T(TH_i,TH_f))];

R_om_T = @(om_T) [cos(om_T), sin(om_T), 0;
                 -sin(om_T), cos(om_T), 0;
                     0,         0,      1];

% Rotation matix ECI --> PF @(initial and final true , pericenter argument)
T_ECItoPF = @(TH_i,TH_f,om_T) R_om_T(om_T) * R_i_T(TH_i,TH_f) * R_OM_T(TH_i,TH_f);

RR_i_PF = @(TH_i,TH_f,om_T) T_ECItoPF(TH_i,TH_f,om_T) * RR_i(TH_i);
RR_f_PF = @(TH_i,TH_f,om_T) T_ECItoPF(TH_i,TH_f,om_T) * RR_f(TH_f);

% transfer orbit initial and final true anomalies @(initial and final true , pericenter argument)
CosTh1_T = @(TH_i,TH_f,om_T) dot(RR_i_PF(TH_i,TH_f,om_T),I) / norm(RR_i(TH_i));
SinTh1_T = @(TH_i,TH_f,om_T) dot(RR_i_PF(TH_i,TH_f,om_T),J) / norm(RR_i(TH_i));

CosTh2_T = @(TH_i,TH_f,om_T) dot(RR_f_PF(TH_i,TH_f,om_T),I) / norm(RR_f(TH_f));
SinTh2_T = @(TH_i,TH_f,om_T) dot(RR_f_PF(TH_i,TH_f,om_T),J) / norm(RR_f(TH_f));

Th1_T = @(TH_i,TH_f,om_T) wrapTo2Pi(atan2(SinTh1_T(TH_i,TH_f,om_T), CosTh1_T(TH_i,TH_f,om_T)));

Th2_T = @(TH_i,TH_f,om_T) wrapTo2Pi(atan2(SinTh2_T(TH_i,TH_f,om_T), CosTh2_T(TH_i,TH_f,om_T)));

% Eccentricity @(initial and final true , pericenter argument)
e_T = @(TH_i,TH_f,om_T) (norm(RR_f(TH_f)) - norm(RR_i(TH_i))) / (CosTh1_T(TH_i,TH_f,om_T) * norm(RR_i(TH_i)) - CosTh2_T(TH_i,TH_f,om_T) * norm(RR_f(TH_f)));

% Semi-major axis @(initial and final true , pericenter argument)
a_T = @(TH_i,TH_f,om_T) norm(RR_i(TH_i)) * (1 + e_T(TH_i,TH_f,om_T) * CosTh1_T(TH_i,TH_f,om_T)) / (1 - (e_T(TH_i,TH_f,om_T))^2);

[rr1_T, vv1_T] = par2carFUN3(a_T, e_T, i_T, OM_T, Th1_T, mu_t);
[rr2_T, vv2_T] = par2carFUN3(a_T, e_T, i_T, OM_T, Th2_T, mu_t);

% Delta V @(initial and final true , pericenter argument)
DeltaV = @(x) norm(vv1_T(x(1), x(2), x(3)) - VV_i(x(1))) + norm(VV_f(x(2)) - vv2_T(x(1), x(2), x(3)));

opts_ga = optimoptions('ga', 'Display', 'off', 'PlotFcn', []);        
lb = [0, 0, 0];           % Lower bounds for [TH_i, TH_f, om_T]
ub = [2*pi, 2*pi, 2*pi];  % Upper bounds [TH_i, TH_f, om_T]
RT = 6378 + 100; % Karman line

cond = @(x) deal([
    e_T(x(1), x(2), x(3)) - 1;    % e_T ≤ 1
   -e_T(x(1), x(2), x(3));        % e_T ≥ 0
   RT - a_T(x(1), x(2), x(3)) * (1 - e_T(x(1), x(2), x(3)))  % Pericenter radius > Karman line
], []);

rng(42); % Setting random seed for stochastic methods

% Minimizing Delta V
[OPT, minV] = ga(DeltaV, 3,[],[],[],[], lb, ub, cond, opts_ga);

% Calculating impulse cost
DeltaV1_MINV = norm(vv1_T(OPT(1), OPT(2), OPT(3)) - VV_i(OPT(1)));
DeltaV2_MINV = norm(VV_f(OPT(2)) - vv2_T(OPT(1), OPT(2), OPT(3)));


% Wait first impulse
DeltaT1 = TOF(a_i, e_i, th_i, OPT(1), mu_t);
% Wait second impulse
DeltaT2 = TOF(a_T(OPT(1),OPT(2),OPT(3)), e_T(OPT(1),OPT(2),OPT(3)), Th1_T(OPT(1),OPT(2),OPT(3)),Th2_T(OPT(1),OPT(2),OPT(3)), mu_t);
% Wait arrival at final position
DeltaT3 = TOF(a_f, e_f, OPT(2), th_f, mu_t);


% Delta t
DeltaT_minV = DeltaT1 + DeltaT2 + DeltaT3;

% Saving results
DeltaVV = [DeltaVV, minV];
DeltaTT = [DeltaTT, DeltaT_minV];


figure('Name','A5 - Direct transfer Delta V MIN (ga)','Units', 'normalized', 'OuterPosition', [0 0 1 1]);
hold on
grid on
box on
axis equal
earth_sphere('km');

% Initial orbit and position
plotOrbit(a_i, e_i, i_i, OM_i, om_i, 0, 2*pi, pi/10000, mu_t);
scatter3(rr_i(1),rr_i(2),rr_i(3),'b*');

% Final orbit and position
plotOrbit(a_f, e_f, i_f, OM_f, om_f, 0, 2*pi, pi/10000, mu_t);
scatter3(rr_f(1),rr_f(2),rr_f(3),'*');

% Transfer orbit
plotOrbit(a_T(OPT(1),OPT(2),OPT(3)), e_T(OPT(1),OPT(2),OPT(3)), i_T(OPT(1),OPT(2)), OM_T(OPT(1),OPT(2)), OPT(3), 0, 2*pi, pi/10000, mu_t)

% Wait first impulse
plotOrbit(a_i, e_i, i_i, OM_i, om_i, th_i, OPT(1), pi/10000, mu_t, 3);
% Wait second impulse
plotOrbit(a_T(OPT(1),OPT(2),OPT(3)), e_T(OPT(1),OPT(2),OPT(3)), i_T(OPT(1),OPT(2)), OM_T(OPT(1),OPT(2)), OPT(3), Th1_T(OPT(1),OPT(2),OPT(3)), Th2_T(OPT(1),OPT(2),OPT(3)), pi/10000, mu_t, 3)
% Wait arrival at final destination
plotOrbit(a_f, e_f, i_f, OM_f, om_f, OPT(2), th_f, pi/10000, mu_t, 3);


leg = legend('', 'Inital orbit', 'Initial position', 'Final orbit','Final position', '','Wait first impulse','Wait second impulse', 'Wait final position','Interpreter','latex');
set(leg, 'FontSize',22)

m_p = m_0 *(exp(minV * 1e3/I_s)-1); % Propellant mass

fprintf("================================================================\n"  + ...
        "Table: A5 with ga \n" + ...
        "Optimal Delta V =  %.6f km/s \n" + ...
        "Respective Delta T = %.4f s = %.4f gg \n" + ...
        "Propellant mass (with Is = 300 s & M_empty = 300 kg): %.2f kg\n" + ...
        "================================================================\n\n", minV, DeltaT_minV, DeltaT_minV/(3600*24), m_p);

% Table A5 - Direct transfer Delta V MIN (ga)
A5_ga(1,:) = [0, a_i, e_i, i_i, OM_i, om_i, th_i, 0]; 
A5_ga(2,:) = [DeltaT1, a_i, e_i, i_i, OM_i, om_i, OPT(1), DeltaV1_MINV];
A5_ga(3,:) = [DeltaT1, a_T(OPT(1),OPT(2),OPT(3)), e_T(OPT(1),OPT(2),OPT(3)), i_T(OPT(1),OPT(2)), OM_T(OPT(1),OPT(2)), OPT(3), Th1_T(OPT(1),OPT(2),OPT(3)), 0];
A5_ga(4,:) = [DeltaT2 + DeltaT1, a_T(OPT(1),OPT(2),OPT(3)), e_T(OPT(1),OPT(2),OPT(3)), i_T(OPT(1),OPT(2)), OM_T(OPT(1),OPT(2)), OPT(3), Th2_T(OPT(1),OPT(2),OPT(3)), DeltaV2_MINV];
A5_ga(5,:) = [DeltaT2 + DeltaT1, a_f, e_f, i_f, OM_f, om_f, OPT(2), 0];
A5_ga(6,:) = [DeltaT_minV, a_f, e_f, i_f, OM_f, om_f, th_f, 0];


TAB5_ga = array2table(A5_ga, 'VariableNames', ...
    {'Time [s]', 'a [km]', 'e [-]', 'i [rad]', 'OM [rad]', 'om [rad]', 'theta [rad]', 'DeltaV [km/s]'});
disp(TAB5_ga)