function plotOrbit(a,e,i,OM,om,th0,thf,dth,mu, L)

% 3D orbit plot
%
% plotOrbit(a,e,i,OM,om,th0,thf,dth,mu)
%
%
% Input arguments:
% a        [1x1]  semi-major axis            [km]
% e        [1x1]  eccentricity               [-]
% i        [1x1]  inclination                [rad]
% OM       [1x1]  RAAN                       [rad]
% om       [1x1]  argument of the periapsis  [rad]
% th0      [1x1]  initial true anomaly       [rad]
% thf      [1x1]  final true anomaly         [rad]
% dth      [1x1]  true anomaly step size     [rad]
% mu       [1x1]  gravitational parameter    [km^3/s^2]
% L        [1x1]  (optional) scaling factor 
%                 for hyperbola cutoff or  
%                 the thickness of the conic [-]


% Set default value for L if not provided
if nargin<10 && e>1
    L = 10;
elseif nargin<10 && e<1
    L = 1;
end
     


if e > 1  % hyperbolic orbit
    th_max = acos(-1/e);  % valid values
    theta = -th_max:dth:th_max;
    r_max = L * abs(a * (1 - e));
    L = 1;
elseif thf < th0
    thf = 2 * pi + thf;
    theta = th0:dth:thf;
    r_max = inf;
else
    theta = th0:dth:thf;
    r_max = inf;
end


r = NaN(3, length(theta));

for j = 1:length(theta)
    try
        [r_ij, ~] = par2car(a, e, i, OM, om, theta(j), mu);

        if isreal(r_ij) && all(~isnan(r_ij)) && norm(r_ij) < r_max  % <--- Limit
            r(:, j) = r_ij;
        end

    catch
        
    end
end

hold on
grid on


plot3(r(1,:) , r(2,:), r(3,:), 'LineWidth', L);


xlabel('X','FontSize',15);
ylabel('Y','FontSize',15);
zlabel('Z','FontSize',15);

end



