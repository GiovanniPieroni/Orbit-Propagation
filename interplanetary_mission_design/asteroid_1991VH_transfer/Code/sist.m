function fun = sist(var, alpha, r_H, a_H)

% function fun = sist(var, alfa, r_H, a_H)
% This function contains the system of equation needed to calculate the
% eccentricity e_H, theta_inf and theta_H for the non-coplanar patched
% conics method.
% Input arguments:
% var               [-]   vector of variables         [-]
% alpha             [1x1] angle between r_H and r_inf [rad]
% r_H               [1x1] position on the hyperbola   [km]
% a_H               [1x1] hyperbola semi-major axis   [km]


e_H = var(1);
th_inf = var(2);
th_H = var(3);

fun(1) = cos(th_inf) + 1 / e_H;
fun(2) = th_inf - alpha - th_H;
fun(3) = r_H - (a_H * (1 - e_H^2) )/ (1 + e_H * cos(th_H));
end