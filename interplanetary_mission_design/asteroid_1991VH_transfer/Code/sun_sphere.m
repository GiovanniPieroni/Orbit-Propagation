function sun_sphere
% This function plot a yellow sphere that should represent the Sun

R = 696340; % Sun radius

[X, Y, Z] = sphere(50);
X = R * X;
Y = R * Y;
Z = R * Z;

surf(X, Y, Z, 'FaceColor', 'yellow', 'EdgeColor', 'none') 
alpha(1) 
axis equal
xlabel('X'), ylabel('Y'), zlabel('Z')

camlight 
lighting gouraud
end
