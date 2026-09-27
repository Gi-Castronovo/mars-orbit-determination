function ds = mars_dynamics_stm(~, s, p, n)
%MARS_DYNAMICS_STM Mars orbiter dynamics (monopole gravity + exponential-
%   atmosphere drag) augmented with the State Transition Matrix.
%
%   Estimated state (n = 6): X = [x, y, vx, vy, GM, rho_0]
%   Estimated state (n = 7): X = [x, y, vx, vy, GM, rho_0, C_D]
%   s = [X; Phi(:)], with Phi the n x n STM, dPhi/dt = A(X) * Phi.
%   Parameters (GM, rho_0, C_D) are constant: their derivatives are zero.

X = s(1:n);
[a, ~] = mars_acceleration(X, p, n);

Phi = reshape(s(n+1:end), n, n);
dPhi = mars_dynamics_jacobian(X, p, n) * Phi;

ds = [X(3); X(4); a; zeros(n-4, 1); dPhi(:)];
end

function [a, rho] = mars_acceleration(X, p, n)
x = X(1);  y = X(2);  vx = X(3);  vy = X(4);  GM = X(5);  rho_0 = X(6);
if n == 7, Cd = X(7); else, Cd = p.C_d; end

r   = sqrt(x^2 + y^2);
v   = sqrt(vx^2 + vy^2);
rho = rho_0*exp(-(r - p.R_m - p.z_0)/p.H);     % exponential atmosphere
k   = 0.5*rho*Cd*p.area/p.m_sc;

a = [-GM/r^3*x - k*v*vx;
     -GM/r^3*y - k*v*vy];
end
