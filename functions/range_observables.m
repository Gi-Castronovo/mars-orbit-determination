function [z, Ht] = range_observables(X, Xs, n)
%RANGE_OBSERVABLES Computed range / range-rate and their partial derivatives.
%   X  : spacecraft state (first 4 components = position and velocity)
%   Xs : ground-station state [x; y; vx; vy]
%   z  = [range; range_rate],  Ht = dz/dX (2 x n, zero for the parameters)

dr = X(1:2) - Xs(1:2);  dr = dr(:);
dv = X(3:4) - Xs(3:4);  dv = dv(:);

rho     = norm(dr);
rho_dot = dr'*dv/rho;
z = [rho; rho_dot];

Ht = zeros(2, n);
Ht(1, 1:2) = dr'/rho;
Ht(2, 1:2) = dv'/rho - rho_dot*dr'/rho^2;
Ht(2, 3:4) = dr'/rho;
end
