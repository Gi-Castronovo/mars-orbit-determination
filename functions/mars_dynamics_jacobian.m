function A = mars_dynamics_jacobian(X, p, n)
%MARS_DYNAMICS_JACOBIAN Analytical Jacobian A = dF/dX of the orbiter dynamics.
%   n = 6 or 7 (see MARS_DYNAMICS_STM). Rows 5..n are zero (constant parameters).

x = X(1);  y = X(2);  vx = X(3);  vy = X(4);  GM = X(5);  rho_0 = X(6);
if n == 7, Cd = X(7); else, Cd = p.C_d; end

r  = sqrt(x^2 + y^2);
v  = sqrt(vx^2 + vy^2);
e  = exp(-(r - p.R_m - p.z_0)/p.H);            % rho = rho_0 * e
c  = 0.5*p.area/p.m_sc;                        % drag constant: a_drag = -c*Cd*rho*v*[vx; vy]
k  = c*Cd*rho_0*e;

A = zeros(n);
A(1,3) = 1;
A(2,4) = 1;

% Gravity (monopole) + drag through d(rho)/dr = -rho/H
A(3,1) = -GM/r^3 + 3*GM*x^2/r^5 + k*v*vx*x/(p.H*r);
A(3,2) =           3*GM*x*y/r^5 + k*v*vx*y/(p.H*r);
A(4,1) =           3*GM*x*y/r^5 + k*v*vy*x/(p.H*r);
A(4,2) = -GM/r^3 + 3*GM*y^2/r^5 + k*v*vy*y/(p.H*r);

% Drag: derivatives with respect to the velocity
A(3,3) = -k*(2*vx^2 + vy^2)/v;
A(3,4) = -k*vx*vy/v;
A(4,3) = -k*vx*vy/v;
A(4,4) = -k*(vx^2 + 2*vy^2)/v;

% Parameters
A(3,5) = -x/r^3;
A(4,5) = -y/r^3;
A(3,6) = -c*Cd*e*v*vx;
A(4,6) = -c*Cd*e*v*vy;
if n == 7
    A(3,7) = -c*rho_0*e*v*vx;
    A(4,7) = -c*rho_0*e*v*vy;
end
end
