function Xs = ground_station_state(t, p)
%GROUND_STATION_STATE Position and velocity of the ground station at time t.
%   Mars-centred frame, 2D: the Earth centre is fixed at (d_earth, 0) and
%   the station rotates with the Earth (sidereal day), starting at the
%   angle p.station_angle0. Returns Xs = [x; y; vx; vy] in km and km/s.

w     = 2*pi/p.t_sd;
theta = p.station_angle0 + w*t;
Xs = [p.d_earth + p.R_e*cos(theta);
      p.R_e*sin(theta);
     -p.R_e*w*sin(theta);
      p.R_e*w*cos(theta)];
end
