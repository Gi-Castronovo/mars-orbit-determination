function est = batch_least_squares(X0, sigma0, obs, p, n_iter)
%BATCH_LEAST_SQUARES Batch weighted least-squares estimator with a priori.
%   Iterates the linearised normal equations (Tapley, Schutz & Born):
%       x_hat = (H' W H + W_ap) \ (H' W y + W_ap x_bar)
%       P     = (H' W H + W_ap)^-1
%   where H maps deviations of the epoch state to the observations through
%   the STM, W contains the observation weights and W_ap = P_0^-1.
%
%   X0     : a priori state (1 x n)          sigma0 : a priori 1-sigma (1 x n)
%   obs    : struct with fields t, range, rate (column vectors)
%   p      : model parameters                n_iter : number of iterations
%
%   est.X, est.P            : estimate and covariance at t0
%   est.res_range/res_rate  : residuals of every iteration (N x n_iter)
%   est.H                   : observation matrix of the last iteration

n  = numel(X0);
N  = numel(obs.t);
X  = X0(:);
x_bar = zeros(n, 1);
W_ap  = diag(1./sigma0(:).^2);
w     = repmat([1/p.sigma_range^2; 1/p.sigma_rate^2], N, 1);   % diagonal of W (interleaved)

opts  = odeset('RelTol', 5e-14, 'AbsTol', 5e-14);
tspan = [0; obs.t(:)];

est.res_range = zeros(N, n_iter);
est.res_rate  = zeros(N, n_iter);

for it = 1:n_iter
    [~, S] = ode113(@(t, s) mars_dynamics_stm(t, s, p, n), tspan, [X; reshape(eye(n), [], 1)], opts);

    H = zeros(2*N, n);
    y = zeros(2*N, 1);
    for i = 1:N
        s   = S(i+1, :)';                          % row 1 is t0
        Phi = reshape(s(n+1:end), n, n);
        [z, Ht] = range_observables(s(1:n), ground_station_state(obs.t(i), p), n);
        rows = 2*i-1 : 2*i;
        H(rows, :) = Ht*Phi;
        y(rows)    = [obs.range(i); obs.rate(i)] - z;
    end

    Lambda = H'*(w.*H) + W_ap;                    % information matrix
    x_hat  = Lambda \ (H'*(w.*y) + W_ap*x_bar);
    P      = inv(Lambda);

    X     = X + x_hat;                             % update the reference trajectory
    x_bar = x_bar - x_hat;                         % a priori deviation w.r.t. the new reference

    est.res_range(:, it) = y(1:2:end);
    est.res_rate(:, it)  = y(2:2:end);
    fprintf('  iteration %d: RMS range = %.3e m, RMS range-rate = %.3e mm/s, |x_hat| = %.2e\n', it, ...
            1e3*sqrt(mean(y(1:2:end).^2)), 1e6*sqrt(mean(y(2:2:end).^2)), norm(x_hat));
end

est.X = X';
est.P = P;
est.H = H;
end
