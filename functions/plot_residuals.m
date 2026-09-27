function fig = plot_residuals(t, est, iterations, p, name)
%PLOT_RESIDUALS Range and range-rate residuals for the selected iterations.
%   Statistics: mean residual, RMS and normalised chi-square
%   chi2/N = mean((residual/sigma)^2), with sigma the declared accuracy.
%   chi2/N ~ 1 means residuals consistent with the declared measurement noise.

fig = figure('Name', name);
nI  = numel(iterations);
tiledlayout(nI, 2, 'TileSpacing', 'compact');

for j = 1:nI
    k  = iterations(j);
    rr = est.res_range(:, k);
    rv = est.res_rate(:, k);

    nexttile;
    plot(t/3600, rr*1e3, 'k.'); grid on;
    xlabel('time [h]'); ylabel('range residual [m]');
    title(sprintf('Iteration %d - mean %.2e m, RMS %.2e m, \\chi^2/N %.3g', ...
          k, 1e3*mean(rr), 1e3*sqrt(mean(rr.^2)), mean((rr/p.sigma_range).^2)));

    nexttile;
    plot(t/3600, rv*1e6, '.', 'Color', [0.85 0.33 0.10]); grid on;
    xlabel('time [h]'); ylabel('range-rate residual [mm/s]');
    title(sprintf('Iteration %d - mean %.2e mm/s, RMS %.2e mm/s, \\chi^2/N %.3g', ...
          k, 1e6*mean(rv), 1e6*sqrt(mean(rv.^2)), mean((rv/p.sigma_rate).^2)));
end
end
