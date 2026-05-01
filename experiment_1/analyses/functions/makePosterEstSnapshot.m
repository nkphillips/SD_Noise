function snap = makePosterEstSnapshot(rb, sd, rb_ci, sd_ci, sd_ci_cluster, sd_boot_cluster, perf_ci, delta_theta_centers)
% makePosterEstSnapshot
%
% Builds an in-memory estimates snapshot in the same shape that loadEstimates
% returns for a single n_back. Used to drive poster plotting functions
% inline (without reloading saved .mat files).

    snap = struct();
    snap.rb              = rb;
    snap.sd              = sd;
    snap.rb_ci           = rb_ci;
    snap.sd_ci           = sd_ci;
    snap.sd_ci_cluster   = sd_ci_cluster;
    snap.sd_boot_cluster = sd_boot_cluster;
    snap.perf_ci         = perf_ci;
    snap.delta_theta_centers = delta_theta_centers;
    snap.meta            = struct('delta_theta_centers', delta_theta_centers);
end
