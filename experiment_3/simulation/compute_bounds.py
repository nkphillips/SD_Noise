import math

# Contrast
contrast_max = 0.9
contrast_min = 0.10
guess_rate = 0.5
lambd = 0.01
target_levels = [0.65, 0.75, 0.85]
neg_log_1mK = [-math.log(1 - (t - guess_rate) / (1 - guess_rate - lambd)) for t in target_levels]

n_grid = 300
alpha_range = [0.01 + i*(0.80-0.01)/(n_grid-1) for i in range(n_grid)]
beta_range = [0.5 + i*(5.0-0.5)/(n_grid-1) for i in range(n_grid)]

weibull_valid_coords = []

for a in alpha_range:
    for b in beta_range:
        valid = True
        for t in range(3):
            c_level = a * (neg_log_1mK[t] ** (1 / b))
            if not (contrast_min <= c_level <= contrast_max):
                valid = False
                break
        
        if valid:
            P_at_cmin = guess_rate + (1 - guess_rate - lambd) * (1 - math.exp(-(contrast_min / a) ** b))
            P_at_cmax = guess_rate + (1 - guess_rate - lambd) * (1 - math.exp(-(contrast_max / a) ** b))
            if (P_at_cmin > 0.55) and (P_at_cmin < 0.90) and (P_at_cmax > 0.85):
                weibull_valid_coords.append((a, b))

alpha_feasible_min = min(x[0] for x in weibull_valid_coords)
alpha_feasible_max = max(x[0] for x in weibull_valid_coords)
beta_feasible_min = min(x[1] for x in weibull_valid_coords)
beta_feasible_max = max(x[1] for x in weibull_valid_coords)

print('=== WEIBULL (Contrast) Feasible Region ===')
print(f'alpha: [{alpha_feasible_min:.4f}, {alpha_feasible_max:.4f}]')
print(f'beta:  [{beta_feasible_min:.4f}, {beta_feasible_max:.4f}]')

# Filter Width
filter_width_max = 80
filter_width_min = 10
prec_min = 1 / filter_width_max
prec_max = 1 / filter_width_min

alpha_fw_range = [0.005 + i*(0.15-0.005)/(n_grid-1) for i in range(n_grid)]
beta_fw_range = [0.5 + i*(5.0-0.5)/(n_grid-1) for i in range(n_grid)]

fw_valid_coords = []

for a in alpha_fw_range:
    for b in beta_fw_range:
        valid = True
        for t in range(3):
            prec_level = a * (neg_log_1mK[t] ** (1 / b))
            fw_level = 1 / prec_level
            if not (filter_width_min <= fw_level <= filter_width_max):
                valid = False
                break
        
        if valid:
            P_at_prec_min = guess_rate + (1 - guess_rate - lambd) * (1 - math.exp(-(prec_min / a) ** b))
            P_at_prec_max = guess_rate + (1 - guess_rate - lambd) * (1 - math.exp(-(prec_max / a) ** b))
            if (P_at_prec_max > 0.75) and (P_at_prec_min > 0.52) and (P_at_prec_min < 0.80):
                fw_valid_coords.append((a, b))

alpha_fw_feasible_min = min(x[0] for x in fw_valid_coords)
alpha_fw_feasible_max = max(x[0] for x in fw_valid_coords)
beta_fw_feasible_min = min(x[1] for x in fw_valid_coords)
beta_fw_feasible_max = max(x[1] for x in fw_valid_coords)

print('\n=== WEIBULL on PRECISION (Filter Width) Feasible Region ===')
print(f'alpha_fw (precision threshold): [{alpha_fw_feasible_min:.4f}, {alpha_fw_feasible_max:.4f}]')
print(f'beta_fw  (slope):               [{beta_fw_feasible_min:.4f}, {beta_fw_feasible_max:.4f}]')

# RECOMMENDED simulate_experiment.m PARAMETERS
rec_alpha_mean = (alpha_feasible_min + alpha_feasible_max) / 2
rec_alpha_std  = (alpha_feasible_max - alpha_feasible_min) / 6
rec_beta_mean  = (beta_feasible_min + beta_feasible_max) / 2
rec_beta_std   = (beta_feasible_max - beta_feasible_min) / 6

rec_alpha_fw_mean = (alpha_fw_feasible_min + alpha_fw_feasible_max) / 2
rec_alpha_fw_std  = (alpha_fw_feasible_max - alpha_fw_feasible_min) / 6
rec_beta_fw_mean  = (beta_fw_feasible_min + beta_fw_feasible_max) / 2
rec_beta_fw_std   = (beta_fw_feasible_max - beta_fw_feasible_min) / 6

print('\n=== RECOMMENDED simulate_experiment.m PARAMETERS ===')
print('\n--- Contrast Weibull ---')
print(f'gt.contrast_alpha = {rec_alpha_mean:.4f} + {rec_alpha_std:.4f} * randn();')
print(f'gt.contrast_alpha = max({alpha_feasible_min:.4f}, min({alpha_feasible_max:.4f}, gt.contrast_alpha));')
print(f'gt.contrast_beta  = {rec_beta_mean:.4f} + {rec_beta_std:.4f} * randn();')
print(f'gt.contrast_beta  = max({beta_feasible_min:.4f}, min({beta_feasible_max:.4f}, gt.contrast_beta));')

print('\n--- Filter Width Weibull (precision = 1/fw) ---')
print(f'gt.filter_alpha = {rec_alpha_fw_mean:.4f} + {rec_alpha_fw_std:.4f} * randn();')
print(f'gt.filter_alpha = max({alpha_fw_feasible_min:.4f}, min({alpha_fw_feasible_max:.4f}, gt.filter_alpha));')
print(f'gt.filter_beta  = {rec_beta_fw_mean:.4f} + {rec_beta_fw_std:.4f} * randn();')
print(f'gt.filter_beta  = max({beta_fw_feasible_min:.4f}, min({beta_fw_feasible_max:.4f}, gt.filter_beta));')

print('\n--- fit_calibration.m bounds ---')
lb_c = [max(0.001, alpha_feasible_min * 0.5), max(0.5, beta_feasible_min * 0.5)]
ub_c = [min(1.0, alpha_feasible_max * 1.5), min(5.0, beta_feasible_max * 1.5)]
lb_fw = [max(0.001, alpha_fw_feasible_min * 0.5), max(0.5, beta_fw_feasible_min * 0.5)]
ub_fw = [min(0.2, alpha_fw_feasible_max * 1.5), min(5.0, beta_fw_feasible_max * 1.5)]

print(f'lb_c = [{lb_c[0]:.4f}, {lb_c[1]:.4f}];  ub_c = [{ub_c[0]:.4f}, {ub_c[1]:.4f}];')
print(f'lb_fw = [{lb_fw[0]:.4f}, {lb_fw[1]:.4f}];  ub_fw = [{ub_fw[0]:.4f}, {ub_fw[1]:.4f}];')