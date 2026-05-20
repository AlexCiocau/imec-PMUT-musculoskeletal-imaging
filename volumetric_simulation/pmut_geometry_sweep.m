% pmut_geometry_sweep.m
%
% Sweep alternative array geometries and tabulate focal-quality metrics.
%
% Constraints assumed:
%   - 64 addressable channels (per-channel transmit delay)
%   - Roughly comparable overall aperture to the current chip
%   - Cells per channel, cell pitch, and operating frequency may vary
%
% Output:
%   - A table comparing peak depth, astigmatic separation, FWHM, focal
%     gain, and dB gain vs the baseline.
%   - On-axis intensity curves for visual comparison.
%   - Predicted physical focal pressure for each geometry, using the user's
%     measured baseline (e.g. 4 kPa/V at 50 mm) as a calibration anchor.
%
% Requires:  pmut_metrics.m, pressure_budget.m, alpha_cal.mat
%            (alpha_cal.mat is produced by calibrate_from_planewave.m)

clear; clc; close all;

%% ---------- 1. Configuration -------------------------------------------
% Set to true to compute the absolute kPa column using your measured
% plane-wave pressure curve as an anchor. Set to false for pure-simulation
% comparisons (recommended - all dB ratios are unaffected by calibration).
USE_MEASURED_CALIBRATION = false;

V_drive_max = 30;       % maximum drive voltage you can apply (V)

if USE_MEASURED_CALIBRATION
    load('alpha_cal.mat', 'alpha_cal', 'rmse_pct', ...
                          'z_sim', 'I_sim_z', 'z_meas', 'p_meas');
    fprintf('Loaded calibration: alpha_cal = %.3e (validation RMS %.1f%%)\n', ...
            alpha_cal, rmse_pct);
end

V_drive_max = 30;       % maximum drive voltage you can apply (V)

%% ---------- 2. Common simulation parameters ---------------------------
c              = 1480;
z_focus_target = 30e-3;     % focal depth used for ALL configs
mode_focus     = 'plane';   % 'focus' or 'plane'

% Imaging volume - keep modest for sweep speed; refine after the headline
% configs are chosen.
Nx = 41; Ny = 41; Nz = 120;
x_grid = single(linspace(-6e-3,  6e-3, Nx));
y_grid = single(linspace(-6e-3,  6e-3, Ny));
z_grid = single(linspace( 0.5e-3, 100e-3, Nz));

%% ---------- 3. Configurations to sweep --------------------------------
% Each config keeps 64 channels and roughly the same lateral footprint.
% Add/remove rows as needed.
configs = struct('name', {}, 'fc', {}, 'N_ch', {}, 'ch_pitch', {}, ...
                 'N_cell_x', {}, 'N_cell_y', {}, ...
                 'cell_pitch_x', {}, 'cell_pitch_y', {});

configs(end+1) = struct( ...
    'name',         'Baseline 1x80, 75um, fc=10.5MHz', ...
    'fc',           10.5e6, ...
    'N_ch',         64, ...
    'ch_pitch',     75e-6, ...
    'N_cell_x',     1, ...
    'N_cell_y',     80, ...
    'cell_pitch_x', 75e-6, ...
    'cell_pitch_y', 75e-6);

% configs(end+1) = struct( ...
%     'name',         'Square 1x64, 75um, fc=10.5MHz', ...
%     'fc',           10.5e6, ...
%     'N_ch',         64, ...
%     'ch_pitch',     75e-6, ...
%     'N_cell_x',     1, ...
%     'N_cell_y',     64, ...
%     'cell_pitch_x', 75e-6, ...
%     'cell_pitch_y', 75e-6);
% 
% configs(end+1) = struct( ...
%     'name',         'Square optimized, 75um, fc=10.5MHz', ...
%     'fc',           10.5e6, ...
%     'N_ch',         64, ...
%     'ch_pitch',     77.241e-6, ...
%     'N_cell_x',     1, ...
%     'N_cell_y',     64, ...
%     'cell_pitch_x', 65.3e-6, ...
%     'cell_pitch_y', 65.3e-6);


% configs(end+1) = struct( ...
%     'name',         'Taller 1x128, 75um, fc=10.5MHz', ...
%     'fc',           10.5e6, ...
%     'N_ch',         64, ...
%     'ch_pitch',     75e-6, ...
%     'N_cell_x',     1, ...
%     'N_cell_y',     128, ...
%     'cell_pitch_x', 75e-6, ...
%     'cell_pitch_y', 75e-6);
% 
% configs(end+1) = struct( ...
%     'name',         'Baseline geom, fc=8MHz', ...
%     'fc',           8e6, ...
%     'N_ch',         64, ...
%     'ch_pitch',     75e-6, ...
%     'N_cell_x',     1, ...
%     'N_cell_y',     80, ...
%     'cell_pitch_x', 75e-6, ...
%     'cell_pitch_y', 75e-6);
% 
% configs(end+1) = struct( ...
%     'name',         'Wide y-pitch 100um, 1x80', ...
%     'fc',           10.5e6, ...
%     'N_ch',         64, ...
%     'ch_pitch',     75e-6, ...
%     'N_cell_x',     1, ...
%     'N_cell_y',     80, ...
%     'cell_pitch_x', 75e-6, ...
%     'cell_pitch_y', 100e-6);

%% ---------- 4. Run the sweep ------------------------------------------
results = cell(numel(configs), 1);
for ic = 1:numel(configs)
    fprintf('\n=== Config %d: %s ===\n', ic, configs(ic).name);
    [I3, N_cells, lambda] = run_pmut_sim(configs(ic), c, z_focus_target, ...
                                          mode_focus, x_grid, y_grid, z_grid);
    m = pmut_metrics(I3, x_grid, y_grid, z_grid, ...
                     'N_cells', N_cells, 'z_target', z_focus_target);
    m.config  = configs(ic);
    m.lambda  = lambda;
    m.N_cells = N_cells;
    results{ic} = m;
end

%% ---------- 5. Pre-compute the measured/baseline reference -----------
% Only needed if USE_MEASURED_CALIBRATION is true. Otherwise, all
% comparisons are made in dB relative to the baseline configuration.
V_drive_meas = 1;

[~, ix0] = min(abs(double(x_grid)));
[~, iy0] = min(abs(double(y_grid)));

if USE_MEASURED_CALIBRATION
    p_meas_on_zgrid    = interp1(z_meas, p_meas,   double(z_grid), 'linear', NaN);
    I_PW_base_on_zgrid = interp1(z_sim,  I_sim_z,  double(z_grid), 'linear', NaN);
end

%% ---------- 6. Build comparison table --------------------------------
T = table;
for ic = 1:numel(results)
    r = results{ic};

    % Pure-sim focal-gain ratio vs baseline AT THE TARGET DEPTH
    % (most useful for comparing configurations at a fixed imaging depth)
    gain_at_target_dB = 10*log10(r.focal_gain_at_target / ...
                                  results{1}.focal_gain_at_target);
    % Also keep the peak-based number for reference
    gain_at_peak_dB   = 10*log10(r.focal_gain_eff / results{1}.focal_gain_eff);

    T.Name(ic,1)               = string(r.config.name);
    T.fc_MHz(ic,1)             = r.config.fc / 1e6;
    T.N_cells(ic,1)            = r.N_cells;
    T.D_x_mm(ic,1)             = ((r.config.N_ch-1)*r.config.ch_pitch + ...
                                   r.config.N_cell_x*r.config.cell_pitch_x) * 1e3;
    T.D_y_mm(ic,1)             = (r.config.N_cell_y-1) * r.config.cell_pitch_y * 1e3;
    T.zNFx_mm(ic,1)            = (T.D_x_mm(ic)*1e-3)^2 / (4*r.lambda) * 1e3;
    T.zNFy_mm(ic,1)            = (T.D_y_mm(ic)*1e-3)^2 / (4*r.lambda) * 1e3;
    T.z_peak_mm(ic,1)          = r.z_peak * 1e3;
    T.astig_sep_mm(ic,1)       = r.astig_separation * 1e3;
    T.FWHM_x_mm(ic,1)          = r.FWHM_x * 1e3;
    T.FWHM_y_mm(ic,1)          = r.FWHM_y * 1e3;
    T.DOF_mm(ic,1)             = r.DOF_6dB * 1e3;
    T.GainAtTarget_dB(ic,1)    = gain_at_target_dB;
    T.GainAtPeak_dB(ic,1)      = gain_at_peak_dB;

    % Optional kPa prediction (anchored to measurement; carries calibration uncertainty)
    if USE_MEASURED_CALIBRATION
        I_focused_on_axis = double(r.I_peak) * 10.^(r.on_axis_dB / 10);
        p_focused_profile = p_meas_on_zgrid(:) .* ...
                            sqrt(I_focused_on_axis(:) ./ I_PW_base_on_zgrid(:)) ...
                            .* (V_drive_max / V_drive_meas);
        [~, iz_t]   = min(abs(double(z_grid) - z_focus_target));
        T.p_focal_kPa(ic,1) = p_focused_profile(iz_t) / 1e3;
        results{ic}.p_focused_profile = p_focused_profile;
    end
end
fprintf('\n=== Geometry comparison ===\n');
disp(T);

%% ---------- 7. Pressure-budget overlay (pure dB) ---------------------
% In dB above the baseline configuration's focal pressure at z_focus_target.
% No absolute pressure values needed - this is the depth-extension argument
% for the design team in its cleanest form.
params.alpha_dB_per_cm_per_MHz = 0.5;       % CIRS 040 GSE
params.f_MHz                   = 10.5;
params.z_baseline_m            = z_focus_target;
params.p_focus_baseline_Pa     = 1;          % dummy - only ratios are used

z_query = linspace(z_focus_target, 80e-3, 200);
budget  = pressure_budget(z_query, params);

figure('Color','w','Position',[80 80 760 480]);
plot(z_query*1e3, budget.deficit_dB, 'k-', 'LineWidth', 1.8, ...
     'DisplayName','dB needed (CIRS attenuation + 1/z spreading)'); hold on;
yline(0, 'r-', 'LineWidth', 1.4, ...
      'Label','Baseline performance (0 dB)', 'DisplayName','Baseline');

% Overlay each candidate's available focal-gain improvement
colors = lines(numel(results));
for ic = 2:numel(results)
    yline(T.GainAtTarget_dB(ic), '--', 'Color', colors(ic,:), 'LineWidth',1.3, ...
          'Label', sprintf('%s: %+.1f dB', T.Name(ic), T.GainAtTarget_dB(ic)), ...
          'Interpreter','none', 'LabelHorizontalAlignment','left');
end
xlabel('Depth z (mm)');
ylabel('dB above baseline focal pressure');
title(sprintf('Depth-extension budget (baseline focuses at z = %.0f mm in CIRS 040 GSE)', ...
              z_focus_target*1e3));
grid on; legend('Location','southeast');

%% ---------- 8. On-axis comparison plot --------------------------------
figure('Color','w','Position',[80 580 760 420]);
hold on;
for ic = 1:numel(results)
    plot(double(results{ic}.z_grid)*1e3, results{ic}.on_axis_dB, ...
         'LineWidth', 1.4, 'DisplayName', results{ic}.config.name);
end
xlabel('Depth z (mm)');
ylabel('On-axis intensity (dB rel each config''s peak)');
title(sprintf('On-axis intensity, focused at z = %.0f mm', z_focus_target*1e3));
legend('Location','best','Interpreter','none');
grid on; ylim([-25 1]);
xline(z_focus_target*1e3, 'k--', 'LineWidth', 1, 'Label','target focus');
hold off;


%% ====================== local helper function =========================
function [I3, N_cells, lambda] = run_pmut_sim(cfg, c, z_focus, mode, ...
                                               x_grid, y_grid, z_grid)
% Single-shot Huygens-Fresnel coherent sum for one geometry.
% Uses 1/r monopole weighting (corrected from 1/sqrt(r)).

lambda = c / cfg.fc;
k      = 2*pi/lambda;

% --- cell positions ---
[ix, iy, ich]    = ndgrid(1:cfg.N_cell_x, 1:cfg.N_cell_y, 1:cfg.N_ch);
ich              = ich(:);
ch_centre_x      = (ich - (cfg.N_ch+1)/2) * cfg.ch_pitch;
dx_in_ch         = (ix(:) - (cfg.N_cell_x+1)/2) * cfg.cell_pitch_x;
dy_in_ch         = (iy(:) - (cfg.N_cell_y+1)/2) * cfg.cell_pitch_y;
x_cell           = ch_centre_x + dx_in_ch;
y_cell           = dy_in_ch;
N_cells          = numel(x_cell);

% --- per-channel transmit delays ---
ch_x = ((1:cfg.N_ch) - (cfg.N_ch+1)/2) * cfg.ch_pitch;
switch mode
    case 'plane'
        ch_tau = zeros(1, cfg.N_ch);
    case 'focus'
        d_ch    = sqrt(ch_x.^2 + z_focus^2);
        ch_tau  = (max(d_ch) - d_ch) / c;
    otherwise
        error('mode must be ''plane'' or ''focus''.');
end
cell_tau = ch_tau(ich);

% --- coherent sum ---
[X3, Y3, Z3] = ndgrid(x_grid, y_grid, z_grid);
xc = single(x_cell); yc = single(y_cell); tc = single(cell_tau);
ks = single(k);      ws = single(2*pi*cfg.fc);

P3 = complex(zeros(size(X3), 'single'));
fprintf('  %d cells over %dx%dx%d voxels ...\n', N_cells, ...
        size(X3,1), size(X3,2), size(X3,3));
tic;
for i = 1:N_cells
    if mod(i, 1024) == 0
        fprintf('   cell %5d / %5d  (%.1fs)\n', i, N_cells, toc);
    end
    r     = sqrt((X3 - xc(i)).^2 + (Y3 - yc(i)).^2 + Z3.^2);
    phase = -(ks*r + ws*tc(i));
    P3    = P3 + exp(1i*phase) ./ r;          % 1/r monopole weighting
end
fprintf('  done (%.1fs).\n', toc);
I3 = abs(P3).^2;
end