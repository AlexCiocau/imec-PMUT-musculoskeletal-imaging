%% Multi-configuration depth sweep
%
% Define any number of layout configurations (each with its own channel
% pitch if desired) and sweep the channel focus across a list of depths
% for each. Produces:
%   - One depth-sweep overlay figure per configuration (sweep curves
%     overlaid, gradient-coloured by intended z_focus)
%   - One CAPABILITY plot with all configs overlaid: intensity AT the
%     intended z_focus vs the intended z_focus. The single most
%     diagnostic comparison plot.
%   - One SNAPSHOT plot showing all configs' on-axis profiles when
%     channel focus = z_focus_demo. Quick visual comparison of profile
%     shape at a single steered depth.
%
% Each config can specify its own ch_pitch -- so e.g. baseline at 75um
% and baseline at 150um can be compared directly.
%
% Common normalisation: every dB value is referenced to the GLOBAL max
% across all configs and all sweep depths -- so the curves' relative
% heights reflect real intensity differences.

clear; clc; close all;

%% ---------- Acoustic + array params (common to all configs) -------
c        = 1480;
fc       = 10.5e6;
lambda   = c / fc;
k        = 2*pi / lambda;
N_ch     = 64;          % fixed channel count

%% ---------- Define configurations (USER EDITS THIS BLOCK) ---------
% Each config struct sets cfg.type plus type-specific parameters and
% (optionally) cfg.ch_pitch. Defaults: ch_pitch=75um, N_cell_y=80,
% p_min_y=75um, cell_pitch_x=75um, N_cell_x=1.

% configs(1).name         = 'Pure chirp ratio 4 @ 150um channel pitch';
% configs(1).type         = 'chirp';
% configs(1).chirp_ratio  = 4;
% configs(1).ch_pitch     = 150e-6;

configs(1).name         = 'Equal pitch baseline @ 75um channel pitch';
configs(1).type         = 'equal';
configs(1).cell_pitch_y = 75e-6;
configs(1).ch_pitch     = 75e-6;

configs(2).name         = 'Square @ 75um channel pitch';
configs(2).type         = 'equal';
configs(2).cell_pitch_y = 75.4e-6;
configs(2).ch_pitch     = 89.2e-6;
% 
% configs(2).name         = 'Equal pitch baseline @ 150um channel pitch';
% configs(2).type         = 'equal';
% configs(2).cell_pitch_y = 75e-6;
% configs(2).ch_pitch     = 150e-6;
% 
% configs(3).name         = 'Grouped Fresnel @ 35mm, 150um channel pitch';
% configs(3).type         = 'grouped';
% configs(3).z_focus_geom = 40e-3;
% configs(3).dr_frac      = 1/8;
% configs(3).ch_pitch     = 75e-6;

configs(3).name         = 'Fresnel @ 35mm, 75um channel pitch';
configs(3).type         = 'parabolic';
configs(3).z_focus_geom = 40e-3;
configs(3).ch_pitch     = 75e-6;

configs(4).name         = 'Grouped Fresnel @ 35mm, 150um channel pitch';
configs(4).type         = 'grouped';
configs(4).z_focus_geom = 40e-3;
configs(4).dr_frac      = 1/6;
configs(4).ch_pitch     = 75e-6;

configs(5).name         = 'Pure chirp ratio 4 @ 150um channel pitch';
configs(5).type         = 'chirp';
configs(5).chirp_ratio  = 4;
configs(5).ch_pitch     = 150e-6;

% configs(6).name         = 'Pure chirp ratio 4 @ 150um channel pitch';
% configs(6).type         = 'chirp';
% configs(6).chirp_ratio  = 4;
% configs(6).ch_pitch     = 150e-6;

% configs(6).name         = 'Aggressive Chirp (Ratio 5, 60 cells)';
% configs(6).type         = 'chirp';
% configs(6).p_min_y      = 75e-6;   % Starts at 75um in the center
% configs(6).chirp_ratio  = 5;       % Outer cells grow to 375um
% configs(6).N_cell_y     = 60;      % Yields an aperture D_y around 6.5 mm
% configs(6).ch_pitch     = 75e-6;   % Keeps your azimuth steering tight

% configs(1).name         = 'Equal pitch baseline @ 75um channel pitch';
% configs(1).type         = 'equal';
% configs(1).cell_pitch_y = 75e-6;
% configs(1).ch_pitch     = 75e-6;
% 
% configs(2).name         = 'Grouped Fresnel @ 35mm, 75um channel pitch';
% configs(2).type         = 'grouped';
% configs(2).z_focus_geom = 35e-3;
% configs(2).dr_frac      = 1/8;
% configs(2).ch_pitch     = 75e-6;
% 
% configs(3).name         = 'Pure chirp ratio 4 @ 150um channel pitch';
% configs(3).type         = 'chirp';
% configs(3).chirp_ratio  = 1;
% configs(3).ch_pitch     = 150e-6;

% Fill in defaults
for i = 1:numel(configs)
    if ~isfield(configs(i), 'ch_pitch')     || isempty(configs(i).ch_pitch),     configs(i).ch_pitch     = 75e-6; end
    if ~isfield(configs(i), 'N_cell_x')     || isempty(configs(i).N_cell_x),     configs(i).N_cell_x     = 1;     end
    if ~isfield(configs(i), 'N_cell_y')     || isempty(configs(i).N_cell_y),     configs(i).N_cell_y     = 80;    end
    if ~isfield(configs(i), 'cell_pitch_x') || isempty(configs(i).cell_pitch_x), configs(i).cell_pitch_x = 75e-6; end
    if ~isfield(configs(i), 'p_min_y')      || isempty(configs(i).p_min_y),      configs(i).p_min_y      = 75e-6; end
end


%% ---------- Channel focus sweep -----------------------------------
z_focus_sweep = linspace(10e-3, 40e-3, 11);   % 10, 15, ..., 60 mm
z_focus_demo  = 35e-3;                          % for the snapshot plot
N_sweep       = numel(z_focus_sweep);

%% ---------- Common 1D depth grid ----------------------------------
z_grid = single(linspace(1e-3, 100e-3, 400));
Nz     = numel(z_grid);

%% ---------- Process every configuration ---------------------------
N_cfg = numel(configs);
results(N_cfg) = struct('cfg',[],'profiles',[],'D_y',0,'y_in_ch',[]);

fprintf('Processing %d configurations, %d sweep depths each.\n', N_cfg, N_sweep);
total_tic = tic;
for c_idx = 1:N_cfg
    cfg = configs(c_idx);
    fprintf('  %d/%d: %s\n', c_idx, N_cfg, cfg.name);

    [y_in_ch, ~]          = local_y_in_channel(cfg, lambda);
    [x_cell, y_cell, ich] = local_build_cells(y_in_ch, ...
                                cfg.N_cell_x, cfg.cell_pitch_x, ...
                                N_ch, cfg.ch_pitch);
    ch_x_centres = ((1:N_ch) - (N_ch+1)/2) * cfg.ch_pitch;

    profiles = zeros(N_sweep, Nz, 'single');
    for j = 1:N_sweep
        z_f      = z_focus_sweep(j);
        d_ch     = sqrt(ch_x_centres.^2 + z_f^2);
        ch_tau   = (max(d_ch) - d_ch) / c;
        cell_tau = ch_tau(ich);
        profiles(j, :) = local_on_axis_profile(x_cell, y_cell, cell_tau, ...
                                                z_grid, k, fc);
    end

    results(c_idx).cfg      = cfg;
    results(c_idx).profiles = profiles;
    results(c_idx).D_y      = max(y_cell) - min(y_cell);
    results(c_idx).y_in_ch  = y_in_ch;
end
fprintf('Total time: %.1f s.\n', toc(total_tic));

%% ---------- Global max for common normalisation -------------------
maxI_global = 0;
for i = 1:N_cfg, maxI_global = max(maxI_global, max(results(i).profiles(:))); end

%% ---------- Cell layouts (one figure with rows per config) --------
clr_cfg = lines(N_cfg);
figure('Color','w','Position',[60 60 1300 60+50*N_cfg]);
hold on; grid on;
for i = 1:N_cfg
    plot(results(i).y_in_ch*1e3, i*ones(size(results(i).y_in_ch)), 'o', ...
         'MarkerSize',5,'MarkerFaceColor', clr_cfg(i,:),'Color', clr_cfg(i,:));
end
xline(0,'k:'); yticks(1:N_cfg);
yticklabels(arrayfun(@(c) c.cfg.name, results, 'UniformOutput', false));
ylim([0.5 N_cfg+0.5]);
xlabel('Cell y position within one channel (mm)');
title('Cell layouts (one row per configuration)','Interpreter','none');
set(gca,'TickLabelInterpreter','none');

%% ---------- Per-config depth-sweep overlay (one figure per config) -
clr_sweep = local_make_colormap(N_sweep);
for c_idx = 1:N_cfg
    cfg      = results(c_idx).cfg;
    profiles = results(c_idx).profiles;

    figure('Color','w','Position',[80 + (c_idx-1)*30, 100 + (c_idx-1)*20, 1150, 480]);
    hold on; grid on;
    for j = 1:N_sweep
        plot(double(z_grid)*1e3, ...
             10*log10(double(profiles(j,:))/maxI_global + eps), ...
             'LineWidth', 1.5, 'Color', clr_sweep(j,:), ...
             'DisplayName', sprintf('z_{focus} = %.0f mm', z_focus_sweep(j)*1e3));
        [~, iz] = min(abs(z_grid - z_focus_sweep(j)));
        plot(z_focus_sweep(j)*1e3, ...
             10*log10(double(profiles(j, iz))/maxI_global + eps), ...
             'v', 'Color', clr_sweep(j,:), 'MarkerFaceColor', clr_sweep(j,:), ...
             'MarkerSize', 7, 'HandleVisibility','off');
    end
    if isfield(cfg, 'z_focus_geom') && ~isempty(cfg.z_focus_geom)
        xline(cfg.z_focus_geom*1e3, 'k--', 'LineWidth', 1.2, ...
              'DisplayName', sprintf('z_{geom} = %.0f mm', cfg.z_focus_geom*1e3));
    end
    xlabel('Depth z (mm)');
    ylabel('On-axis intensity (dB, common scale across all configs)');
    ylim([-30 1]);
    title({sprintf('Depth sweep: %s', cfg.name), ...
           sprintf('ch\\_pitch = %.0f um, D_y = %.1f mm', ...
                   cfg.ch_pitch*1e6, results(c_idx).D_y*1e3)}, ...
          'Interpreter','tex');
    legend('Location','eastoutside');
end

%% ---------- KEY PLOT: capability curves overlaid for ALL configs --
figure('Color','w','Position',[80 700 1300 540]);
hold on; grid on;
% Vertical dashed lines at any design depth used in any config
geom_depths = [];
for c_idx = 1:N_cfg
    cfg = results(c_idx).cfg;
    if isfield(cfg,'z_focus_geom') && ~isempty(cfg.z_focus_geom)
        geom_depths(end+1) = cfg.z_focus_geom; %#ok<SAGROW>
    end
end
geom_depths = unique(geom_depths);
for zg = geom_depths
    xline(zg*1e3, ':', 'Color',[0.45 0.45 0.45], 'LineWidth', 1.2, ...
          'DisplayName', sprintf('z_{geom} = %.0f mm', zg*1e3));
end
for c_idx = 1:N_cfg
    profiles = results(c_idx).profiles;
    cap = zeros(1, N_sweep);
    for j = 1:N_sweep
        [~, iz] = min(abs(z_grid - z_focus_sweep(j)));
        cap(j)  = double(profiles(j, iz));
    end
    plot(z_focus_sweep*1e3, ...
         10*log10(cap/maxI_global + eps), ...
         'o-', 'LineWidth', 2.0, 'Color', clr_cfg(c_idx,:), ...
         'MarkerSize', 7, 'MarkerFaceColor', clr_cfg(c_idx,:), ...
         'DisplayName', results(c_idx).cfg.name);
end
xlabel('Intended z_{focus} (mm)');
ylabel('Intensity AT the intended z_{focus} (dB, common scale)');
ylim([-25 5]);
title({'Capability comparison across configs', ...
       'Higher = better focus quality at that intended depth'});
legend('Location','eastoutside','Interpreter','none');

%% ---------- SNAPSHOT: all configs at one channel-focus depth ------
[~, j_demo] = min(abs(z_focus_sweep - z_focus_demo));
figure('Color','w','Position',[80 1300 1200 480]);
hold on; grid on;
xline(z_focus_sweep(j_demo)*1e3, 'k--', 'LineWidth', 1.2, ...
      'DisplayName','channel focus');
for c_idx = 1:N_cfg
    profile = results(c_idx).profiles(j_demo, :);
    plot(double(z_grid)*1e3, ...
         10*log10(double(profile)/maxI_global + eps), ...
         'LineWidth', 1.8, 'Color', clr_cfg(c_idx,:), ...
         'DisplayName', results(c_idx).cfg.name);
end
xlabel('Depth z (mm)');
ylabel('On-axis intensity (dB, common scale)');
ylim([-30 1]);
title({sprintf('Snapshot: all configs with channel focus at z = %.0f mm', ...
               z_focus_sweep(j_demo)*1e3), ...
       'Direct shape comparison at one steered depth'});
legend('Location','eastoutside','Interpreter','none');

%% ---------- Summary table -----------------------------------------
fprintf('\n%s\n', repmat('=', 1, 130));
fprintf('SUMMARY (intensity in dB, common scale, global max = 0 dB)\n');
fprintf('%-50s | %-9s | %-8s | %-9s | %-9s | %-9s | %-13s\n', ...
        'Config', 'ch_pitch', 'D_y(mm)', 'Best dB', 'Worst dB', 'Range dB', 'f_x# @ 35mm');
fprintf('%s\n', repmat('-', 1, 130));
for c_idx = 1:N_cfg
    cfg      = results(c_idx).cfg;
    profiles = results(c_idx).profiles;
    cap = zeros(1, N_sweep);
    for j = 1:N_sweep
        [~, iz] = min(abs(z_grid - z_focus_sweep(j)));
        cap(j)  = double(profiles(j, iz));
    end
    best_dB  = 10*log10(max(cap)/maxI_global + eps);
    worst_dB = 10*log10(max(min(cap), eps)/maxI_global + eps);
    rng_dB   = best_dB - worst_dB;
    fnum_x   = 35e-3 / ((N_ch-1)*cfg.ch_pitch);
    fprintf('%-50s | %6.0f um | %7.2f | %+9.1f | %+9.1f | %9.1f | %12.2f\n', ...
            cfg.name, cfg.ch_pitch*1e6, results(c_idx).D_y*1e3, ...
            best_dB, worst_dB, rng_dB, fnum_x);
end
fprintf('\nLegend: Best/Worst dB are the highest/lowest "I at intended z_focus" across the sweep.\n');
fprintf('        Range = Best - Worst (smaller = more uniform behavior across depths).\n');


%% ====================== HELPER FUNCTIONS =============================
function [y_in_ch, y_label] = local_y_in_channel(cfg, lambda)
    y_label = [];
    switch cfg.type
        case 'equal'
            y_in_ch = ((1:cfg.N_cell_y) - (cfg.N_cell_y+1)/2) * cfg.cell_pitch_y;
        case 'parabolic'
            y_in_ch = local_compute_parabolic_y_strict(cfg.N_cell_y, ...
                cfg.z_focus_geom, lambda, cfg.p_min_y);
        case 'grouped'
            y_in_ch = local_compute_grouped_y(cfg.N_cell_y, cfg.p_min_y, ...
                cfg.z_focus_geom, lambda, cfg.dr_frac);
        case 'chirp'
            y_in_ch = local_compute_chirped_pitch_y(cfg.N_cell_y, ...
                cfg.p_min_y, cfg.chirp_ratio);
        case 'hybrid'
            y_in_ch = local_compute_hybrid_y(cfg.N_cell_y, cfg.p_min_y, ...
                cfg.z_focus_geom, lambda, cfg.N_fresnel, cfg.chirp_ratio);
        case 'multifocal_hybrid'
            [y_in_ch, y_label] = local_compute_multifocal_hybrid_y( ...
                cfg.N_cell_y, cfg.p_min_y, lambda, ...
                cfg.z_focus_design_list, cfg.n_F_per_depth, cfg.chirp_ratio);
        otherwise
            error('Unknown cfg.type: %s', cfg.type);
    end
end

function y_pos = local_compute_chirped_pitch_y(N_cells, p_min, chirp_ratio)
    if N_cells <= 0, y_pos = []; return; end
    if N_cells == 1, y_pos = 0; return; end
    is_odd = mod(N_cells, 2) == 1;
    n_side = floor(N_cells / 2);
    if n_side == 1, pitches = p_min;
    else,           pitches = p_min * linspace(1, chirp_ratio, n_side); end
    if is_odd
        y_side = cumsum(pitches);
        y_pos  = [-fliplr(y_side), 0, y_side];
    else
        y_side = (pitches(1)/2) + [0, cumsum(pitches(2:end))];
        y_pos  = [-fliplr(y_side), y_side];
    end
end

function y_pos = local_compute_parabolic_y_strict(N_cells, z_f, lambda, p_min)
    if mod(N_cells, 2) == 1
        idx = -(N_cells-1)/2 : (N_cells-1)/2;
    else
        idx = -(N_cells/2 - 0.5) : (N_cells/2 - 0.5);
    end
    y_pos = sign(idx) .* sqrt(2*abs(idx)*lambda*z_f + (abs(idx)*lambda).^2);
    [y_sorted, perm] = sort(y_pos);
    [~, i_anchor]    = min(abs(idx));
    i_anchor_sorted  = find(perm == i_anchor, 1);
    for i = i_anchor_sorted+1 : numel(y_sorted)
        if y_sorted(i) - y_sorted(i-1) < p_min, y_sorted(i) = y_sorted(i-1) + p_min; end
    end
    for i = i_anchor_sorted-1 : -1 : 1
        if y_sorted(i+1) - y_sorted(i) < p_min, y_sorted(i) = y_sorted(i+1) - p_min; end
    end
    y_pos(perm) = y_sorted;
end

function y_pos = local_compute_grouped_y(N_cells, p_min, z_focus, lambda, dr_frac)
    if dr_frac < 0 || dr_frac >= 0.5
        error('dr_frac must be in [0, 0.5)');
    end
    is_odd = mod(N_cells, 2) == 1;
    if is_odd
        n_side = (N_cells - 1) / 2;  y_curr = p_min;
    else
        n_side = N_cells / 2;        y_curr = p_min/2;
    end
    y_side = [];
    while numel(y_side) < n_side && y_curr <= 50e-3
        delta_r        = sqrt(y_curr^2 + z_focus^2) - z_focus;
        zone_real      = delta_r / lambda;
        dist_from_zone = abs(zone_real - round(zone_real));
        if dist_from_zone <= dr_frac
            y_side(end+1) = y_curr;  %#ok<AGROW>
        end
        y_curr = y_curr + p_min;
    end
    if is_odd, y_pos = sort([-y_side, 0, y_side]);
    else,      y_pos = sort([-y_side, y_side]);
    end
end

function y_pos = local_compute_hybrid_y(N_total, p_min, z_focus, lambda, ...
                                        N_fresnel, chirp_ratio)
    if mod(N_fresnel, 2) ~= 0, error('N_fresnel must be even.'); end
    if N_fresnel == 0
        y_pos = local_compute_chirped_pitch_y(N_total, p_min, chirp_ratio); return;
    end
    if N_fresnel >= N_total
        y_pos = local_compute_parabolic_y_strict(N_total, z_focus, lambda, p_min); return;
    end
    n_F_side = N_fresnel / 2;
    N_chirp  = N_total - N_fresnel;
    zones    = 1:n_F_side;
    y_F_side = sqrt(2*zones*lambda*z_focus + (zones*lambda).^2);
    y_chirp_max_allowed = y_F_side(1) - p_min;
    if y_chirp_max_allowed <= 0, error('No central space for chirp cells.'); end
    n_side_c = floor(N_chirp / 2);
    if N_chirp == 0
        y_chirp = [];
    else
        if mod(N_chirp, 2) == 0
            chirp_max_fit = 1 + 2*(y_chirp_max_allowed - (n_side_c-0.5)*p_min) / ...
                                 (n_side_c*p_min);
        else
            chirp_max_fit = 2*y_chirp_max_allowed/(n_side_c*p_min) - 1;
        end
        chirp_max_fit = max(1, chirp_max_fit);
        actual_chirp_ratio = min(chirp_ratio, chirp_max_fit);
        y_chirp = local_compute_chirped_pitch_y(N_chirp, p_min, actual_chirp_ratio);
    end
    y_pos = sort([y_chirp, -y_F_side, y_F_side]);
end

function [y_pos, y_label] = local_compute_multifocal_hybrid_y(N_total, p_min, ...
    lambda, z_focus_design_list, n_F_per_depth, chirp_ratio)
    if mod(n_F_per_depth, 2) ~= 0, error('n_F_per_depth must be even.'); end
    K = numel(z_focus_design_list);
    N_F_total = K * n_F_per_depth;
    if N_F_total > N_total
        error('Too many Fresnel cells (%d) for total cells (%d).', N_F_total, N_total);
    end
    n_F_side = n_F_per_depth / 2;
    y_F_all = []; label_F = [];
    for kk = 1:K
        z_k = z_focus_design_list(kk);
        zones = 1:n_F_side;
        y_F_k_side = sqrt(2*zones*lambda*z_k + (zones*lambda).^2);
        y_F_all = [y_F_all, -y_F_k_side, y_F_k_side]; %#ok<AGROW>
        label_F = [label_F, kk*ones(1, 2*n_F_side)];   %#ok<AGROW>
    end
    y_F_inner = min(abs(y_F_all));
    y_chirp_max_allowed = y_F_inner - p_min;
    if y_chirp_max_allowed <= 0, error('No central space for chirp.'); end
    N_chirp = N_total - N_F_total;
    if N_chirp == 0
        y_chirp = []; label_chirp = [];
    else
        n_side_c = floor(N_chirp/2);
        if mod(N_chirp, 2) == 0
            chirp_max_fit = 1 + 2*(y_chirp_max_allowed - (n_side_c-0.5)*p_min) / ...
                                 (n_side_c*p_min);
        else
            chirp_max_fit = 2*y_chirp_max_allowed/(n_side_c*p_min) - 1;
        end
        chirp_max_fit = max(1, chirp_max_fit);
        actual_chirp_ratio = min(chirp_ratio, chirp_max_fit);
        y_chirp = local_compute_chirped_pitch_y(N_chirp, p_min, actual_chirp_ratio);
        label_chirp = zeros(1, N_chirp);
    end
    y_combined     = [y_chirp, y_F_all];
    label_combined = [label_chirp, label_F];
    [y_pos, perm]  = sort(y_combined);
    y_label        = label_combined(perm);
    [~, ic] = min(abs(y_pos)); n = numel(y_pos);
    for i = ic+1:n
        if y_pos(i) - y_pos(i-1) < p_min, y_pos(i) = y_pos(i-1) + p_min; end
    end
    for i = ic-1:-1:1
        if y_pos(i+1) - y_pos(i) < p_min, y_pos(i) = y_pos(i+1) - p_min; end
    end
end

function on_axis = local_on_axis_profile(x_cell, y_cell, cell_tau, z_grid, k, fc)
    Nz = numel(z_grid);
    x_cell_s   = single(x_cell);
    y_cell_s   = single(y_cell);
    cell_tau_s = single(cell_tau);
    k_s        = single(k);
    omega_s    = single(2*pi*fc);
    P  = complex(zeros(1, Nz, 'single'));
    Nc = numel(x_cell);
    for i = 1:Nc
        r     = sqrt(x_cell_s(i)^2 + y_cell_s(i)^2 + z_grid.^2);
        phase = -(k_s * r + omega_s * cell_tau_s(i));
        P     = P + exp(1i*phase) ./ r;
    end
    on_axis = abs(P).^2;
end

function colors = local_make_colormap(N)
    t = linspace(0, 1, N)';
    colors = [0.15 + 0.80*t, 0.30 + 0.30*(1-abs(2*t-1)), 0.85 - 0.70*t];
end

% function [x_cell, y_cell, ich] = local_build_cells(y_in_ch, N_cell_x, ...
%     cell_pitch_x, N_ch, ch_pitch)
%     x_in_ch = ((1:N_cell_x) - (N_cell_x+1)/2) * cell_pitch_x;
%     N_cell_y = numel(y_in_ch);
%     [ix, iy, ich_grid] = ndgrid(1:N_cell_x, 1:N_cell_y, 1:N_ch);
%     ich = ich_grid(:);
%     ch_centre_x = (ich - (N_ch+1)/2) * ch_pitch;
%     dx_in_ch = x_in_ch(ix(:))';
%     dy_in_ch = y_in_ch(iy(:))';
%     x_cell   = ch_centre_x + dx_in_ch;
%     y_cell   = dy_in_ch;
% end

function [x_cell, y_cell, ich] = local_build_cells(y_in_ch, N_cell_x, ...
    cell_pitch_x, N_ch, ch_pitch)
    x_in_ch = ((1:N_cell_x) - (N_cell_x+1)/2) * cell_pitch_x;
    N_cell_y = numel(y_in_ch);
    [ix, iy, ich_grid] = ndgrid(1:N_cell_x, 1:N_cell_y, 1:N_ch);
    ich = ich_grid(:);
    ch_centre_x = (ich - (N_ch+1)/2) * ch_pitch;
    
    % FIXED: Removed the transpose (') from the next two lines
    dx_in_ch = x_in_ch(ix(:)); 
    dy_in_ch = y_in_ch(iy(:));
    
    x_cell   = ch_centre_x + dx_in_ch;
    y_cell   = dy_in_ch;
end
