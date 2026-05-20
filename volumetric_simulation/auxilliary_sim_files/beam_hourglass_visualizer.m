%% Beam Profile Visualizer  –  Hourglass, Fresnel Zones & Resolution Tradeoff
%
%  For each layout configuration this script computes and plots:
%    Fig 1 (per config) : 2-D elevational beam map (the "hourglass") with
%                         coloured zone boundaries, FWHM envelope, and two
%                         side panels: on-axis intensity & FWHM vs depth.
%    Fig 2              : Multi-config overlay of FWHM and intensity.
%    Fig 3              : Intensity vs FWHM scatter (one point per depth)
%                         — the key tradeoff diagram.
%    Fig 4              : Cell layouts (one row per config).
%
%  Physics recap
%  -------------
%  Electronic x-focusing is applied DYNAMICALLY at every depth z (simulating
%  real-time delay-and-sum in the azimuthal direction).  The elevational beam
%  is therefore governed purely by the cell arrangement within each channel,
%  which is the design variable you are optimising.
%
%  Zone boundaries (vertical dashed lines on every plot):
%    N_y = D_y^2 / (4*lambda)   — elevational near-field limit
%    N_x = D_x^2 / (4*lambda)   — azimuthal   near-field limit
%
%    z < N_y          Symmetric  zone  — both axes in near-field, tight spot
%    N_y < z < N_x    Astigmatic zone  — y diverges, x still focused
%    z > N_x          Divergent  zone  — both axes diverge, focus lost

clear; clc; close all;

%% ===================================================================
%%  ACOUSTIC & ARRAY PARAMETERS
%% ===================================================================
c_sound = 1480;              % speed of sound (m/s)
fc      = 10.5e6;            % centre frequency (Hz)
lambda  = c_sound / fc;      % wavelength (m)  ≈ 140.95 µm
k       = 2*pi / lambda;     % wavenumber (rad/m)
N_ch    = 64;                % number of addressable channels

%% ===================================================================
%%  CONFIGURATIONS  ── edit this block freely
%% ===================================================================
configs(1).name         = 'Equal pitch @ 75 µm ch pitch';
configs(1).type         = 'equal';
configs(1).cell_pitch_y = 75e-6;
configs(1).ch_pitch     = 75e-6;

configs(2).name         = 'Equal pitch @ 150 µm ch pitch';
configs(2).type         = 'equal';
configs(2).cell_pitch_y = 75e-6;
configs(2).ch_pitch     = 150e-6;

configs(3).name         = 'Grouped Fresnel z_{geom}=40 mm, dr=1/4';
configs(3).type         = 'grouped';
configs(3).z_focus_geom = 40e-3;
configs(3).dr_frac      = 1/4;
configs(3).ch_pitch     = 75e-6;

configs(4).name         = 'Grouped Fresnel z_{geom}=40 mm, dr=1/6';
configs(4).type         = 'grouped';
configs(4).z_focus_geom = 40e-3;
configs(4).dr_frac      = 1/6;
configs(4).ch_pitch     = 75e-6;

configs(5).name         = 'Chirp @ 150 µm ch pitch, ratio 4';
configs(5).type         = 'chirp';
configs(5).chirp_ratio  = 10; % 10 works
configs(5).N_cell_y     = 80; % 60 works
configs(5).ch_pitch     = 150e-6;

% ── fill in universal defaults ──────────────────────────────────────
for i = 1:numel(configs)
    if ~isfield(configs(i),'ch_pitch')     || isempty(configs(i).ch_pitch),     configs(i).ch_pitch     = 75e-6; end
    if ~isfield(configs(i),'N_cell_x')     || isempty(configs(i).N_cell_x),     configs(i).N_cell_x     = 1;     end
    if ~isfield(configs(i),'N_cell_y')     || isempty(configs(i).N_cell_y),     configs(i).N_cell_y     = 80;    end
    if ~isfield(configs(i),'cell_pitch_x') || isempty(configs(i).cell_pitch_x), configs(i).cell_pitch_x = 75e-6; end
    if ~isfield(configs(i),'p_min_y')      || isempty(configs(i).p_min_y),      configs(i).p_min_y      = 75e-6; end
end

%% ===================================================================
%%  COMPUTATION GRIDS  (reduce Nz/Ny if runtime is too long)
%% ===================================================================
z_max      = 80e-3;                        % maximum depth (m)
z_grid     = linspace(2e-3, z_max, 150);   % depth axis (150 pts)
y_range    = 8e-3;                         % ± half-width of elevation axis (m)
y_obs_grid = linspace(-y_range, y_range, 151);  % elevational observation axis
Nz         = numel(z_grid);
Ny         = numel(y_obs_grid);

% Pre-cast to single column for broadcasting
y_obs_col  = single(reshape(y_obs_grid, [], 1));  % Ny × 1

%% ===================================================================
%%  MAIN PROCESSING LOOP
%% ===================================================================
N_cfg   = numel(configs);
results = struct([]);

fprintf('Computing 2-D beam fields (%d configs × %d depths)…\n', N_cfg, Nz);
t_all = tic;

for c_idx = 1:N_cfg
    cfg = configs(c_idx);
    fprintf('  [%d/%d] %s\n', c_idx, N_cfg, cfg.name);
    t_cfg = tic;

    % ── build cell layout ──────────────────────────────────────────
    [y_in_ch, ~]          = local_y_in_channel(cfg, lambda);
    [x_cell, y_cell, ich] = local_build_cells(y_in_ch, cfg.N_cell_x, ...
                                cfg.cell_pitch_x, N_ch, cfg.ch_pitch);

    x_cell_s = single(x_cell(:)');   % 1 × Nc
    y_cell_s = single(y_cell(:)');   % 1 × Nc
    Nc       = numel(x_cell_s);

    % Channel x-centres for delay law
    ch_x_ctrs = single(((1:N_ch) - (N_ch+1)/2) * cfg.ch_pitch);  % 1 × N_ch

    % ── aperture dimensions & near-field distances ──────────────────
    D_y    = max(y_cell) - min(y_cell);
    D_x    = max(x_cell) - min(x_cell);
    N_y_nf = D_y^2 / (4*lambda);   % elevational natural near-field (m)
    N_x_nf = D_x^2 / (4*lambda);   % azimuthal   natural near-field (m)

    % ── 2-D field: dynamically focused in x, natural in y ──────────
    %  For each depth z:
    %    1. Compute channel x-delays that focus to (x=0, z).
    %    2. Propagate from all cells to every (x_obs=0, y_obs, z) point.
    %  r(iy, cell) = sqrt( x_cell^2 + (y_cell − y_obs)^2 + z^2 )
    field2D  = zeros(Ny, Nz, 'single');
    k_s      = single(k);
    omega_s  = single(2*pi*fc);
    c_s      = single(c_sound);

    for iz = 1:Nz
        z = single(z_grid(iz));

        % Azimuthal (x) dynamic delay-and-sum ── focuses to (0, 0, z)
        d_ch_x  = sqrt(ch_x_ctrs.^2 + z^2);           % 1 × N_ch
        ch_tau  = (max(d_ch_x) - d_ch_x) / c_s;        % 1 × N_ch
        tau_s   = single(ch_tau(ich(:)'));               % 1 × Nc

        % Propagation distances — broadcast: (Ny×1) − (1×Nc) → Ny×Nc
        dy  = y_obs_col - y_cell_s;                     % Ny × Nc
        r   = sqrt(dy.^2 + (x_cell_s.^2 + z^2));       % Ny × Nc

        % Coherent pressure field and intensity
        phase        = k_s .* r + omega_s .* tau_s;    % Ny × Nc
        P            = sum(exp(1i .* phase) ./ r, 2);   % Ny × 1
        field2D(:,iz) = abs(P).^2;
    end

    % ── on-axis intensity (y_obs ≈ 0) ──────────────────────────────
    [~, iy0]  = min(abs(y_obs_grid));
    on_axis_I = double(field2D(iy0, :));   % 1 × Nz

    % ── elevational FWHM at each depth (robust peak-centred search) ─
    fwhm_mm = nan(1, Nz);
    for iz = 1:Nz
        prof = double(field2D(:, iz));
        [peak_val, pk_idx] = max(prof);
        if peak_val == 0, continue; end
        half = peak_val * 0.5;                     % −3 dB (half-power)
        L = find(prof(1:pk_idx) < half, 1, 'last');      if isempty(L), L = 1;  end
        R = find(prof(pk_idx:end) < half, 1, 'first');
        if isempty(R), R = Ny; else, R = R + pk_idx - 1; end
        fwhm_mm(iz) = (y_obs_grid(R) - y_obs_grid(L)) * 1e3;
    end

    % ── store ───────────────────────────────────────────────────────
    results(c_idx).cfg     = cfg;
    results(c_idx).field2D = field2D;
    results(c_idx).on_axis = on_axis_I;
    results(c_idx).fwhm_mm = fwhm_mm;
    results(c_idx).D_y     = D_y;
    results(c_idx).D_x     = D_x;
    results(c_idx).N_y_nf  = N_y_nf;
    results(c_idx).N_x_nf  = N_x_nf;
    results(c_idx).y_in_ch = y_in_ch;
    fprintf('      Nc = %d  |  D_y = %.2f mm  |  N_y = %.1f mm  |  N_x = %.1f mm  [%.1f s]\n', ...
            Nc, D_y*1e3, N_y_nf*1e3, N_x_nf*1e3, toc(t_cfg));
end
fprintf('Total time: %.1f s\n\n', toc(t_all));

%% ===================================================================
%%  GLOBAL NORMALISATION (common dB reference)
%% ===================================================================
maxI_global = 0;
for i = 1:N_cfg
    maxI_global = max(maxI_global, max(results(i).on_axis));
end

%% ===================================================================
%%  COLOUR PALETTE
%% ===================================================================
clr_cfg   = lines(N_cfg);
clr_sym   = [0.10 0.72 0.10];   % green  – symmetric zone
clr_astig = [0.95 0.58 0.00];   % amber  – astigmatic zone
clr_div   = [0.85 0.12 0.12];   % red    – divergent zone
db_floor  = -30;                 % dB floor for colour maps
z_mm      = z_grid * 1e3;        % depth in mm (plot axis)
y_mm      = y_obs_grid * 1e3;    % elevation in mm (plot axis)

%% ===================================================================
%%  FIGURE 1 (per config): HOURGLASS MAP + SIDE PANELS
%% ===================================================================
for c_idx = 1:N_cfg
    r       = results(c_idx);
    cfg     = r.cfg;
    fdb     = 10*log10(double(r.field2D) / maxI_global + eps);
    Ny_mm   = r.N_y_nf * 1e3;
    Nx_mm   = r.N_x_nf * 1e3;
    z_max_mm = z_max * 1e3;

    fig = figure('Color','w','Position',[50 50 1400 600], ...
                 'Name', ['Hourglass | ' cfg.name]);

    % ── Left: 2-D hourglass map ──────────────────────────────────
    ax1 = axes('Position',[0.05 0.11 0.53 0.78]);
    imagesc(ax1, z_mm, y_mm, fdb);
    set(ax1,'YDir','normal');
    colormap(ax1, local_hot_transparent());  % dark = low, bright = high
    clim(ax1, [db_floor 0]);
    cb = colorbar(ax1);
    cb.Label.String = 'Intensity  (dB, common ref)';
    cb.Label.FontSize = 9;
    hold(ax1,'on');

    % ── Zone shading (semi-transparent rectangles painted over image) ─
    yl = [y_mm(1) y_mm(end)];
    local_zone_patch(ax1, z_mm(1),           min(Ny_mm, z_max_mm), yl, clr_sym,   0.13);
    local_zone_patch(ax1, min(Ny_mm,z_max_mm), min(Nx_mm,z_max_mm), yl, clr_astig, 0.13);
    local_zone_patch(ax1, min(Nx_mm,z_max_mm), z_max_mm,             yl, clr_div,   0.13);

    % ── Zone boundary lines ──────────────────────────────────────
    if Ny_mm < z_max_mm
        xline(ax1, Ny_mm, '--', 'Color', clr_sym,   'LineWidth', 2.2, ...
              'DisplayName', sprintf('N_y = %.1f mm', Ny_mm));
    end
    if Nx_mm < z_max_mm
        xline(ax1, Nx_mm, '--', 'Color', clr_div,   'LineWidth', 2.2, ...
              'DisplayName', sprintf('N_x = %.1f mm', Nx_mm));
    end

    % ── FWHM envelope (beam edge contour, −3 dB) ────────────────
    hw   = r.fwhm_mm / 2;
    good = ~isnan(hw);
    plot(ax1, z_mm(good),  hw(good), 'c-', 'LineWidth', 2.2, ...
         'DisplayName', '±FWHM/2  (−3 dB)');
    plot(ax1, z_mm(good), -hw(good), 'c-', 'LineWidth', 2.2, ...
         'HandleVisibility','off');

    % ── Zone labels ──────────────────────────────────────────────
    txt_y = y_mm(end) * 0.87;
    local_zone_text(ax1, z_mm(1),          Ny_mm,    z_max_mm, txt_y, 'Sym.',      clr_sym,   10);
    local_zone_text(ax1, Ny_mm,            Nx_mm,    z_max_mm, txt_y, 'Astig.',    clr_astig, 10);
    local_zone_text(ax1, Nx_mm,            z_max_mm, z_max_mm, txt_y, 'Divergent', clr_div,   10);

    legend(ax1,'Location','northeast','FontSize',7,'Box','off');
    xlabel(ax1,'Depth  z  (mm)');
    ylabel(ax1,'Elevational position  y_{elev}  (mm)');
    title(ax1, {sprintf('Elevational Beam Map — %s', cfg.name), ...
                sprintf('D_y=%.2f mm  |  D_x=%.2f mm  |  N_y=%.1f mm  |  N_x=%.1f mm', ...
                        r.D_y*1e3, r.D_x*1e3, Ny_mm, Nx_mm)}, ...
          'Interpreter','none','FontSize',9);

    % ── Top-right: on-axis intensity ────────────────────────────
    ax2 = axes('Position',[0.63 0.55 0.34 0.35]);
    on_dB = 10*log10(r.on_axis / maxI_global + eps);
    plot(ax2, z_mm, on_dB, 'k-', 'LineWidth', 1.8);
    hold(ax2,'on'); grid(ax2,'on');
    local_zone_xlines(ax2, Ny_mm, Nx_mm, z_max_mm, clr_sym, clr_div, 1.5);
    ylim(ax2,[db_floor 3]); xlim(ax2,[z_mm(1) z_max_mm]);
    xlabel(ax2,'Depth  (mm)'); ylabel(ax2,'Intensity  (dB)');
    title(ax2,'On-axis intensity  ↑ = more power','FontSize',9);
    legend(ax2,'Location','southwest','FontSize',7,'Box','off');

    % ── Bottom-right: FWHM (resolution proxy) ──────────────────
    ax3 = axes('Position',[0.63 0.11 0.34 0.35]);
    plot(ax3, z_mm, r.fwhm_mm, '-','Color',[0.12 0.47 0.87],'LineWidth',1.8);
    hold(ax3,'on'); grid(ax3,'on');
    local_zone_xlines(ax3, Ny_mm, Nx_mm, z_max_mm, clr_sym, clr_div, 1.5);
    xlim(ax3,[z_mm(1) z_max_mm]);
    xlabel(ax3,'Depth  (mm)'); ylabel(ax3,'FWHM  (mm)');
    title(ax3,'Elevational FWHM  ↓ = better resolution','FontSize',9);
end

%% ===================================================================
%%  FIGURE 2: Multi-config FWHM and Intensity overlay
%% ===================================================================
fig2 = figure('Color','w','Position',[80 80 1300 520], ...
              'Name','Comparison — FWHM & Intensity (all configs)');

ax_fw = subplot(1,2,1,'Parent',fig2);
hold(ax_fw,'on'); grid(ax_fw,'on');
for i = 1:N_cfg
    plot(ax_fw, z_mm, results(i).fwhm_mm, '-', ...
         'Color', clr_cfg(i,:), 'LineWidth', 2.0, ...
         'DisplayName', results(i).cfg.name);
end
% find the tightest common N_y_nf for a reference line
Ny_ref = results(1).N_y_nf*1e3; Nx_ref = results(1).N_x_nf*1e3;
xline(ax_fw, Ny_ref, ':k','LineWidth',1.2,'HandleVisibility','off');
xline(ax_fw, Nx_ref, ':k','LineWidth',1.2,'HandleVisibility','off');
xlabel(ax_fw,'Depth  z  (mm)');
ylabel(ax_fw,'Elevational FWHM  (mm)');
title(ax_fw,{'Resolution comparison — all configs', ...
             '↓ smaller FWHM = tighter focus = better resolution'},'FontSize',9);
legend(ax_fw,'Interpreter','none','Location','northwest','FontSize',8);

ax_in = subplot(1,2,2,'Parent',fig2);
hold(ax_in,'on'); grid(ax_in,'on');
for i = 1:N_cfg
    on_dB = 10*log10(results(i).on_axis / maxI_global + eps);
    plot(ax_in, z_mm, on_dB, '-', ...
         'Color', clr_cfg(i,:), 'LineWidth', 2.0, ...
         'DisplayName', results(i).cfg.name);
end
xline(ax_in, Ny_ref, ':k','LineWidth',1.2,'HandleVisibility','off');
xline(ax_in, Nx_ref, ':k','LineWidth',1.2,'HandleVisibility','off');
ylim(ax_in,[db_floor 5]);
xlabel(ax_in,'Depth  z  (mm)');
ylabel(ax_in,'On-axis intensity  (dB, common ref)');
title(ax_in,{'Intensity comparison — all configs', ...
             '↑ higher = more acoustic power at that depth'},'FontSize',9);
legend(ax_in,'Interpreter','none','Location','southwest','FontSize',8);

%% ===================================================================
%%  FIGURE 3: Tradeoff scatter — Intensity vs FWHM (coloured by depth)
%%
%%  Each point represents one steered depth.
%%  Upper-left corner = high power AND fine resolution = ideal.
%%  A config that traces further upper-left is strictly better.
%% ===================================================================
fig3 = figure('Color','w','Position',[100 100 980 620], ...
              'Name','Tradeoff: Intensity vs. Resolution');
ax_td = axes(fig3);
hold(ax_td,'on'); grid(ax_td,'on');

depth_cmap = parula(256);
z_norm     = (z_mm - z_mm(1)) / (z_mm(end) - z_mm(1));  % 0‥1

for c_idx = 1:N_cfg
    r     = results(c_idx);
    on_dB = 10*log10(r.on_axis / maxI_global + eps);
    valid = ~isnan(r.fwhm_mm);

    % thin connecting line
    plot(ax_td, r.fwhm_mm(valid), on_dB(valid), '-', ...
         'Color', [clr_cfg(c_idx,:) 0.35], 'LineWidth', 1.2, ...
         'HandleVisibility','off');

    % scatter points coloured by depth
    scatter(ax_td, r.fwhm_mm(valid), on_dB(valid), 55, ...
            z_mm(valid), 'filled', ...
            'MarkerEdgeColor', clr_cfg(c_idx,:), 'LineWidth', 1.2, ...
            'DisplayName', r.cfg.name);
end

% Annotate with depth markers at a few representative depths
anno_depths_mm = [10 20 30 40 60];
for c_idx = 1
    r     = results(c_idx);
    on_dB = 10*log10(r.on_axis / maxI_global + eps);
    for d = anno_depths_mm
        [~,iz] = min(abs(z_mm - d));
        if ~isnan(r.fwhm_mm(iz))
            text(ax_td, r.fwhm_mm(iz)+0.05, on_dB(iz), ...
                 sprintf('%g mm', d), 'FontSize', 7, 'Color', [0.4 0.4 0.4]);
        end
    end
end

cb3 = colorbar(ax_td);
cb3.Label.String = 'Depth z  (mm)';
colormap(ax_td, parula);
clim(ax_td, [z_mm(1) z_mm(end)]);
xlabel(ax_td, 'Elevational FWHM  (mm)   ←  narrower = better resolution');
ylabel(ax_td, 'On-axis intensity  (dB)   ↑  higher = more power');
title(ax_td, {'Intensity–Resolution Tradeoff  (each point = one steered depth)', ...
              'Upper-left corner = best of both worlds'}, 'FontSize',10);
legend(ax_td,'Interpreter','none','Location','southeast','FontSize',8);

%% ===================================================================
%%  FIGURE 4: Cell layouts within one channel
%% ===================================================================
fig4 = figure('Color','w','Position',[60 60 1200 80+65*N_cfg], ...
              'Name','Cell layouts within one channel');
ax_lay = axes(fig4);
hold(ax_lay,'on'); grid(ax_lay,'on');
for i = 1:N_cfg
    plot(ax_lay, results(i).y_in_ch*1e3, ...
         i*ones(size(results(i).y_in_ch)), 'o', ...
         'MarkerSize',5,'MarkerFaceColor',clr_cfg(i,:),'Color',clr_cfg(i,:));
    % annotate aperture span
    D = (max(results(i).y_in_ch) - min(results(i).y_in_ch))*1e3;
    text(ax_lay, max(results(i).y_in_ch)*1e3+0.1, i, ...
         sprintf('D_y=%.2f mm', D), 'FontSize',7, 'Color', clr_cfg(i,:));
end
xline(ax_lay,0,'k:','LineWidth',1.2);
yticks(ax_lay,1:N_cfg);
yticklabels(ax_lay, arrayfun(@(r) r.cfg.name, results, 'UniformOutput',false));
ylim(ax_lay,[0.5 N_cfg+0.5]);
xlabel(ax_lay,'Cell y-position within channel  (mm)');
title(ax_lay,'Cell layouts — one row per configuration','Interpreter','none');
set(ax_lay,'TickLabelInterpreter','none');

fprintf('All figures generated.\n');

%% ===================================================================
%%  HELPER FUNCTIONS
%% ===================================================================

% ── Zone shading patch (drawn OVER imagesc; uistack to bottom) ────
function local_zone_patch(ax, x1, x2, yl, clr, alpha)
    if x2 <= x1, return; end
    h = fill(ax, [x1 x2 x2 x1], [yl(1) yl(1) yl(2) yl(2)], ...
             clr, 'FaceAlpha', alpha, 'EdgeColor', 'none', ...
             'HandleVisibility','off');
    uistack(h,'bottom');          % push behind image
end

% ── Zone boundary xlines (for side panels) ────────────────────────
function local_zone_xlines(ax, Ny_mm, Nx_mm, z_max, c1, c2, lw)
    if Ny_mm < z_max
        xline(ax, Ny_mm, '--', 'Color', c1, 'LineWidth', lw, ...
              'DisplayName', sprintf('N_y=%.1fmm', Ny_mm));
    end
    if Nx_mm < z_max
        xline(ax, Nx_mm, '--', 'Color', c2, 'LineWidth', lw, ...
              'DisplayName', sprintf('N_x=%.1fmm', Nx_mm));
    end
end

% ── Zone label (centred in its zone band if visible) ─────────────
function local_zone_text(ax, x_start, x_end, x_max, y_pos, lbl, clr, fsz)
    x1 = min(x_start, x_max);
    x2 = min(x_end,   x_max);
    xm = (x1 + x2) / 2;
    if x2 > x1 && xm > ax.XLim(1)
        text(ax, xm, y_pos, lbl, 'Color', clr, 'FontWeight','bold', ...
             'FontSize', fsz, 'HorizontalAlignment','center', ...
             'VerticalAlignment','top');
    end
end

% ── Custom hot colourmap (black→red→orange→white) ─────────────────
function cmap = local_hot_transparent()
    cmap = hot(256);
    % darken the bottom end slightly so the background stays black
    cmap(1:20,:) = cmap(1:20,:) .* linspace(0, 1, 20)';
end

% ── Cell layout builder (fixed version from sweep script) ─────────
function [x_cell, y_cell, ich] = local_build_cells(y_in_ch, N_cell_x, ...
        cell_pitch_x, N_ch, ch_pitch)
    x_in_ch  = ((1:N_cell_x) - (N_cell_x+1)/2) * cell_pitch_x;
    N_cell_y = numel(y_in_ch);
    [ix, iy, ich_grid] = ndgrid(1:N_cell_x, 1:N_cell_y, 1:N_ch);
    ich = ich_grid(:);
    ch_centre_x = (ich - (N_ch+1)/2) * ch_pitch;
    x_cell = ch_centre_x + x_in_ch(ix(:));
    y_cell = y_in_ch(iy(:));
end

% ── Cell y-positions dispatcher ───────────────────────────────────
function [y_in_ch, y_label] = local_y_in_channel(cfg, lambda)
    y_label = [];
    switch cfg.type
        case 'equal'
            y_in_ch = ((1:cfg.N_cell_y) - (cfg.N_cell_y+1)/2) * cfg.cell_pitch_y;
        case 'grouped'
            y_in_ch = local_compute_grouped_y(cfg.N_cell_y, cfg.p_min_y, ...
                          cfg.z_focus_geom, lambda, cfg.dr_frac);
        case 'chirp'
            y_in_ch = local_compute_chirped_pitch_y(cfg.N_cell_y, ...
                          cfg.p_min_y, cfg.chirp_ratio);
        case 'parabolic'
            y_in_ch = local_compute_parabolic_y_strict(cfg.N_cell_y, ...
                          cfg.z_focus_geom, lambda, cfg.p_min_y);
        case 'hybrid'
            y_in_ch = local_compute_hybrid_y(cfg.N_cell_y, cfg.p_min_y, ...
                          cfg.z_focus_geom, lambda, cfg.N_fresnel, cfg.chirp_ratio);
        case 'multifocal_hybrid'
            [y_in_ch, y_label] = local_compute_multifocal_hybrid_y( ...
                          cfg.N_cell_y, cfg.p_min_y, lambda, ...
                          cfg.z_focus_design_list, cfg.n_F_per_depth, cfg.chirp_ratio);
        otherwise
            error('Unknown cfg.type: ''%s''', cfg.type);
    end
end

% ── Grouped Fresnel layout ────────────────────────────────────────
function y_pos = local_compute_grouped_y(N_cells, p_min, z_focus, lambda, dr_frac)
    if dr_frac < 0 || dr_frac >= 0.5, error('dr_frac must be in [0, 0.5)'); end
    is_odd = mod(N_cells,2)==1;
    if is_odd, n_side = (N_cells-1)/2; y_curr = p_min;
    else,       n_side = N_cells/2;     y_curr = p_min/2; end
    y_side = [];
    while numel(y_side) < n_side && y_curr <= 50e-3
        delta_r        = sqrt(y_curr^2 + z_focus^2) - z_focus;
        dist_from_zone = abs(delta_r/lambda - round(delta_r/lambda));
        if dist_from_zone <= dr_frac, y_side(end+1) = y_curr; end  %#ok<AGROW>
        y_curr = y_curr + p_min;
    end
    if is_odd, y_pos = sort([-y_side, 0, y_side]);
    else,      y_pos = sort([-y_side,    y_side]); end
end

% ── Chirped pitch layout ──────────────────────────────────────────
function y_pos = local_compute_chirped_pitch_y(N_cells, p_min, chirp_ratio)
    if N_cells<=0, y_pos=[]; return; end
    if N_cells==1, y_pos=0;  return; end
    is_odd = mod(N_cells,2)==1;
    n_side = floor(N_cells/2);
    if n_side==1, pitches = p_min;
    else,          pitches = p_min * linspace(1, chirp_ratio, n_side); end
    if is_odd
        y_side = cumsum(pitches);
        y_pos  = [-fliplr(y_side), 0, y_side];
    else
        y_side = (pitches(1)/2) + [0, cumsum(pitches(2:end))];
        y_pos  = [-fliplr(y_side), y_side];
    end
end

% ── Parabolic (strict Fresnel) layout ─────────────────────────────
function y_pos = local_compute_parabolic_y_strict(N_cells, z_f, lambda, p_min)
    if mod(N_cells,2)==1, idx=-(N_cells-1)/2:(N_cells-1)/2;
    else, idx=-(N_cells/2-0.5):(N_cells/2-0.5); end
    y_pos = sign(idx) .* sqrt(2*abs(idx)*lambda*z_f + (abs(idx)*lambda).^2);
    [ys,perm] = sort(y_pos);
    [~,ia] = min(abs(idx)); ias = find(perm==ia,1);
    for i=ias+1:numel(ys), if ys(i)-ys(i-1)<p_min, ys(i)=ys(i-1)+p_min; end; end
    for i=ias-1:-1:1,       if ys(i+1)-ys(i)<p_min, ys(i)=ys(i+1)-p_min; end; end
    y_pos(perm)=ys;
end

% ── Hybrid (Fresnel outer + chirp inner) layout ───────────────────
function y_pos = local_compute_hybrid_y(N_total, p_min, z_focus, lambda, N_fresnel, chirp_ratio)
    if mod(N_fresnel,2)~=0, error('N_fresnel must be even.'); end
    if N_fresnel==0,      y_pos=local_compute_chirped_pitch_y(N_total,p_min,chirp_ratio); return; end
    if N_fresnel>=N_total, y_pos=local_compute_parabolic_y_strict(N_total,z_focus,lambda,p_min); return; end
    n_F_side   = N_fresnel/2;
    y_F_side   = sqrt(2*(1:n_F_side)*lambda*z_focus + ((1:n_F_side)*lambda).^2);
    y_max_c    = y_F_side(1)-p_min;
    if y_max_c<=0, error('No central space for chirp cells.'); end
    N_chirp    = N_total - N_fresnel;
    if N_chirp==0, y_chirp=[]; else
        nc = floor(N_chirp/2);
        if mod(N_chirp,2)==0, cmf=1+2*(y_max_c-(nc-0.5)*p_min)/(nc*p_min);
        else, cmf=2*y_max_c/(nc*p_min)-1; end
        y_chirp = local_compute_chirped_pitch_y(N_chirp,p_min,min(chirp_ratio,max(1,cmf)));
    end
    y_pos = sort([y_chirp,-y_F_side,y_F_side]);
end

% ── Multi-focal hybrid layout ─────────────────────────────────────
function [y_pos, y_label] = local_compute_multifocal_hybrid_y(N_total, p_min, lambda, ...
        z_list, n_F_per_depth, chirp_ratio)
    if mod(n_F_per_depth,2)~=0, error('n_F_per_depth must be even.'); end
    K=numel(z_list); NF=K*n_F_per_depth;
    if NF>N_total, error('Too many Fresnel cells (%d) for total cells (%d).',NF,N_total); end
    nfs=n_F_per_depth/2; y_F_all=[]; label_F=[];
    for kk=1:K
        yFk=sqrt(2*(1:nfs)*lambda*z_list(kk)+((1:nfs)*lambda).^2);
        y_F_all=[y_F_all,-yFk,yFk]; label_F=[label_F,kk*ones(1,2*nfs)]; %#ok<AGROW>
    end
    ymax_c=min(abs(y_F_all))-p_min;
    if ymax_c<=0, error('No central space for chirp.'); end
    NC=N_total-NF;
    if NC==0, y_chirp=[]; label_chirp=[];
    else
        nc=floor(NC/2);
        if mod(NC,2)==0, cmf=1+2*(ymax_c-(nc-0.5)*p_min)/(nc*p_min);
        else, cmf=2*ymax_c/(nc*p_min)-1; end
        y_chirp=local_compute_chirped_pitch_y(NC,p_min,min(chirp_ratio,max(1,cmf)));
        label_chirp=zeros(1,NC);
    end
    [y_pos,perm]=sort([y_chirp,y_F_all]);
    y_label=[label_chirp,label_F]; y_label=y_label(perm);
    [~,ic]=min(abs(y_pos)); n=numel(y_pos);
    for i=ic+1:n,    if y_pos(i)-y_pos(i-1)<p_min, y_pos(i)=y_pos(i-1)+p_min; end; end
    for i=ic-1:-1:1, if y_pos(i+1)-y_pos(i)<p_min, y_pos(i)=y_pos(i+1)-p_min; end; end
end
