%% 3D volumetric beam visualization for PMUT array
%
% Renders the volumetric intensity field for any layout type:
%   - 'equal'     : fixed-pitch
%   - 'parabolic' : strict Fresnel (1 cell/zone)
%   - 'grouped'   : grouped Fresnel (relaxed-coherence)
%
% Each configuration produces ONE 3D figure containing:
%   1. Translucent orthogonal slice planes at x = 0, y = 0, z = z_focus
%      (these reveal the internal structure of the beam)
%   2. Nested isosurfaces at -3, -6, -12 dB
%      (these show the 3D shape and taper of the main lobe)
%   3. The array footprint plotted at z = 0
%      (so you can see geometry-vs-beam at a glance)
%   4. A green crosshair at the requested focus
%
% Grid order note: the sim uses ndgrid (Nx, Ny, Nz). MATLAB's isosurface
% and slice expect meshgrid (Ny, Nx, Nz). The viz function permutes once.

clear; clc; close all;

%% ---------- 1. Acoustic parameters ------------------------------------
c       = 1480;
fc      = 10.5e6;
lambda  = c / fc;
k       = 2*pi / lambda;

%% ---------- 2. Common channel layout ----------------------------------
N_ch        = 64;
ch_pitch    = 75e-6;

%% ---------- 3. Transmit (electronic) focus law -----------------------
mode    = 'focus';
x_focus = 0;
z_focus = 40e-3;

ch_x_centres = ((1:N_ch) - (N_ch+1)/2) * ch_pitch;
switch mode
    case 'plane'
        ch_tau = zeros(1, N_ch);
    case 'focus'
        d_ch    = sqrt((x_focus - ch_x_centres).^2 + z_focus^2);
        ch_tau  = (max(d_ch) - d_ch) / c;
    otherwise
        error('mode must be ''plane'' or ''focus''.');
end

%% ---------- 4. Imaging volume ----------------------------------------
% Higher resolution and a tighter z-range than the benchmark scripts:
% smooth isosurfaces want fine voxels, and we don't need 220 mm of depth
% for visualization (most of the structure lives within ~2 z_focus).
Nx = 81;  Ny = 81;  Nz = 160;
x_grid = single(linspace(-6e-3,  6e-3, Nx));
y_grid = single(linspace(-6e-3,  6e-3, Ny));
z_grid = single(linspace( 1e-3,  90e-3, Nz));   % 1 mm to ~2.25 * z_focus
[X3, Y3, Z3] = ndgrid(x_grid, y_grid, z_grid);

%% ---------- 5. Configurations to render ------------------------------
% Comment out / add freely. Each one becomes its own 3D figure.

cfg(1).name         = 'Equal pitch 75um, 80 cells/ch';
cfg(1).type         = 'equal';
cfg(1).N_cell_x     = 1;
cfg(1).N_cell_y     = 80;
cfg(1).cell_pitch_x = 75e-6;
cfg(1).cell_pitch_y = 75e-6;

% cfg(2).name         = 'Strict Fresnel, 80 cells/ch, focus = 40 mm';
% cfg(2).type         = 'parabolic';
% cfg(2).N_cell_x     = 1;
% cfg(2).N_cell_y     = 80;
% cfg(2).cell_pitch_x = 75e-6;
% cfg(2).p_min_y      = 75e-6;
% cfg(2).z_focus_geom = z_focus;
% 
% cfg(3).name         = 'Grouped Fresnel, dr = lambda/4, 80 cells/ch';
% cfg(3).type         = 'grouped';
% cfg(3).N_cell_x     = 1;
% cfg(3).N_cell_y     = 80;
% cfg(3).cell_pitch_x = 75e-6;
% cfg(3).p_min_y      = 75e-6;
% cfg(3).z_focus_geom = z_focus;
% cfg(3).dr_frac      = 1/4;

N_cfg = numel(cfg);

%% ---------- 6. Visualisation options ---------------------------------
viz_opts.iso_levels = [-3, -6, -12];                  % dB
viz_opts.iso_colors = [1.00 0.90 0.30;                % yellow (innermost)
                       1.00 0.45 0.10;                % orange (middle)
                       0.80 0.15 0.10];               % red    (outermost)
viz_opts.iso_alpha  = [0.85, 0.45, 0.18];
viz_opts.slice_alpha = 0.55;
viz_opts.colormap   = parula(256);
viz_opts.clim_dB    = [-25 0];

%% ---------- 7. Loop: build cells, run sim, render --------------------
for j = 1:N_cfg
    fprintf('\n===== Config %d/%d: %s =====\n', j, N_cfg, cfg(j).name);

    [x_cell, y_cell, ich, group_y] = local_build_cells(cfg(j), N_ch, ch_pitch, lambda);
    cell_tau = ch_tau(ich);
    fprintf('  N_cells=%d   y range = [%.2f, %.2f] mm\n', ...
            numel(x_cell), min(y_cell)*1e3, max(y_cell)*1e3);

    P3   = local_compute_field(x_cell, y_cell, cell_tau, X3, Y3, Z3, k, fc);
    I3   = abs(P3).^2;
    I3dB = 10*log10(I3 / max(I3(:)) + eps);

    local_viz_beam_3d(I3dB, x_grid, y_grid, z_grid, ...
                      x_cell, y_cell, group_y, ...
                      x_focus, z_focus, cfg(j).name, viz_opts);
end


%% ====================== HELPER FUNCTIONS =============================
function local_viz_beam_3d(I3dB, xg, yg, zg, x_cell, y_cell, group_y, ...
                            x_focus, z_focus, name, opts)
    % Render one 3D figure: slices + isosurfaces + array + focal point.

    % ndgrid (Nx,Ny,Nz)  ->  meshgrid (Ny,Nx,Nz). isosurface/slice want meshgrid.
    V  = permute(double(I3dB), [2 1 3]);
    xg = double(xg)*1e3;          % all axes in mm
    yg = double(yg)*1e3;
    zg = double(zg)*1e3;
    zfm = z_focus*1e3;
    xfm = x_focus*1e3;

    figure('Color','w','Position',[80 80 1180 820]);
    hold on;

    % --- (a) Translucent slice planes through the focus ---
    hs = slice(xg, yg, zg, V, xfm, 0, zfm);
    set(hs, 'EdgeColor','none', 'FaceAlpha', opts.slice_alpha);
    colormap(opts.colormap);
    clim(opts.clim_dB);

    % --- (b) Nested isosurfaces (outer first so inner draws on top) ---
    iso_levels = opts.iso_levels;
    iso_colors = opts.iso_colors;
    iso_alpha  = opts.iso_alpha;
    [iso_levels, ord] = sort(iso_levels, 'ascend');     % most-negative first
    iso_colors = iso_colors(ord, :);
    iso_alpha  = iso_alpha(ord);

    iso_handles = gobjects(0);
    iso_labels  = {};
    for ii = 1:numel(iso_levels)
        lvl = iso_levels(ii);
        fv  = isosurface(xg, yg, zg, V, lvl);
        if isempty(fv.faces), continue; end
        h = patch(fv, 'FaceColor', iso_colors(ii,:), 'EdgeColor','none', ...
                  'FaceAlpha', iso_alpha(ii), ...
                  'AmbientStrength', 0.4, 'DiffuseStrength', 0.7, ...
                  'SpecularStrength', 0.15);
        iso_handles(end+1) = h;                         %#ok<AGROW>
        iso_labels{end+1}  = sprintf('%g dB', lvl);     %#ok<AGROW>
    end

    % --- (c) Array elements at z = 0 ---
    plot3(x_cell*1e3, y_cell*1e3, zeros(size(x_cell)), '.', ...
          'Color', [0.25 0.25 0.25], 'MarkerSize', 3);

    % --- (d) Focal point ---
    plot3(xfm, 0, zfm, 'g+', 'MarkerSize', 22, 'LineWidth', 2.5);
    plot3(xfm, 0, zfm, 'go', 'MarkerSize', 14, 'LineWidth', 1.8);

    % --- Cosmetic axes ---
    xlabel('Azimuth x (mm)');
    ylabel('Elevation y (mm)');
    zlabel('Depth z (mm)');
    title(name, 'Interpreter','none', 'FontSize', 12);

    cb = colorbar('Location','eastoutside');
    ylabel(cb, 'Intensity (dB, slice rendering)');

    if ~isempty(iso_handles)
        legend(iso_handles, iso_labels, 'Location','northeast', ...
               'AutoUpdate','off');
    end

    axis tight;
    grid on;
    box on;
    view(-37.5, 22);

    % Lighting: makes isosurfaces 3D-looking
    camlight headlight;
    camlight(-30, -30);
    lighting gouraud;
    material dull;

    drawnow;
end

function [y_pos, group_idx] = local_compute_grouped_fresnel_y( ...
                                  N_cells, z_f, lambda, p_min, dr_max)
    y_list = []; g_list = [];
    g = 0; max_groups = 10000;
    while numel(y_list) < N_cells
        r_anchor  = z_f + g*lambda;
        r_max     = r_anchor + dr_max;
        y_start_g = sqrt(r_anchor^2 - z_f^2);
        y_end_g   = sqrt(r_max^2    - z_f^2);
        if g == 0
            if mod(N_cells,2) == 1
                y_pos_side = 0:p_min:y_end_g;
                if isempty(y_pos_side), y_pos_side = 0; end
            else
                y_pos_side = (p_min/2):p_min:y_end_g;
                if isempty(y_pos_side), y_pos_side = p_min/2; end
            end
        else
            y_pos_side = y_start_g:p_min:y_end_g;
            if isempty(y_pos_side), y_pos_side = y_start_g; end
        end
        for i = 1:numel(y_pos_side)
            yi = y_pos_side(i);
            if abs(yi) < 1e-12
                if numel(y_list) >= N_cells, break; end
                y_list(end+1)=0;   g_list(end+1)=g; %#ok<AGROW>
            else
                if numel(y_list) >= N_cells, break; end
                y_list(end+1)=+yi; g_list(end+1)=g; %#ok<AGROW>
                if numel(y_list) >= N_cells, break; end
                y_list(end+1)=-yi; g_list(end+1)=g; %#ok<AGROW>
            end
        end
        g = g + 1;
        if g > max_groups
            error('Too many groups; check parameters.');
        end
    end
    [y_pos, perm] = sort(y_list);
    group_idx = g_list(perm);
end

function y_pos = local_compute_parabolic_y_strict(N_cells, z_f, lambda, p_min)
    if mod(N_cells,2) == 1
        idx = -(N_cells-1)/2 : (N_cells-1)/2;
    else
        idx = -(N_cells/2 - 0.5) : (N_cells/2 - 0.5);
    end
    y_pos = sign(idx) .* sqrt(2*abs(idx)*lambda*z_f + (abs(idx)*lambda).^2);
    [y_sorted, perm] = sort(y_pos);
    [~, i_anchor]    = min(abs(idx));
    i_anchor_sorted  = find(perm == i_anchor, 1);
    for i = i_anchor_sorted+1 : numel(y_sorted)
        if y_sorted(i) - y_sorted(i-1) < p_min
            y_sorted(i) = y_sorted(i-1) + p_min;
        end
    end
    for i = i_anchor_sorted-1 : -1 : 1
        if y_sorted(i+1) - y_sorted(i) < p_min
            y_sorted(i) = y_sorted(i+1) - p_min;
        end
    end
    y_pos(perm) = y_sorted;
end

function [y_in_ch, group] = local_y_in_channel(cfg, lambda)
    if strcmpi(cfg.type, 'equal')
        y_in_ch = ((1:cfg.N_cell_y) - (cfg.N_cell_y+1)/2) * cfg.cell_pitch_y;
        group = zeros(size(y_in_ch));
    elseif strcmpi(cfg.type, 'parabolic')
        y_in_ch = local_compute_parabolic_y_strict( ...
            cfg.N_cell_y, cfg.z_focus_geom, lambda, cfg.p_min_y);
        group = zeros(size(y_in_ch));
    elseif strcmpi(cfg.type, 'grouped')
        [y_in_ch, group] = local_compute_grouped_fresnel_y( ...
            cfg.N_cell_y, cfg.z_focus_geom, lambda, cfg.p_min_y, ...
            cfg.dr_frac * lambda);
    else
        error('Unknown cfg.type: %s', cfg.type);
    end
end

function [x_cell, y_cell, ich, group_y] = local_build_cells(cfg, N_ch, ch_pitch, lambda)
    [y_in_ch, group_in_ch] = local_y_in_channel(cfg, lambda);
    x_in_ch = ((1:cfg.N_cell_x) - (cfg.N_cell_x+1)/2) * cfg.cell_pitch_x;
    [ix, iy, ich_grid] = ndgrid(1:cfg.N_cell_x, 1:cfg.N_cell_y, 1:N_ch);
    ich = ich_grid(:);
    ch_centre_x = (ich - (N_ch+1)/2) * ch_pitch;
    x_in_ch_col = x_in_ch(:);
    y_in_ch_col = y_in_ch(:);
    g_in_ch_col = group_in_ch(:);
    dx_in_ch = x_in_ch_col(ix(:));
    dy_in_ch = y_in_ch_col(iy(:));
    group_y  = g_in_ch_col(iy(:));
    x_cell   = ch_centre_x + dx_in_ch;
    y_cell   = dy_in_ch;
end

function P3 = local_compute_field(x_cell, y_cell, cell_tau, X3, Y3, Z3, k, fc)
    [Nx, Ny, Nz] = size(X3);
    x_cell_s   = single(x_cell);
    y_cell_s   = single(y_cell);
    cell_tau_s = single(cell_tau);
    k_s        = single(k);
    omega_s    = single(2*pi*fc);
    P3      = complex(zeros(Nx, Ny, Nz, 'single'));
    N_cells = numel(x_cell);
    fprintf('  Field over %dx%dx%d voxels and %d cells...\n', Nx, Ny, Nz, N_cells);
    tic;
    for i = 1:N_cells
        if mod(i, 1024) == 0
            fprintf('    cell %5d / %d  (%.1f s)\n', i, N_cells, toc);
        end
        r     = sqrt((X3 - x_cell_s(i)).^2 + (Y3 - y_cell_s(i)).^2 + Z3.^2);
        phase = -(k_s*r + omega_s*cell_tau_s(i));
        P3    = P3 + exp(1i*phase) ./ r;
    end
    fprintf('  Done (%.1f s).\n', toc);
end
