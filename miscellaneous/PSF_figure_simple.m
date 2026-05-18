%% PSF COMPARISON FIGURE  --  thesis figure (heatmaps only)
% Reconstructs three single-target water-tank acquisitions with the same
% pipeline as the main reconstruction script, then plots the three
% B-modes side-by-side with:
%   * one shared dB-normalisation across all three (so target-to-target
%     reflectivity differences are visible in the colour scale);
%   * a zoom window around each target peak (so sidelobes are actually
%     visible at print resolution).
%
% Output: ./out/psf_comparison.fig + .pdf + .png

clear; close all; clc;

%% ======================================================================
%  USER INPUT  --  the three targets to compare
%  ======================================================================
%  peak_*_mm   : leave [] for automatic peak detection
%  z_search_mm : optional [z_min z_max] (mm) to constrain the auto search
%                (useful if a PMUT-surface reverb is brighter than the
%                actual target)

% ---------------- Single Targets (different materials) -------------------
% targets(1).name        = 'Nylon String';
% targets(1).file        = 'C:\Users\cioca100\Desktop\TX water tank\PulseEcho_Nylon\Focus 35mm (actual depth ~33p8)\NylonString_Rcv_35.mat';
% targets(1).peak_x_mm   = [];
% targets(1).peak_z_mm   = [];
% targets(1).z_search_mm = [];     % e.g. [25 45] to look only in this band
% 
% targets(2).name        = 'Al rod';
% targets(2).file        = 'C:\Users\cioca100\Desktop\TX water tank\PulseEcho_SmallRod\SmallRod_Rcv.mat';
% targets(2).peak_x_mm   = [];
% targets(2).peak_z_mm   = [];
% targets(2).z_search_mm = [];
% 
% targets(3).name        = 'Tin rod';
% targets(3).file        = 'C:\Users\cioca100\Desktop\TX water tank\PulseEcho_BigRod\BigRod_Rcv.mat';
% targets(3).peak_x_mm   = [];
% targets(3).peak_z_mm   = [];
% targets(3).z_search_mm = [];

% --------------------------- PHANTOM IMAGING -----------------------------
% targets(1).file        = "C:\Users\cioca100\Desktop\TX water tank\Phantom_Imaging_ChipMUT\9p5Sc\Phantom_NF_1_Rcv.mat";
targets(1).name        = 'a) 9.5% Sc - 30.8V';
targets(1).file        = "C:\Users\cioca100\Desktop\TX water tank\Phantom_Imaging_ChipMUT\9p5Sc\5p3VA\Phantom_NF_1_9p5Sc_5p3VA_Rcv.mat";
targets(1).peak_x_mm   = [];
targets(1).peak_z_mm   = [];
targets(1).z_search_mm = [12 12.5];     % e.g. [25 45] to look only in this band

targets(2).name        = 'b) 9.5% Sc - 90.2V';
targets(2).file        = "C:\Users\cioca100\Desktop\TX water tank\Phantom_Imaging_ChipMUT\9p5Sc\Phantom_NF_2_Rcv.mat";
targets(2).peak_x_mm   = [];
targets(2).peak_z_mm   = [];
targets(2).z_search_mm = [];

targets(3).name        = 'c) 30 % Sc - 30.8V';
targets(3).file        = "C:\Users\cioca100\Desktop\TX water tank\Phantom_Imaging_ChipMUT\30Sc\5p3VA\Phantom_NF_1_30Sc_Rcv.mat";
targets(3).fc          = 9e6; % override center frequency [MHz]
targets(3).peak_x_mm   = [];
targets(3).peak_z_mm   = [];
targets(3).z_search_mm = [];

% ------------------------- Res targets -----------------------------------
% targets(1).name        = 'Axial Resolution';
% targets(1).file        = "C:\Users\cioca100\Desktop\WaterTank_Imaging\flat_axial_res_test\Axial_PlaneCustom_2Targets_Rcv.mat";
% targets(1).peak_x_mm   = [];
% targets(1).peak_z_mm   = [];
% targets(1).z_search_mm = []; 
% 
% targets(2).name        = 'Lateral Resolution';
% targets(2).file        = "C:\Users\cioca100\Desktop\WaterTank_Imaging\flat_axial_res_test\Lateral_PlaneCustom_alt_Rcv.mat";
% targets(2).peak_x_mm   = [];
% targets(2).peak_z_mm   = [];
% targets(2).fc          = 10.5e6;
% targets(2).z_search_mm = [10 50];

%% ======================================================================
%  USER INPUT  --  reconstruction parameters (mirror your main script)
%  ======================================================================
params.fs              = 42e6;
params.fc              = 10.5e6;
params.c               = 1497;
params.start_depth_wl  = 18;         % 18 works % match Receive(1).startDepth
params.pitch           = 0.075e-3;
params.element_width   = 0.051e-3;
params.f_number        = 4.0;
params.cf_weight       = 0.75;
params.num_frames_avg  = 30;
params.x_axis          = linspace(-10e-3, 10e-3,  256);
params.z_axis          = linspace( 8e-3, 80e-3, 1024);

%% ======================================================================
%  USER INPUT  --  display parameters
%  ======================================================================
zoom_lat_mm    = 3;             % +/- mm laterally around target peak
zoom_ax_mm     = 3;             % +/- mm axially
heatmap_dB     = [-50 0];       % heatmap colour-axis limits (dB)
cmap           = parula(256);   % colormap (try hot/gray/parula/turbo)
fig_size_cm    = [20 10];       % figure size [W H] in centimetres
out_dir        = 'out';

%% ======================================================================
%  RECONSTRUCT EACH TARGET
%  ======================================================================
N    = numel(targets);
x_mm = params.x_axis * 1000;
z_mm = params.z_axis * 1000;

for k = 1:N
    fprintf('\n=== Reconstructing %s ===\n', targets(k).name);

    % --- per-target fc override --------------------------------------
    p_k = params;
    if ~isempty(targets(k).fc)
        p_k.fc = targets(k).fc;
        fprintf('   Using per-target fc = %.1f MHz\n', p_k.fc/1e6);
    end

    targets(k).BMode = reconstruct_psf(targets(k).file, params);

    % --- peak detection (auto unless overridden) -----------------------
    if isempty(targets(k).peak_x_mm) || isempty(targets(k).peak_z_mm)
        if isempty(targets(k).z_search_mm)
            z_mask = true(size(z_mm));
        else
            zr = targets(k).z_search_mm;
            z_mask = (z_mm >= zr(1)) & (z_mm <= zr(2));
        end
        Bs = targets(k).BMode;
        Bs(~z_mask, :) = -inf;
        [~, ind] = max(Bs(:));
        [iz, ix] = ind2sub(size(Bs), ind);
        if isempty(targets(k).peak_x_mm), targets(k).peak_x_mm = x_mm(ix); end
        if isempty(targets(k).peak_z_mm), targets(k).peak_z_mm = z_mm(iz); end
    end

    % linear-amplitude peak value (used for the shared 0-dB reference)
    [~, ix_pk] = min(abs(x_mm - targets(k).peak_x_mm));
    [~, iz_pk] = min(abs(z_mm - targets(k).peak_z_mm));
    targets(k).peak_lin = targets(k).BMode(iz_pk, ix_pk);

    fprintf('   peak: x = %.2f mm, z = %.2f mm  (linear amp = %.3g)\n', ...
        targets(k).peak_x_mm, targets(k).peak_z_mm, targets(k).peak_lin);
end

% Global maximum across all three targets -> shared 0-dB reference
global_max = max([targets.peak_lin]);
fprintf('\nGlobal max (shared 0-dB reference) = %.3g\n', global_max);
for k = 1:N
    fprintf('   %-15s  peak = %+5.1f dB rel. global max\n', ...
        targets(k).name, 20*log10(targets(k).peak_lin/global_max));
end

%% ======================================================================
%  BUILD FIGURE
%  ======================================================================
fig = figure('Color','w','Units','centimeters', ...
             'Position',[2 2 fig_size_cm], 'Name','PSF comparison');

% one tile per heatmap
tl = tiledlayout(fig, 1, N, ...
                 'TileSpacing','compact', 'Padding','compact');

ax_heat = gobjects(1, N);

for k = 1:N
    B_dB = 20*log10(targets(k).BMode / global_max + 1e-12);

    z_pk  = targets(k).peak_z_mm;
    x_pk  = targets(k).peak_x_mm;
    z_lim = [z_pk - zoom_ax_mm,  z_pk + zoom_ax_mm];
    x_lim = [x_pk - zoom_lat_mm, x_pk + zoom_lat_mm];

    ax_h = nexttile(tl);
    % imagesc(ax_h, x_mm, z_mm, B_dB);
    imagesc(ax_h, x_mm, z_mm, targets(k).BMode);
    set(ax_h, 'YDir','reverse');
    axis(ax_h, 'image');
    xlim(ax_h, x_lim);  ylim(ax_h, z_lim);
    % clim(ax_h, heatmap_dB);
    colormap(ax_h, cmap);
    title(ax_h, targets(k).name, 'FontWeight','bold');
    xlabel(ax_h, 'Lateral (mm)');
    if k == 1
        ylabel(ax_h, 'Axial (mm)');
    else
        ax_h.YTickLabel = [];
    end
    set(ax_h, 'Box','on', 'Layer','top');
    ax_heat(k) = ax_h;
end

% --- shared colour bar ---------------------------------------------
cb = colorbar(ax_heat(end));
cb.Layout.Tile      = 'east';
cb.Label.String     = 'Normalised intensity (dB)';
cb.Label.FontSize   = 10;

%% ======================================================================
%  SAVE
%  ======================================================================
if ~exist(out_dir,'dir'), mkdir(out_dir); end
savefig(fig, fullfile(out_dir, 'psf_comparison.fig'));
exportgraphics(fig, fullfile(out_dir, 'psf_comparison.pdf'), ...
               'ContentType','vector');
exportgraphics(fig, fullfile(out_dir, 'psf_comparison.png'), ...
               'Resolution', 300);
fprintf('\nSaved to %s/psf_comparison.{fig,pdf,png}\n', out_dir);


%% ======================================================================
%  LOCAL FUNCTION  --  reconstruction (encapsulates your main script)
%  ======================================================================
function B_Mode_Lin = reconstruct_psf(filename, p)
% Pipeline: load -> channel reorder -> matched filter -> IQ demod ->
% DAS+CF beamforming on (p.x_axis, p.z_axis). Returns linear-amplitude
% B-mode (size: numel(p.z_axis) x numel(p.x_axis)).

% --- load ----------------------------------------------------------
L = load(filename);
if     isfield(L,'RcvData'),  rawCell = L.RcvData;
elseif isfield(L,'MyMatrix'), rawCell = L.MyMatrix;
else
    f = fieldnames(L);  rawCell = L.(f{1});
end
if iscell(rawCell), RF = rawCell{1}; else, RF = rawCell; end

[~, ~, frames] = size(RF);
Frame_unavg = RF(:, 65:128, :);
nfa   = min(p.num_frames_avg, frames);
Frame = mean(Frame_unavg(:, :, 1:nfa), 3);

% --- channel reorder (interleaved -> physical sequential) ---------
reorder_idx = zeros(1, 64);
for k = 1:32
    reorder_idx(2*k - 1) = 32 + k;
    reorder_idx(2*k)     = k;
end
Frame = Frame(:, reorder_idx);

% --- matched filter -----------------------------------------------
fs = p.fs;  fc = p.fc;
T_pulse   = 2 / (2 * fc);
t_pulse   = 0 : 1/fs : T_pulse;
duty      = 0.67;
ref_pulse = double(mod(t_pulse * fc, 1) < duty) .* sin(2*pi*fc*t_pulse);
ref_pulse = ref_pulse / norm(ref_pulse);
mf_kernel = fliplr(conj(ref_pulse));
Frame_MF  = zeros(size(Frame));
for ch = 1:size(Frame, 2)
    Frame_MF(:, ch) = conv(Frame(:, ch), mf_kernel, 'same');
end
Frame = Frame_MF;

% --- IQ demodulation ----------------------------------------------
[nsamples, nchannels] = size(Frame);
t_vec = (0 : nsamples-1).' / fs;
demod = exp(-1i * 2 * pi * fc * t_vec);
IQ    = Frame .* demod;
[b_lp, a_lp] = butter(3, (fc * 0.8) / (fs/2), 'low');
IQ    = filtfilt(b_lp, a_lp, double(IQ));

% --- DAS + CF beamforming -----------------------------------------
c        = p.c;
lambda   = c / fc;
t_start  = p.start_depth_wl / fc;
pitch    = p.pitch;
elem_w   = p.element_width;
f_num    = p.f_number;

probe_x       = ((0:nchannels-1) - (nchannels-1)/2) * pitch;
half_aperture = (nchannels * pitch) / 2;

[X_pix, Z_pix] = meshgrid(p.x_axis, p.z_axis);
dist_tx = Z_pix;     % plane-wave TX

Sum_DAS = zeros(size(X_pix));
Sum_CF  = zeros(size(X_pix));
Sum_E   = zeros(size(X_pix));
Cnt     = zeros(size(X_pix));

for i = 1:nchannels
    dx       = X_pix - probe_x(i);
    dist_rx  = sqrt(dx.^2 + Z_pix.^2);
    tau      = (dist_tx + dist_rx) / c;
    idx_ex   = (tau - t_start) * fs + 1;

    max_radius = Z_pix / (2 * f_num);
    lat_dist   = abs(dx);
    valid      = (lat_dist <= max_radius) & ...
                 (idx_ex >= 1) & (idx_ex <= nsamples);

    iq_d         = zeros(size(X_pix));
    iq_d(valid)  = interp1(1:nsamples, IQ(:, i), idx_ex(valid), 'linear');

    phase_rot    = exp(1i * 2 * pi * fc * (tau(valid) - t_start));
    iq_a         = zeros(size(X_pix));
    iq_a(valid)  = iq_d(valid) .* phase_rot;

    weights      = zeros(size(X_pix));
    denom_r      = max(max_radius, half_aperture);
    han_w        = 0.5 * (1 + cos(pi * lat_dist(valid) ./ denom_r(valid)));
    sin_th       = dx(valid) ./ dist_rx(valid);
    direct       = sinc((elem_w * sin_th) / lambda);
    weights(valid) = han_w .* direct;

    Sum_DAS = Sum_DAS + iq_a .* weights;
    Sum_CF  = Sum_CF  + iq_a;
    Sum_E   = Sum_E   + abs(iq_a).^2;
    Cnt     = Cnt     + double(valid);
end

Num = abs(Sum_CF).^2;
Den = Cnt .* Sum_E;
Den(Den < eps) = eps;
CF_R = Num ./ Den;
CF_S = imgaussfilt(CF_R, 1.0);

IQ_Final   = Sum_DAS .* (CF_S .^ p.cf_weight);
B_Mode_Lin = abs(IQ_Final);
end






















%==========================================================================

% 
% 
% 
% %% PSF COMPARISON FIGURE  --  thesis figure (heatmaps only)
% % Reconstructs two single-target water-tank acquisitions with the same
% % pipeline as the main reconstruction script, then plots the two
% % B-modes side-by-side.
% %
% % KEY CHANGES vs. previous version:
% %   1. Matched filter now uses a TUKEY-WINDOWED reference pulse
% %      (mismatched filter) to suppress axial time-sidelobes.
% %      alpha = 0  -> rectangular (original behaviour, ~-13 dB sidelobes)
% %      alpha = 0.5-> half-Tukey  (good compromise,   ~-27 dB sidelobes)
% %      alpha = 1  -> full Hanning (maximum suppression,~-42 dB sidelobes,
% %                                  slight main-lobe broadening)
% %   2. T_pulse length corrected to 8 half-cycles (was 4 in the old PSF
% %      script, inconsistent with the main reconstruction script).
% %   3. params.mf_window_alpha added as a top-level tunable parameter.
% %
% % Output: ./out/psf_comparison.fig + .pdf + .png
% 
% clear; close all; clc;
% 
% %% ======================================================================
% %  USER INPUT  --  the targets to compare
% %  ======================================================================
% targets(1).name        = 'Axial Resolution';
% targets(1).file        = "C:\Users\cioca100\Desktop\WaterTank_Imaging\flat_axial_res_test\Axial_PlaneCustom_2Targets_Rcv.mat";
% targets(1).peak_x_mm   = [];
% targets(1).peak_z_mm   = [];
% targets(1).z_search_mm = [];
% 
% targets(2).name        = 'Lateral Resolution';
% targets(2).file        = "C:\Users\cioca100\Desktop\WaterTank_Imaging\flat_axial_res_test\Lateral_PlaneCustom_alt_Rcv.mat";
% targets(2).peak_x_mm   = [];
% targets(2).peak_z_mm   = [];
% targets(2).fc          = 10.5e6;
% targets(2).z_search_mm = [10 50];
% 
% %% ======================================================================
% %  USER INPUT  --  reconstruction parameters
% %  ======================================================================
% params.fs              = 42e6;
% params.fc              = 10.5e6;
% params.c               = 1497;
% params.start_depth_wl  = 18;
% params.pitch           = 0.075e-3;
% params.element_width   = 0.051e-3;
% params.f_number        = 4.0;
% params.cf_weight       = 0.75;
% params.num_frames_avg  = 30;
% params.x_axis          = linspace(-10e-3, 10e-3,  256);
% params.z_axis          = linspace( 8e-3,  80e-3, 1024);
% 
% % -----------------------------------------------------------------------
% % <<< CHANGED: Tukey-window parameter for the matched filter
% %   0   = rectangular envelope (original, highest sidelobes ~-13 dB)
% %   0.5 = half-Tukey            (recommended, ~-27 dB sidelobes)
% %   1.0 = full Hanning          (lowest sidelobes ~-42 dB, slight
% %                                broadening of main lobe)
% % -----------------------------------------------------------------------
% params.mf_window_alpha = 0;   % <<< CHANGED (was implicitly 0 before)
% 
% %% ======================================================================
% %  USER INPUT  --  display parameters
% %  ======================================================================
% zoom_lat_mm = 3;
% zoom_ax_mm  = 3;
% heatmap_dB  = [-50 0];
% cmap        = parula(256);
% fig_size_cm = [20 10];
% out_dir     = 'out';
% 
% %% ======================================================================
% %  RECONSTRUCT EACH TARGET
% %  ======================================================================
% N    = numel(targets);
% x_mm = params.x_axis * 1000;
% z_mm = params.z_axis * 1000;
% 
% for k = 1:N
%     fprintf('\n=== Reconstructing %s ===\n', targets(k).name);
% 
%     p_k = params;
%     if isfield(targets(k),'fc') && ~isempty(targets(k).fc)
%         p_k.fc = targets(k).fc;
%         fprintf('   Using per-target fc = %.1f MHz\n', p_k.fc/1e6);
%     end
% 
%     targets(k).BMode = reconstruct_psf(targets(k).file, p_k);
% 
%     % --- auto peak detection (unless manually overridden) --------------
%     if isempty(targets(k).peak_x_mm) || isempty(targets(k).peak_z_mm)
%         if isempty(targets(k).z_search_mm)
%             z_mask = true(size(z_mm));
%         else
%             zr     = targets(k).z_search_mm;
%             z_mask = (z_mm >= zr(1)) & (z_mm <= zr(2));
%         end
%         Bs = targets(k).BMode;
%         Bs(~z_mask, :) = -inf;
%         [~, ind] = max(Bs(:));
%         [iz, ix] = ind2sub(size(Bs), ind);
%         if isempty(targets(k).peak_x_mm), targets(k).peak_x_mm = x_mm(ix); end
%         if isempty(targets(k).peak_z_mm), targets(k).peak_z_mm = z_mm(iz); end
%     end
% 
%     [~, ix_pk] = min(abs(x_mm - targets(k).peak_x_mm));
%     [~, iz_pk] = min(abs(z_mm - targets(k).peak_z_mm));
%     targets(k).peak_lin = targets(k).BMode(iz_pk, ix_pk);
% 
%     fprintf('   peak: x = %.2f mm, z = %.2f mm  (linear amp = %.3g)\n', ...
%         targets(k).peak_x_mm, targets(k).peak_z_mm, targets(k).peak_lin);
% end
% 
% global_max = max([targets.peak_lin]);
% fprintf('\nGlobal max (shared 0-dB reference) = %.3g\n', global_max);
% for k = 1:N
%     fprintf('   %-20s  peak = %+5.1f dB rel. global max\n', ...
%         targets(k).name, 20*log10(targets(k).peak_lin/global_max));
% end
% 
% %% ======================================================================
% %  BUILD FIGURE
% %  ======================================================================
% fig = figure('Color','w','Units','centimeters', ...
%              'Position',[2 2 fig_size_cm], 'Name','PSF comparison');
% 
% tl       = tiledlayout(fig, 1, N, 'TileSpacing','compact','Padding','compact');
% ax_heat  = gobjects(1, N);
% 
% for k = 1:N
%     B_dB = 20*log10(targets(k).BMode / global_max + 1e-12);
% 
%     z_pk  = targets(k).peak_z_mm;
%     x_pk  = targets(k).peak_x_mm;
%     z_lim = [z_pk - zoom_ax_mm,  z_pk + zoom_ax_mm];
%     x_lim = [x_pk - zoom_lat_mm, x_pk + zoom_lat_mm];
% 
%     ax_h = nexttile(tl);
%     imagesc(ax_h, x_mm, z_mm, B_dB);
%     set(ax_h, 'YDir','reverse');
%     axis(ax_h, 'image');
%     xlim(ax_h, x_lim);  ylim(ax_h, z_lim);
%     clim(ax_h, heatmap_dB);
%     colormap(ax_h, cmap);
%     title(ax_h, targets(k).name, 'FontWeight','bold');
%     xlabel(ax_h, 'Lateral (mm)');
%     if k == 1
%         ylabel(ax_h, 'Axial (mm)');
%     else
%         ax_h.YTickLabel = [];
%     end
%     set(ax_h, 'Box','on', 'Layer','top');
%     ax_heat(k) = ax_h;
% end
% 
% cb = colorbar(ax_heat(end));
% cb.Layout.Tile     = 'east';
% cb.Label.String    = 'Normalised intensity (dB)';
% cb.Label.FontSize  = 10;
% 
% % ---- optional: add alpha annotation to figure title -------------------
% sgtitle(fig, sprintf('PSF Comparison  |  MF Tukey \\alpha = %.1f', ...
%     params.mf_window_alpha), 'FontSize', 12);  % <<< CHANGED (new)
% 
% %% ======================================================================
% %  SAVE
% %  ======================================================================
% if ~exist(out_dir,'dir'), mkdir(out_dir); end
% savefig(fig, fullfile(out_dir, 'psf_comparison.fig'));
% exportgraphics(fig, fullfile(out_dir, 'psf_comparison.pdf'), ...
%                'ContentType','vector');
% exportgraphics(fig, fullfile(out_dir, 'psf_comparison.png'), ...
%                'Resolution', 300);
% fprintf('\nSaved to %s/psf_comparison.{fig,pdf,png}\n', out_dir);
% 
% 
% %% ======================================================================
% %  LOCAL FUNCTION  --  reconstruction
% %  ======================================================================
% function B_Mode_Lin = reconstruct_psf(filename, p)
% 
% % --- load --------------------------------------------------------------
% L = load(filename);
% if     isfield(L,'RcvData'),  rawCell = L.RcvData;
% elseif isfield(L,'MyMatrix'), rawCell = L.MyMatrix;
% else,  f = fieldnames(L);     rawCell = L.(f{1});
% end
% if iscell(rawCell), RF = rawCell{1}; else, RF = rawCell; end
% 
% [~, ~, frames] = size(RF);
% Frame_unavg = RF(:, 65:128, :);
% nfa   = min(p.num_frames_avg, frames);
% Frame = mean(Frame_unavg(:, :, 1:nfa), 3);
% 
% % --- channel reorder ---------------------------------------------------
% reorder_idx = zeros(1, 64);
% for k = 1:32
%     reorder_idx(2*k - 1) = 32 + k;
%     reorder_idx(2*k)     = k;
% end
% Frame = Frame(:, reorder_idx);
% 
% % --- matched filter (Tukey-windowed) -----------------------------------
% fs = p.fs;  fc = p.fc;
% 
% % <<< CHANGED: 8 half-cycles (was 4 in old PSF script — now matches main script)
% T_pulse   = 2 / (2 * fc);
% t_pulse   = 0 : 1/fs : T_pulse;
% duty      = 0.67;
% 
% % Rectangular-envelope reference pulse (models the actual TX waveform)
% ref_pulse = double(mod(t_pulse * fc, 1) < duty) .* sin(2*pi*fc*t_pulse);
% 
% % <<< CHANGED: apply Tukey window to the MF kernel
% %   This is a "mismatched filter": it no longer perfectly matches the TX
% %   pulse, but the tapered autocorrelation has far lower time-sidelobes.
% %   Cost: ~1-2 dB SNR loss and very slight main-lobe broadening.
% N_pulse  = length(ref_pulse);
% alpha    = p.mf_window_alpha;           % 0=rect, 0.5=half-Tukey, 1=Hanning
% 
% % Build Tukey window of length N_pulse
% %   The tukeywin() function is available in Signal Processing Toolbox.
% %   If you don't have it, the manual fallback below is used instead.
% if license('test','Signal_Toolbox') && exist('tukeywin','file')
%     t_win = tukeywin(N_pulse, alpha)';
% else
%     % Manual Tukey window (identical result, no toolbox required)
%     t_win = ones(1, N_pulse);
%     n_taper = floor(alpha * N_pulse / 2);
%     if n_taper > 0
%         taper_idx = 1 : n_taper;
%         taper     = 0.5 * (1 - cos(pi * (taper_idx - 1) / (n_taper - 1)));
%         t_win(taper_idx)                    = taper;           % leading ramp
%         t_win(N_pulse - n_taper + 1 : end)  = fliplr(taper);  % trailing ramp
%     end
% end
% 
% % Apply window and normalise
% ref_pulse_win = ref_pulse .* t_win;
% ref_pulse_win = ref_pulse_win / norm(ref_pulse_win);  % <<< CHANGED (normalise windowed version)
% 
% mf_kernel = fliplr(conj(ref_pulse_win));              % time-reversed conjugate
% Frame_MF  = zeros(size(Frame));
% for ch = 1:size(Frame, 2)
%     Frame_MF(:, ch) = conv(Frame(:, ch), mf_kernel, 'same');
% end
% Frame = Frame_MF;
% 
% % --- IQ demodulation ---------------------------------------------------
% [nsamples, nchannels] = size(Frame);
% t_vec = (0 : nsamples-1).' / fs;
% demod = exp(-1i * 2 * pi * fc * t_vec);
% IQ    = Frame .* demod;
% [b_lp, a_lp] = butter(3, (fc * 0.8) / (fs/2), 'low');
% IQ    = filtfilt(b_lp, a_lp, double(IQ));
% 
% % --- DAS + CF beamforming ----------------------------------------------
% c        = p.c;
% lambda   = c / fc;
% t_start  = p.start_depth_wl / fc;
% pitch    = p.pitch;
% elem_w   = p.element_width;
% f_num    = p.f_number;
% 
% probe_x       = ((0:nchannels-1) - (nchannels-1)/2) * pitch;
% half_aperture = (nchannels * pitch) / 2;
% 
% [X_pix, Z_pix] = meshgrid(p.x_axis, p.z_axis);
% dist_tx = Z_pix;
% 
% Sum_DAS = zeros(size(X_pix));
% Sum_CF  = zeros(size(X_pix));
% Sum_E   = zeros(size(X_pix));
% Cnt     = zeros(size(X_pix));
% 
% for i = 1:nchannels
%     dx       = X_pix - probe_x(i);
%     dist_rx  = sqrt(dx.^2 + Z_pix.^2);
%     tau      = (dist_tx + dist_rx) / c;
%     idx_ex   = (tau - t_start) * fs + 1;
% 
%     max_radius = Z_pix / (2 * f_num);
%     lat_dist   = abs(dx);
%     valid      = (lat_dist <= max_radius) & ...
%                  (idx_ex >= 1) & (idx_ex <= nsamples);
% 
%     iq_d        = zeros(size(X_pix));
%     iq_d(valid) = interp1(1:nsamples, IQ(:,i), idx_ex(valid), 'linear');
% 
%     phase_rot   = exp(1i * 2 * pi * fc * (tau(valid) - t_start));
%     iq_a        = zeros(size(X_pix));
%     iq_a(valid) = iq_d(valid) .* phase_rot;
% 
%     weights       = zeros(size(X_pix));
%     denom_r       = max(max_radius, half_aperture);
%     han_w         = 0.5 * (1 + cos(pi * lat_dist(valid) ./ denom_r(valid)));
%     sin_th        = dx(valid) ./ dist_rx(valid);
%     direct        = sinc((elem_w * sin_th) / lambda);
%     weights(valid)= han_w .* direct;
% 
%     Sum_DAS = Sum_DAS + iq_a .* weights;
%     Sum_CF  = Sum_CF  + iq_a;
%     Sum_E   = Sum_E   + abs(iq_a).^2;
%     Cnt     = Cnt     + double(valid);
% end
% 
% Num = abs(Sum_CF).^2;
% Den = Cnt .* Sum_E;
% Den(Den < eps) = eps;
% CF_R = Num ./ Den;
% CF_S = imgaussfilt(CF_R, 1.0);
% 
% IQ_Final   = Sum_DAS .* (CF_S .^ p.cf_weight);
% B_Mode_Lin = abs(IQ_Final);
% end