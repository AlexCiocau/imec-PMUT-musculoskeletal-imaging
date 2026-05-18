%% =========================================================================
%  L7-4 RECONSTRUCTION & SNR CALCULATION
%  Adapted from PMUT pipeline — IQ demod + DAS + CF beamforming
%  Key differences from PMUT:
%    - Buffer stitching replaces channel reordering
%    - No interleaving: channels are sequential after stitching
%    - Different fc, fs, pitch, element width, pulse definition
% =========================================================================

%% 1. BUFFER SETUP & STITCHING
clear all; close all; clc;

% -------------------- December 2025 - Measurements ----------------------
% filename = 'MultiPointRFData.mat';
% filename = "C:\Users\cioca100\Desktop\Verasonics Matlab code\L7_4\AxialLateral_Res_RFData.mat";
% filename = 'SensCutoff_axlat_Res.mat';
% filename = 'EchoicChambers';
% ------------------------------------------------------------------------


% ----------------------- March 2026 - Measurements ----------------------
% filename = "C:\Users\cioca100\Desktop\TX water tank\Phantom_Imaging_L7_4\L7_4_NF_Vertical_deep_Rcv.mat";
filename = "C:\Users\cioca100\Desktop\TX water tank\Phantom_Imaging_L7_4\L7_4_ResGroup_deep_Rcv.mat";
% filename = "C:\Users\cioca100\Desktop\TX water tank\Phantom_Imaging_L7_4\L7_4_Anechoic_Rcv.mat";
% ------------------------------------------------------------------------

fprintf('Loading RF data from: %s ...\n', filename);
loadedData = load(filename);

if isfield(loadedData, 'RcvData')
    rawCell = loadedData.RcvData;
elseif isfield(loadedData, 'MyMatrix')
    rawCell = loadedData.MyMatrix;
else
    vars = fieldnames(loadedData);
    rawCell = loadedData.(vars{1});
end

if iscell(rawCell)
    RF_Matrix = rawCell{1};
else
    RF_Matrix = rawCell;
end

[rows, cols, frames] = size(RF_Matrix);
fprintf('Buffer size: %d samples x %d channels x %d frames\n', rows, cols, frames);

% ── Stitch the two acquisition halves ────────────────────────────────────
% The L7-4 acquires in two passes due to Verasonics channel limitations:
%   Pass 1: channels  1–64,  rows 1:4096
%   Pass 2: channels 65–128, rows offset by start index (depth-dependent)
% Both halves are aligned to the same time origin after stitching.

framesToView = 1;

FrameData_1 = RF_Matrix(1:4096,    1:64,  framesToView);
FrameData_2 = RF_Matrix(2690:6785, 65:128, framesToView);  % max depth — tune if needed

CompleteFrame = zeros(4096, 128);
CompleteFrame(1:4096,    1:64)  = FrameData_1;
CompleteFrame(1:4096,   65:128) = FrameData_2;

% No channel reordering needed — L7-4 channels are sequential after stitching
fprintf('Buffer stitched: 128 channels x 4096 samples.\n');

% figure(9); clf;
% imagesc(abs(hilbert(CompleteFrame)))
% % colormap gray; 
% colorbar;
% title('L7-4: Averaged RF data (before beamforming)');
% xlabel('Lateral (mm)'); ylabel('Axial (mm)');
% axis image;
% hold on;
% xline(-array_edge_mm, 'Color','w', 'LineStyle',':', 'LineWidth',2);
% xline( array_edge_mm, 'Color','w', 'LineStyle',':', 'LineWidth',2);
% hold off;

%% 2. MATCHED FILTER
% ── Probe & acquisition parameters — verify against your Verasonics script ──
fs = 20.832e6;    % <── sampling frequency (Hz): typically 4x fc in Verasonics
fc =  5.208e6;    % <── centre frequency (Hz): L7-4 nominal ~5.208 MHz

fprintf('Applying matched filter (fc = %.3f MHz, fs = %.3f MHz)...\n', ...
    fc/1e6, fs/1e6);

% L7-4 transmit pulse: adjust n_cycles and duty to match your Verasonics
% sequence (check numTxCycles and TW(1).Duty in your .m script)
n_cycles = 2;        % <── number of transmit cycles (commonly 2–3 for L7-4)
duty     = 0.67;     % <── duty cycle

T_pulse   = n_cycles / fc;
t_pulse   = 0 : 1/fs : T_pulse;
ref_pulse = double(mod(t_pulse * fc, 1) < duty) .* sin(2*pi*fc*t_pulse);
ref_pulse = ref_pulse / norm(ref_pulse);

mf_kernel        = fliplr(conj(ref_pulse));
CompleteFrame_MF = zeros(size(CompleteFrame));
for ch = 1:size(CompleteFrame, 2)
    CompleteFrame_MF(:, ch) = conv(CompleteFrame(:, ch), mf_kernel, 'same');
end
CompleteFrame = CompleteFrame_MF;
fprintf('Matched filter applied (pulse length: %d samples).\n', length(ref_pulse));

%% 3. IQ DEMODULATION
fprintf('Demodulating to IQ baseband...\n');

[nsamples, nchannels] = size(CompleteFrame);
t_vector = (0:nsamples-1)' / fs;

demod_carrier = exp(-1i * 2 * pi * fc * t_vector);
IQ_Matrix     = CompleteFrame .* demod_carrier;

% Low-pass cutoff at 80% of fc — same as PMUT pipeline
[b_lp, a_lp]  = butter(3, (fc * 0.8) / (fs/2), 'low');
IQ_Matrix     = filtfilt(b_lp, a_lp, double(IQ_Matrix));

%% 4. IQ BASEBAND BEAMFORMING (PLANE WAVE)
fprintf('Beamforming setup...\n');

c      = 1540;        % sound speed in tissue (m/s) — use 1480 for water tank
lambda = c / fc;

% ── Probe geometry — verify against your Verasonics Resource.Parameters ──
pitch         = 0.298e-3;    % L7-4 element pitch (m)
element_width = 0.208e-3;    % L7-4 element width (~0.7 * pitch)
f_number      = 1.5;         % <── tune: L7-4 is a focused probe, F# ~1–2 typical

% startDepth in your Verasonics Receive struct (in wavelengths)
start_depth_wl = 5;          % <── verify from Receive(1).startDepth
t_start        = start_depth_wl / fc;

probe_x       = ((0:nchannels-1) - (nchannels-1)/2) * pitch;
half_aperture = (nchannels * pitch) / 2;

% Imaging grid — adjust depth range to your target depth
z_axis = linspace(1e-3, 95e-3, 4096);
x_axis = linspace(-20e-3, 20e-3, 512);   % L7-4 has wider aperture than PMUT
[X_pix, Z_pix] = meshgrid(x_axis, z_axis);

fprintf('Beamforming (%d channels, %.1f mm aperture)...\n', ...
    nchannels, nchannels*pitch*1000);

dist_tx = Z_pix;   % plane wave transmit

Sum_DAS_Weighted = zeros(size(X_pix));
Sum_CF_Signal    = zeros(size(X_pix));
Sum_CF_Energy    = zeros(size(X_pix));
Count_Active     = zeros(size(X_pix));

for i = 1:nchannels
    dx      = X_pix - probe_x(i);
    dist_rx = sqrt(dx.^2 + Z_pix.^2);
    tau     = (dist_tx + dist_rx) / c;

    idx_exact = (tau - t_start) * fs + 1;

    max_radius = Z_pix / (2 * f_number);
    lat_dist   = abs(dx);
    in_cone    = lat_dist <= max_radius;
    in_time    = (idx_exact >= 1) & (idx_exact <= nsamples);
    valid_mask = in_cone & in_time;

    iq_col            = IQ_Matrix(:, i);
    iq_delayed_smooth = zeros(size(X_pix));
    iq_delayed_smooth(valid_mask) = interp1(1:nsamples, iq_col, ...
        idx_exact(valid_mask), 'linear');

    phase_rotation = exp(1i * 2 * pi * fc * (tau(valid_mask) - t_start));

    iq_aligned = zeros(size(X_pix));
    iq_aligned(valid_mask) = iq_delayed_smooth(valid_mask) .* phase_rotation;

    weights      = zeros(size(X_pix));
    denom_radius = max(max_radius, half_aperture);
    hanning_weight = 0.5 * (1 + cos(pi * lat_dist(valid_mask) ./ ...
        denom_radius(valid_mask)));

    sin_theta   = dx(valid_mask) ./ dist_rx(valid_mask);
    directivity = sinc((element_width * sin_theta) / lambda);

    weights(valid_mask) = hanning_weight .* directivity;

    Sum_DAS_Weighted = Sum_DAS_Weighted + (iq_aligned .* weights);
    Sum_CF_Signal    = Sum_CF_Signal    + iq_aligned;
    Sum_CF_Energy    = Sum_CF_Energy    + (abs(iq_aligned).^2);
    Count_Active     = Count_Active     + double(valid_mask);
end

% Coherence Factor
Numerator   = abs(Sum_CF_Signal).^2;
Denominator = Count_Active .* Sum_CF_Energy;
Denominator(Denominator < eps) = eps;
CF_Raw    = Numerator ./ Denominator;
CF_Smooth = imgaussfilt(CF_Raw, 1.0);

cf_weight = 0.75;
% cf_weight = 0;
IQ_Final  = Sum_DAS_Weighted .* (CF_Smooth .^ cf_weight);

%% 5. DISPLAY
B_Mode_Lin = abs(IQ_Final);
B_Mode_dB  = 20 * log10(B_Mode_Lin / max(B_Mode_Lin(:)) + 1e-12);

array_edge_mm = (nchannels * pitch * 1000) / 2;   % ~19.1 mm for L7-4


% ── Figure 1: Full reconstruction ────────────────────────────────────────
figure(1); clf;
imagesc(x_axis*1000, z_axis*1000, B_Mode_dB);
clim([-40 0]);
colormap parula; 
colorbar;
title('L7-4: Full Reconstruction');
xlabel('Lateral (mm)'); ylabel('Axial (mm)');
axis image;
% hold on;
% xline(-array_edge_mm, 'Color','w', 'LineStyle',':', 'LineWidth',2);
% xline( array_edge_mm, 'Color','w', 'LineStyle',':', 'LineWidth',2);
% hold off;


% -------------------------------------------------------------------------

% %% 5. PREPARE IMAGES FOR COMPARISON
% fprintf('Preparing images for display...\n');
% 
% % 1. Raw RF (Envelope & dB)
% % Map samples to depth (mm) using round-trip time, and channels to lateral (mm)
% z_raw_mm = (0:nsamples-1)' / fs * c / 2 * 1000; 
% x_raw_mm = probe_x * 1000;
% Raw_Env  = abs(hilbert(CompleteFrame));
% Raw_dB   = 20 * log10(Raw_Env / max(Raw_Env(:)) + 1e-12);
% 
% % 2. DAS Only
% DAS_Lin = abs(Sum_DAS_Weighted);
% DAS_dB  = 20 * log10(DAS_Lin / max(DAS_Lin(:)) + 1e-12);
% 
% % 3. DAS + Coherence Factor
% % Apply full CF weight (exponent = 1) for the comparison
% IQ_CF  = Sum_DAS_Weighted .* CF_Smooth; 
% CF_Lin = abs(IQ_CF);
% CF_dB  = 20 * log10(CF_Lin / max(CF_Lin(:)) + 1e-12);

%% 6. EXPORT FOR PYTHON RESOLUTION FIGURE
% ── Select target manually ───────────────────────────────────────────────
% Set these to the depth and lateral position of the wire target you want
% to analyse (in mm). Read them off Figure 1 by hovering the cursor.
target_z_mm =  30.0;   % <── edit: axial depth of target (mm)
target_x_mm =   2.0;   % <── edit: lateral position of target (mm)

[~, iz] = min(abs(z_axis*1000 - target_z_mm));
[~, ix] = min(abs(x_axis*1000 - target_x_mm));

lateral_profile = B_Mode_dB(iz, :);   % horizontal slice at target depth
axial_profile   = B_Mode_dB(:, ix);   % vertical slice at target lateral pos

export = struct();
export.B_Mode_dB        = B_Mode_dB;
export.x_axis_mm        = x_axis * 1000;
export.z_axis_mm        = z_axis * 1000;
export.lateral_profile  = lateral_profile;
export.axial_profile    = axial_profile;
export.target_iz        = iz;
export.target_ix        = ix;
export.target_z_mm      = target_z_mm;
export.target_x_mm      = target_x_mm;
export.fc_MHz           = fc / 1e6;
export.fs_MHz           = fs / 1e6;
export.c_mps            = c;

save('resolution_export.mat', '-struct', 'export');
fprintf('Saved resolution_export.mat  (target at x=%.1f mm, z=%.1f mm)\n', ...
        target_x_mm, target_z_mm);

% %% 6. SIDE-BY-SIDE FIGURE
% dyn_range = [-40 0]; % Shared dynamic range for all plots
% 
% % Create a wide figure for the 3 subplots
% figure('Name', 'Reconstruction Comparison', 'Position', [100, 200, 1200, 500]);
% tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
% 
% % --- Plot 1: Raw RF (Unbeamformed) ---
% nexttile;
% imagesc(x_raw_mm, z_raw_mm, Raw_dB);
% clim(dyn_range); % Note: use caxis(dyn_range) if on MATLAB older than R2022a
% colormap parula; 
% title('Raw RF Data (Envelope)');
% xlabel('Lateral (mm)'); ylabel('Axial (mm)');
% axis image;
% 
% % --- Plot 2: DAS Beamforming ---
% nexttile;
% imagesc(x_axis*1000, z_axis*1000, DAS_dB);
% clim(dyn_range); 
% colormap parula;
% title('Delay-and-Sum (DAS)');
% xlabel('Lateral (mm)'); ylabel('Axial (mm)');
% axis image;
% 
% % --- Plot 3: DAS + Coherence Factor ---
% nexttile;
% imagesc(x_axis*1000, z_axis*1000, CF_dB);
% clim(dyn_range); 
% colormap parula;
% title('DAS + Coherence Factor');
% xlabel('Lateral (mm)'); ylabel('Axial (mm)');
% axis image;
% 
% % Add a single shared colorbar on the right
% cb = colorbar;
% cb.Layout.Tile = 'east';
% cb.Label.String = 'Normalized Amplitude (dB)';
% 
% fprintf('Done.\n');


% %% 6. SIDE-BY-SIDE FIGURE
% dyn_range = [-40 0]; % Shared dynamic range for all plots
% 
% % Define common Field of View (FOV) limits based on your imaging grid
% x_limits = [min(x_axis)*1000, max(x_axis)*1000]; % [-20 20] mm
% z_limits = [min(z_axis)*1000, max(z_axis)*1000]; % [1 95] mm
% 
% % Create a wide figure for the 3 subplots
% figure('Name', 'Reconstruction Comparison', 'Position', [100, 200, 1200, 500]);
% tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
% 
% % --- Plot 1: Raw RF (Unbeamformed) ---
% nexttile;
% imagesc(x_raw_mm, z_raw_mm, Raw_dB);
% clim(dyn_range); 
% colormap parula; % Switched back to parula to match your screenshot
% title('Raw RF Data (Envelope)');
% xlabel('Lateral (mm)'); ylabel('Axial (mm)');
% axis image;
% xlim(x_limits); % Force same lateral FOV
% ylim(z_limits); % Force same axial FOV (crops the deep unbeamformed data)
% 
% % --- Plot 2: DAS Beamforming ---
% nexttile;
% imagesc(x_axis*1000, z_axis*1000, DAS_dB);
% clim(dyn_range); 
% colormap parula;
% title('Delay-and-Sum (DAS)');
% xlabel('Lateral (mm)'); ylabel('Axial (mm)');
% axis image;
% xlim(x_limits);
% ylim(z_limits);
% 
% % --- Plot 3: DAS + Coherence Factor ---
% nexttile;
% imagesc(x_axis*1000, z_axis*1000, CF_dB);
% clim(dyn_range); 
% colormap parula;
% title('DAS + Coherence Factor');
% xlabel('Lateral (mm)'); ylabel('Axial (mm)');
% axis image;
% xlim(x_limits);
% ylim(z_limits);
% 
% % Add a single shared colorbar on the right
% cb = colorbar;
% cb.Layout.Tile = 'east';
% cb.Label.String = 'Normalized Amplitude (dB)';
% 
% fprintf('Done.\n');


% % ── Figure 2: Zoomed view around target ──────────────────────────────────
% % Inspect Figure 1 first, then set target_depth_mm accordingly
% target_depth_mm  = 10;     % <── adjust after inspecting Figure 1
% zoom_axial_range = 5;
% zoom_lat_range   = 5;
% 
% figure(2); clf;
% imagesc(x_axis*1000, z_axis*1000, B_Mode_dB);
% clim([-30 0]);
% colormap gray; colorbar;
% title('L7-4: Zoomed View');
% xlabel('Lateral (mm)'); ylabel('Axial (mm)');
% axis image;
% xlim([-zoom_lat_range, zoom_lat_range]);
% ylim([target_depth_mm - zoom_axial_range, target_depth_mm + zoom_axial_range]);
% hold on;
% xline(-array_edge_mm, 'Color','w', 'LineStyle',':', 'LineWidth',2);
% xline( array_edge_mm, 'Color','w', 'LineStyle',':', 'LineWidth',2);
% hold off;


% figure(1); clf;
% imagesc(x_axis*1000, z_axis*1000, B_Mode_dB);
% clim([-60 0]);     % L7-4 has higher dynamic range than PMUT — widen clim
% colormap gray; colorbar;
% title('L7-4: MF + IQ Demod + Directivity + CF');
% xlabel('Lateral (mm)'); ylabel('Axial (mm)');
% axis image;
% hold on;
% xline(-array_edge_mm, 'Color','w', 'LineStyle',':', 'LineWidth',2);
% xline( array_edge_mm, 'Color','w', 'LineStyle',':', 'LineWidth',2);
% hold off;

% %% 6. SNR CALCULATION
% % Identical procedure to PMUT — allows direct comparison
% fprintf('\n--- SNR Calculation (L7-4) ---\n');
% 
% target_depth_mm  = 10;     % <── set to your wire target depth
% zoom_axial_range = 5;
% zoom_lat_range   = 5;      % wider than PMUT due to larger aperture
% 
% x_zoom = [-zoom_lat_range,  zoom_lat_range];
% z_zoom = [target_depth_mm - zoom_axial_range, ...
%           target_depth_mm + zoom_axial_range];
% 
% figure(3); clf;
% imagesc(x_axis*1000, z_axis*1000, B_Mode_dB);
% clim([-40 0]);
% colormap gray; colorbar;
% title('SNR ROI Selection  |  Draw ROI 1: Signal (wire target)');
% xlabel('Lateral (mm)'); ylabel('Axial (mm)');
% axis image;
% xlim(x_zoom); ylim(z_zoom);
% hold on;
% 
% fprintf('Draw ROI 1: Signal region (wire). Double-click to confirm.\n');
% h_roi1 = drawrectangle('Color','r', 'Label','Signal', 'LabelAlpha',0);
% wait(h_roi1);
% 
% pos1     = h_roi1.Position;
% x1_range = [pos1(1), pos1(1)+pos1(3)];
% z1_range = [pos1(2), pos1(2)+pos1(4)];
% 
% x_mm = x_axis * 1000;
% z_mm = z_axis * 1000;
% x1_idx = find(x_mm >= x1_range(1) & x_mm <= x1_range(2));
% z1_idx = find(z_mm >= z1_range(1) & z_mm <= z1_range(2));
% 
% roi1_pixels = B_Mode_Lin(z1_idx, x1_idx);
% roi1_pixels = roi1_pixels(:);
% 
% title('SNR ROI Selection  |  Draw ROI 2: Background (same depth)');
% xlim(x_zoom); ylim(z_zoom);
% 
% fprintf('Draw ROI 2: Background region (same depth, no targets). Double-click to confirm.\n');
% h_roi2 = drawrectangle('Color','b', 'Label','Background', 'LabelAlpha',0);
% wait(h_roi2);
% 
% pos2     = h_roi2.Position;
% x2_range = [pos2(1), pos2(1)+pos2(3)];
% z2_range = [pos2(2), pos2(2)+pos2(4)];
% 
% x2_idx = find(x_mm >= x2_range(1) & x_mm <= x2_range(2));
% z2_idx = find(z_mm >= z2_range(1) & z_mm <= z2_range(2));
% 
% roi2_pixels = B_Mode_Lin(z2_idx, x2_idx);
% roi2_pixels = roi2_pixels(:);
% 
% % ── Depth alignment check ────────────────────────────────────────────────
% z1_centre  = mean(z1_range);
% z2_centre  = mean(z2_range);
% depth_diff = abs(z1_centre - z2_centre);
% 
% if depth_diff > 3.0
%     warning('ROI centres differ by %.1f mm in depth — attenuation bias risk.', depth_diff);
% else
%     fprintf('Depth alignment OK: centres differ by %.2f mm.\n', depth_diff);
% end
% 
% n1 = numel(roi1_pixels);
% n2 = numel(roi2_pixels);
% fprintf('ROI 1 size: %d pixels | ROI 2 size: %d pixels\n', n1, n2);
% if min(n1,n2) < 20
%     warning('ROI too small — increase for reliable statistics.');
% end
% 
% % ── Compute SNR ──────────────────────────────────────────────────────────
% mu_h     = mean(roi1_pixels);
% mu_s     = mean(roi2_pixels);
% sigma2_h = var(roi1_pixels);
% sigma2_s = var(roi2_pixels);
% 
% SNR_linear = abs(mu_h - mu_s) / sqrt((sigma2_h + sigma2_s) / 2);
% SNR_dB     = 20 * log10(SNR_linear);
% 
% fprintf('\n--- Results (L7-4) ---\n');
% fprintf('  mu_h    (signal mean):      %.4e\n', mu_h);
% fprintf('  mu_s    (background mean):  %.4e\n', mu_s);
% fprintf('  sigma_h (signal std):       %.4e\n', sqrt(sigma2_h));
% fprintf('  sigma_s (background std):   %.4e\n', sqrt(sigma2_s));
% fprintf('  SNR (linear):               %.4f\n', SNR_linear);
% fprintf('  SNR (dB):                   %.2f dB\n', SNR_dB);
% 
% title(sprintf('L7-4  |  SNR = %.2f dB  |  Red = Signal, Blue = Background', SNR_dB));
% hold off;
% 
% SNR_result_L74.mu_h       = mu_h;
% SNR_result_L74.mu_s       = mu_s;
% SNR_result_L74.sigma2_h   = sigma2_h;
% SNR_result_L74.sigma2_s   = sigma2_s;
% SNR_result_L74.SNR_linear = SNR_linear;
% SNR_result_L74.SNR_dB     = SNR_dB;
% SNR_result_L74.roi1_pos   = pos1;
% SNR_result_L74.roi2_pos   = pos2;
% % Save with: save('SNR_L74_10mm.mat', 'SNR_result_L74');
% 
% 
% %% 7. ROI HISTOGRAM SANITY CHECK
% figure(4); clf;
% 
% % Use the same bin edges for both so they're directly comparable
% all_pixels = [roi1_pixels; roi2_pixels];
% n_bins     = 40;
% edges      = linspace(min(all_pixels), max(all_pixels), n_bins+1);
% 
% histogram(roi2_pixels, edges, ...
%     'Normalization','probability', ...
%     'FaceColor',[0.2 0.4 0.8], 'FaceAlpha',0.55, ...
%     'DisplayName','Background (ROI 2)');
% hold on;
% histogram(roi1_pixels, edges, ...
%     'Normalization','probability', ...
%     'FaceColor',[0.85 0.2 0.2], 'FaceAlpha',0.55, ...
%     'DisplayName','Signal (ROI 1)');
% 
% % Mark the means
% xline(mu_s, '--', 'Color',[0.1 0.3 0.7], 'LineWidth',1.8, ...
%     'Label',sprintf('\\mu_s = %.2e', mu_s), 'LabelVerticalAlignment','top');
% xline(mu_h, '--', 'Color',[0.7 0.1 0.1], 'LineWidth',1.8, ...
%     'Label',sprintf('\\mu_h = %.2e', mu_h), 'LabelVerticalAlignment','middle');
% 
% xlabel('Pixel intensity (linear, |IQ|)');
% ylabel('Probability');
% title(sprintf('ROI histograms  |  SNR = %.2f dB', SNR_dB));
% legend('Location','northeast');
% grid on;
% hold off;
% 
% % Quick numeric diagnostic
% fprintf('\n--- Histogram diagnostic ---\n');
% fprintf('  Signal     : mean = %.3e, std = %.3e, median = %.3e\n', ...
%     mu_h, sqrt(sigma2_h), median(roi1_pixels));
% fprintf('  Background : mean = %.3e, std = %.3e, median = %.3e\n', ...
%     mu_s, sqrt(sigma2_s), median(roi2_pixels));
% fprintf('  Mean/median ratio (signal):     %.2f  (>> 1 suggests bimodal)\n', ...
%     mu_h / median(roi1_pixels));
% fprintf('  Contrast (mu_h / mu_s):         %.2f  (%.1f dB)\n', ...
%     mu_h/mu_s, 20*log10(mu_h/mu_s));


% %% 8. CNR CALCULATION (hypoechoic target)
% % CNR = |mu_s - mu_c| / sqrt((sigma2_s + sigma2_c)/2)
% %   mu_c, sigma2_c  -> ROI inside the hypoechoic (dark) cylinder
% %   mu_s, sigma2_s  -> ROI in surrounding speckle at the SAME depth
% %
% % Uses the DAS-only envelope (IQ_DAS) — not the CF-weighted image —
% % so Rayleigh speckle statistics are preserved.
% 
% fprintf('\n--- CNR Calculation (L7-4) ---\n');
% 
% % Build DAS-only B-mode if it doesn't already exist from section 7
% if ~exist('B_DAS_lin','var')
%     B_DAS_lin = abs(Sum_DAS_Weighted);
% end
% B_DAS_dB = 20 * log10(B_DAS_lin / max(B_DAS_lin(:)) + 1e-12);
% 
% % ── ROI selection on the DAS image ───────────────────────────────────────
% target_depth_mm_cnr = 30;    % <── adjust to the depth of your hypoechoic cylinder
% zoom_axial_range    = 6;
% zoom_lat_range      = 10;
% 
% x_zoom_cnr = [-zoom_lat_range,  zoom_lat_range];
% z_zoom_cnr = [target_depth_mm_cnr - zoom_axial_range, ...
%               target_depth_mm_cnr + zoom_axial_range];
% 
% figure(5); clf;
% imagesc(x_axis*1000, z_axis*1000, B_DAS_dB);
% clim([-50 0]);
% colormap gray; colorbar;
% title('CNR ROI Selection  |  Draw ROI C: Hypoechoic (dark) cylinder');
% xlabel('Lateral (mm)'); ylabel('Axial (mm)');
% axis image;
% xlim(x_zoom_cnr); ylim(z_zoom_cnr);
% hold on;
% 
% fprintf('Draw ROI C: Hypoechoic cylinder interior. Double-click to confirm.\n');
% h_roiC = drawrectangle('Color',[0.2 0.8 0.2], 'Label','Hypo', 'LabelAlpha',0);
% wait(h_roiC);
% 
% posC     = h_roiC.Position;
% xC_range = [posC(1), posC(1)+posC(3)];
% zC_range = [posC(2), posC(2)+posC(4)];
% 
% x_mm = x_axis * 1000;
% z_mm = z_axis * 1000;
% xC_idx = find(x_mm >= xC_range(1) & x_mm <= xC_range(2));
% zC_idx = find(z_mm >= zC_range(1) & z_mm <= zC_range(2));
% 
% roiC_pixels = B_DAS_lin(zC_idx, xC_idx);
% roiC_pixels = roiC_pixels(:);
% 
% title('CNR ROI Selection  |  Draw ROI S: Background speckle (same depth)');
% xlim(x_zoom_cnr); ylim(z_zoom_cnr);
% 
% fprintf('Draw ROI S: Background speckle at the same depth. Double-click to confirm.\n');
% h_roiS = drawrectangle('Color',[0 0.4 1], 'Label','Background', 'LabelAlpha',0);
% wait(h_roiS);
% 
% posS     = h_roiS.Position;
% xS_range = [posS(1), posS(1)+posS(3)];
% zS_range = [posS(2), posS(2)+posS(4)];
% 
% xS_idx = find(x_mm >= xS_range(1) & x_mm <= xS_range(2));
% zS_idx = find(z_mm >= zS_range(1) & z_mm <= zS_range(2));
% 
% roiS_pixels = B_DAS_lin(zS_idx, xS_idx);
% roiS_pixels = roiS_pixels(:);
% 
% % ── Depth alignment check ────────────────────────────────────────────────
% zC_centre  = mean(zC_range);
% zS_centre  = mean(zS_range);
% depth_diff = abs(zC_centre - zS_centre);
% 
% if depth_diff > 3.0
%     warning('ROI centres differ by %.1f mm in depth — attenuation bias risk.', depth_diff);
% else
%     fprintf('Depth alignment OK: centres differ by %.2f mm.\n', depth_diff);
% end
% 
% nC = numel(roiC_pixels);
% nS = numel(roiS_pixels);
% fprintf('ROI C size: %d pixels | ROI S size: %d pixels\n', nC, nS);
% if min(nC,nS) < 20
%     warning('ROI too small — increase for reliable statistics.');
% end
% 
% % ── Compute CNR ──────────────────────────────────────────────────────────
% mu_c     = mean(roiC_pixels);
% mu_s_cnr = mean(roiS_pixels);
% sigma2_c     = var(roiC_pixels);
% sigma2_s_cnr = var(roiS_pixels);
% 
% CNR_linear = abs(mu_s_cnr - mu_c) / sqrt((sigma2_s_cnr + sigma2_c) / 2);
% CNR_dB     = 20 * log10(CNR_linear);
% 
% % Contrast (useful companion metric)
% contrast_lin = mu_s_cnr / max(mu_c, eps);
% contrast_dB  = 20 * log10(contrast_lin);
% 
% fprintf('\n--- Results (L7-4 CNR) ---\n');
% fprintf('  mu_c     (hypoechoic mean):    %.4e\n', mu_c);
% fprintf('  mu_s     (background mean):    %.4e\n', mu_s_cnr);
% fprintf('  sigma_c  (hypoechoic std):     %.4e\n', sqrt(sigma2_c));
% fprintf('  sigma_s  (background std):     %.4e\n', sqrt(sigma2_s_cnr));
% fprintf('  sigma/mu (hypoechoic):         %.2f   (Rayleigh expects 0.52)\n', sqrt(sigma2_c)/mu_c);
% fprintf('  sigma/mu (background):         %.2f   (Rayleigh expects 0.52)\n', sqrt(sigma2_s_cnr)/mu_s_cnr);
% fprintf('  Contrast (mu_s / mu_c):        %.2f   (%.2f dB)\n', contrast_lin, contrast_dB);
% fprintf('  CNR (linear):                  %.4f\n', CNR_linear);
% fprintf('  CNR (dB):                      %.2f dB\n', CNR_dB);
% 
% title(sprintf('L7-4  |  CNR = %.2f dB  |  Green = Hypo, Blue = Background', CNR_dB));
% hold off;
% 
% % ── Histogram sanity check (same as section 7) ───────────────────────────
% figure(6); clf;
% all_pixels = [roiC_pixels; roiS_pixels];
% edges      = linspace(min(all_pixels), max(all_pixels), 41);
% 
% histogram(roiS_pixels, edges, ...
%     'Normalization','probability', ...
%     'FaceColor',[0.2 0.4 0.8], 'FaceAlpha',0.55, ...
%     'DisplayName','Background (ROI S)');
% hold on;
% histogram(roiC_pixels, edges, ...
%     'Normalization','probability', ...
%     'FaceColor',[0.2 0.7 0.2], 'FaceAlpha',0.55, ...
%     'DisplayName','Hypoechoic (ROI C)');
% 
% xline(mu_c,     '--', 'Color',[0.1 0.5 0.1], 'LineWidth',1.8, ...
%     'Label',sprintf('\\mu_c = %.2e', mu_c));
% xline(mu_s_cnr, '--', 'Color',[0.1 0.3 0.7], 'LineWidth',1.8, ...
%     'Label',sprintf('\\mu_s = %.2e', mu_s_cnr));
% 
% xlabel('Pixel intensity (linear, |IQ|)');
% ylabel('Probability');
% title(sprintf('CNR ROI histograms  |  CNR = %.2f dB', CNR_dB));
% legend('Location','northeast');
% grid on;
% hold off;
% 
% % ── Store result ─────────────────────────────────────────────────────────
% CNR_result_L74.mu_c         = mu_c;
% CNR_result_L74.mu_s         = mu_s_cnr;
% CNR_result_L74.sigma2_c     = sigma2_c;
% CNR_result_L74.sigma2_s     = sigma2_s_cnr;
% CNR_result_L74.CNR_linear   = CNR_linear;
% CNR_result_L74.CNR_dB       = CNR_dB;
% CNR_result_L74.contrast_dB  = contrast_dB;
% CNR_result_L74.roiC_pos     = posC;
% CNR_result_L74.roiS_pos     = posS;
% % Save with: save('CNR_L74_30mm.mat', 'CNR_result_L74');