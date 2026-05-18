%% 1. BUFFER SETUP
clear all; close all; clc;

% ---------------------------- Singule targets ---------------------------
% filename = "C:\Users\cioca100\Desktop\TX water tank\PulseEcho_Nylon\Focus 35mm (actual depth ~33p8)\NylonString_Rcv_35.mat";
% filename = "C:\Users\cioca100\Desktop\TX water tank\PulseEcho_SmallRod\SmallRod_Rcv.mat";
% filename = "C:\Users\cioca100\Desktop\TX water tank\PulseEcho_BigRod\BigRod_Rcv.mat";

% --------------------- Custom targets - Horizontal ----------------------
% filename = "C:\Users\cioca100\Desktop\WaterTank_Imaging\flat_lateral_res_test\1\T1_Rcv.mat";
% filename = "C:\Users\cioca100\Desktop\WaterTank_Imaging\flat_lateral_res_test\2\T2_Rcv.mat";
% filename = "C:\Users\cioca100\Desktop\WaterTank_Imaging\flat_lateral_res_test\3\T3_Rcv.mat";
% filename = "C:\Users\cioca100\Desktop\WaterTank_Imaging\flat_lateral_res_test\4\T4_Rcv.mat";  
% filename = "C:\Users\cioca100\Desktop\WaterTank_Imaging\USB stuff\WaterTank_Imaging\flat_lateral_res_test\5\T5_Rcv.mat";
% filename = "C:\Users\cioca100\Desktop\WaterTank_Imaging\USB stuff\WaterTank_Imaging\flat_lateral_res_test\6\T6_Rcv.mat";
% filename = "C:\Users\cioca100\Desktop\WaterTank_Imaging\USB stuff\WaterTank_Imaging\flat_lateral_res_test\7\T7_Rcv.mat";

% ---------------------- Custom Targets - Axial --------------------------
% filename = "C:\Users\cioca100\Desktop\WaterTank_Imaging\flat_axial_res_test\Axial_PlaneCustom_Rcv.mat";
% filename = "C:\Users\cioca100\Desktop\WaterTank_Imaging\flat_axial_res_test\Lateral_PlaneCustom_d1_Rcv.mat";
% filename = "C:\Users\cioca100\Desktop\WaterTank_Imaging\flat_axial_res_test\Lateral_PlaneCustom_d5_Rcv.mat";
% filename = "C:\Users\cioca100\Desktop\WaterTank_Imaging\flat_axial_res_test\Lateral_PlaneCustom_alt_Rcv.mat";
filename = "C:\Users\cioca100\Desktop\WaterTank_Imaging\flat_axial_res_test\Axial_PlaneCustom_2Targets_Rcv.mat";
% filename = "C:\Users\cioca100\Desktop\WaterTank_Imaging\flat_axial_res_test\Lateral_PlaneCustom_2Targets_Rcv.mat";

% ------------ Custom targets - Axial/Lateral (offset/unshadowed) --------
% filename = "C:\Users\cioca100\Desktop\WaterTank_Imaging\USB stuff\WaterTank_Imaging\offset_axiallat_res_test\1\S1_Rcv.mat";
% filename = "C:\Users\cioca100\Desktop\WaterTank_Imaging\USB stuff\WaterTank_Imaging\offset_axiallat_res_test\2(extra)\S2_Rcv.mat";

% ------------------------ Phantom Imaging -------------------------------
% filename = "C:\Users\cioca100\Desktop\TX water tank\Phantom_Imaging_ChipMUT\9p5Sc\Phantom_NF_2_Rcv.mat";

% filename = "C:\Users\cioca100\Desktop\WaterTank_Imaging\flat_lateral_res_test\4\T4_Rcv.mat";

fprintf('Loading RF data...\n');
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
CompleteFrame_unavg  = RF_Matrix(:, 65:128, :);
num_frames_to_avg    = min(30, frames);
CompleteFrame        = mean(CompleteFrame_unavg(:, :, 1:num_frames_to_avg), 3);

% ── CHANNEL REORDERING ───────────────────────────────────────────────────
% Your PCB routes channels in an interleaved pattern defined by ConnectorES:
%   ch_a = 33:64;  ch_b = 1:32;
%   ConnectorES = 64 + reshape([ch_a;ch_b], 1,[])
%               = [97,65, 98,66, 99,67, ..., 128,96]
%
% This means in our extracted 64 columns (pins 65–128):
%   Our column k = pin (64+k)
%   Physical element 2j-1 (odd)  is on pin 96+j  → our column 32+j
%   Physical element 2j   (even) is on pin 64+j  → our column j
%
% reorder_idx(i) = which extracted column holds physical element i
reorder_idx = zeros(1, 64);
for k = 1:32
    reorder_idx(2*k - 1) = 32 + k;   % odd physical elements: col 33,34,...,64
    reorder_idx(2*k)     = k;         % even physical elements: col 1,2,...,32
end

% Reorder columns so that CompleteFrame(:,i) = physical element i
CompleteFrame = CompleteFrame(:, reorder_idx);
fprintf('Channel reordering applied: interleaved → sequential physical order.\n');
% ─────────────────────────────────────────────────────────────────────────
%% 2. MATCHED FILTER
% Your transmit waveform: 8 half-cycles, 67% duty cycle, fc = 10.5 MHz.
% The matched filter is the time-reversed conjugate of the transmit pulse.
% This compresses the pulse axially and improves SNR vs. Y-diffraction clutter.
fs = 42e6;
fc = 10.5e6;
fprintf('Applying matched filter...\n');

% Construct the reference transmit pulse (8 half-cycles, 67% duty cycle)
T_pulse      = 8 / (2 * fc);               % total pulse duration (s)
t_pulse      = 0 : 1/fs : T_pulse;
duty         = 0.67;
ref_pulse    = double(mod(t_pulse * fc, 1) < duty) .* sin(2*pi*fc*t_pulse);
ref_pulse    = ref_pulse / norm(ref_pulse); % normalise

% Apply matched filter to each channel (correlation = convolution with time-reversed pulse)
mf_kernel    = fliplr(conj(ref_pulse));     % time-reversed conjugate
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

[b_lp, a_lp] = butter(3, (fc * 0.8) / (fs/2), 'low');
IQ_Matrix     = filtfilt(b_lp, a_lp, double(IQ_Matrix));

%% 4. IQ BASEBAND BEAMFORMING (PLANE WAVE)
fprintf('Beamforming Setup...\n');
c      = 1480;
lambda = c / fc;

% Set to Receive(1).startDepth from your Verasonics script (in wavelengths)
start_depth_wl = 18;           % <── TUNE THIS
t_start        = start_depth_wl / fc;   % correct formula (no factor of 2)

pitch         = 0.075e-3;
element_width = 0.051e-3;
f_number      = 4.0; % standard is 4.0

probe_x       = ((0:nchannels-1) - (nchannels-1)/2) * pitch;
half_aperture = (nchannels * pitch) / 2;

z_axis = linspace(10e-3, 80e-3, 1024);
x_axis = linspace(-10e-3, 10e-3, 256);
[X_pix, Z_pix] = meshgrid(x_axis, z_axis);

fprintf('Beamforming (%d channels)...\n', nchannels);

dist_tx = Z_pix;   % plane wave — correct for TX F# = 8.3

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

    % Phase rotation relative to t_start (not absolute tau)
    phase_rotation = exp(1i * 2 * pi * fc * (tau(valid_mask) - t_start));

    iq_aligned = zeros(size(X_pix));
    iq_aligned(valid_mask) = iq_delayed_smooth(valid_mask) .* phase_rotation;

    weights        = zeros(size(X_pix));
    denom_radius   = max(max_radius, half_aperture);
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

array_edge_mm  = 2.4;
line_thickness = 2.0;
line_color     = 'w';
line_style     = ':';

figure(1); clf;
imagesc(x_axis*1000, z_axis*1000, B_Mode_dB);
clim([-30 0]);
colorbar;
title('PMUT: MF + IQ Demod + Directivity + CF');
xlabel('Lateral (mm)'); ylabel('Axial (mm)');
axis image;
hold on;
xline(-array_edge_mm, 'Color', line_color, 'LineStyle', line_style, 'LineWidth', line_thickness);
xline( array_edge_mm, 'Color', line_color, 'LineStyle', line_style, 'LineWidth', line_thickness);
hold off;

figure(2); clf;
imagesc(x_axis*1000, z_axis*1000, B_Mode_Lin);
colorbar;
title('PMUT: MF + IQ Demod + Directivity + CF (Linear)');
xlabel('Lateral (mm)'); ylabel('Axial (mm)');
axis image;