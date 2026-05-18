function compute_pressure(txtFile, calCsv)

% ============================ Settings ==================================
V_to_mV = 1000;  % mV = 1000 | V = 1

useWindow  = true;
removeMean = true;
nfftMode   = "custom"; % "nextpow2" | "signal_length" | "custom"
nfftCustom = 65536;      % 8192 | 16384 | 32768 | 65536

plotMaxMHz = 30;

% ---- Toneburst gating options ----
useGate      = true;
gateMode     = "auto";    % "auto" (envelope threshold) or "manual"
gate_us      = 0.6;       % standardized gate length (us)
gatePad_us   = 0.2;       % padding (us) around detected burst
envThresh    = 0.4;       % threshold as fraction of max envelope (auto)
manual_t0_us = [];        % gate start time (us) if gateMode="manual"

f0 = 10.5e6;
useBandpassForGate = true;
bpFracBW           = 0.35;
minIsland_us       = 0.05;
searchWindow_us    = [];

% NEW: enable/disable time-domain reconstruction
doTimeRecon = true;

% ========================================================================

% ====================== Load raw voltage waveform =======================
[t_s, v_raw] = read_time_series_with_header(txtFile);

v_mV = v_raw * V_to_mV;

dt = median(diff(t_s));
fs = 1/dt;
fprintf("Estimated sampling rate: %.3f GSa/s\n", fs/1e9)

if removeMean
    v_mV = v_mV - mean(v_mV);
end

% ====================== Toneburst gating (IMPROVED) ======================
% CHANGED: initialize for plotting robustness
env = [];
thr = NaN;

if useGate
    Ng_target = round(gate_us*1e-6 * fs);
    pad       = round(gatePad_us*1e-6 * fs);
    minIsland = round(minIsland_us*1e-6 * fs);

    switch gateMode
        case "auto"
            x = v_mV;

            % Optional bandpass before envelope detection
            if useBandpassForGate
                fL = max(0.1e6, (1 - bpFracBW)*f0); % keep >0
                fH = (1 + bpFracBW)*f0;
                x = bandpass(x, [fL fH], fs);
            end

            % Optional search window restriction
            if ~isempty(searchWindow_us)
                t1 = searchWindow_us(1)*1e-6;
                t2 = searchWindow_us(2)*1e-6;
                mask = (t_s >= t1) & (t_s <= t2);
            else
                mask = true(size(t_s));
            end

            env = abs(hilbert(x));             % envelope
            env(~mask) = 0;                    % ignore outside search window

            % Find strongest event (global max)
            [~, imax] = max(env);

            % Local threshold based on peak
            thr = envThresh * env(imax);

            % Grow left/right from imax to find contiguous "island"
            iL = imax;
            while iL > 1 && env(iL) > thr
                iL = iL - 1;
            end
            iR = imax;
            while iR < numel(env) && env(iR) > thr
                iR = iR + 1;
            end

            % Add padding
            iL = max(1, iL - pad);
            iR = min(numel(env), iR + pad);

            % If island too small, fall back to fixed window around imax
            if (iR - iL + 1) < minIsland
                warning("Detected burst island too small; using fixed gate around max.");
                ic = imax;
            else
                ic = round((iL + iR)/2);
            end

            % Final fixed-length gate centered at ic
            i1 = max(1, ic - floor(Ng_target/2));
            i2 = min(numel(v_mV), i1 + Ng_target - 1);

        case "manual"
            if isempty(manual_t0_us)
                error("manual_t0_us must be set when gateMode='manual'.");
            end
            i1 = round((manual_t0_us*1e-6 - t_s(1)) * fs) + 1;
            i1 = max(1, i1);
            i2 = min(numel(v_mV), i1 + Ng_target - 1);

        otherwise
            error("Invalid gateMode. Use 'auto' or 'manual'.");
    end

    v_gate = v_mV(i1:i2);
    t_gate = t_s(i1:i2);

else
    v_gate = v_mV;
    t_gate = t_s;
    i1 = 1; i2 = numel(v_mV);  % NEW: for plotting gate lines consistently
end

Ng = numel(v_gate);
fprintf("FFT segment samples (Ng) = %d (%.3f us)\n", Ng, Ng/fs*1e6);

% Plot raw + envelope + threshold + gate (robust if env/thr empty)
figure('Color','w');
plot(t_s*1e6, v_mV, 'b'); hold on;

if ~isempty(env)
    envScaled = env / max(env) * max(abs(v_mV));
    plot(t_s*1e6, envScaled, 'k'); % scaled env overlay

    if ~isnan(thr) && max(env) > 0
        thrScaled = thr / max(env) * max(abs(v_mV));
        yline(thrScaled, '--k', 'thr'); % scaled threshold
    end
end

xline(t_s(i1)*1e6, '--r', 'Gate start');
xline(t_s(i2)*1e6, '--r', 'Gate end');
grid on;
xlabel('Time (us)'); ylabel('mV');
title('Raw + envelope (scaled) + threshold + gate');
legend('raw','env (scaled)','thr (scaled)','Gate start','Gate end','Location','best');

% ========================================================================

% ===================== Windowing (apply to gated) ========================
if useWindow
    w  = hann(Ng);
    cg = mean(w);
    v_win = v_gate .* w;
else
    cg = 1.0;
    v_win = v_gate;
end

% ======================= FFT of voltage spectrum =========================
switch nfftMode
    case "nextpow2"
        nfft = 2^nextpow2(Ng);      % use Ng not full record length
    case "signal_length"
        nfft = Ng;
    case "custom"
            nfft = nfftCustom;
    otherwise
        if isnumeric(nfftMode)
            nfft = nfftMode;
        else
            error("Invalid nfftMode. Use 'nextpow2', 'signal_length', or numeric.");
        end
end

V = fft(v_win, nfft);
nUnique = floor(nfft/2) + 1;

V1 = V(1:nUnique);
f_Hz = (0:nUnique-1) * (fs/nfft);

% scale by Ng (gated length)
Vmag_mV = abs(V1) * (2/(Ng*cg));
Vmag_mV(1) = Vmag_mV(1)/2;
if mod(nfft,2)==0
    Vmag_mV(end) = Vmag_mV(end)/2;
end

% ===================== Hydrophone Sensitivity ============================
Tcal = readtable(calCsv);
fcal_MHz = Tcal{:,1};
Scal_mV_per_MPa = Tcal{:,2};

if any(Scal_mV_per_MPa <= 0)
    warning("Calibration sensitivity contains non-positive values. Check CSV.");
end

% Interpolate sensitivity onto FFT frequency axis (single-sided, for spectra)
f_MHz = f_Hz / 1e6;
Sinterp_mV_per_MPa = interp1(fcal_MHz, Scal_mV_per_MPa, f_MHz, "linear", "extrap");

% ========================================================================
% NEW: Time-domain reconstruction p(t) from complex spectrum (IFFT)
% ========================================================================
if doTimeRecon
    % Certificate provides sensitivity over 1–30 MHz (typical) [1](https://imecinternational-my.sharepoint.com/personal/cioca100_imec_be/Documents/Microsoft%20Copilot%20Chat%20Files/PA18082%20+%202864%20Sensitivity%2004April24.pdf)
    fcalMin_MHz = min(fcal_MHz);
    fcalMax_MHz = max(fcal_MHz);

    % Use gated but UNWINDOWED voltage segment to avoid tapering peaks
    Vrec = fft(v_gate, nfft);  % complex FFT, length nfft

    % Build full-spectrum frequency axis and map to [0..fs/2]
    k = 0:(nfft-1);
    f_full_Hz = k * (fs/nfft);
    f_pos_Hz  = min(f_full_Hz, fs - f_full_Hz);  % mirror
    f_pos_MHz = f_pos_Hz / 1e6;

    % Interpolate sensitivity magnitude; no extrapolation
    Sfull = interp1(fcal_MHz, Scal_mV_per_MPa, f_pos_MHz, "linear", NaN);

    % Valid only inside calibrated band and positive sensitivity
    valid = (f_pos_MHz >= fcalMin_MHz) & (f_pos_MHz <= fcalMax_MHz) & (Sfull > 0);

    % Outside calibration band: set Inf so division yields 0 (band-limiting)
    Sfull(~valid) = Inf;
    Sfull(1) = Inf;  % DC bin

    % Convert to pressure spectrum (MPa), then IFFT back to p(t)
    Prec_MPa_spec = Vrec ./ Sfull(:);

    p_rec_MPa_full = real(ifft(Prec_MPa_spec, nfft));
    p_rec_MPa = p_rec_MPa_full(1:Ng);  % match gated segment length

    p_rec_kPa = p_rec_MPa * 1000;
    t_rec = t_gate;

    % Time-domain metrics
    p_plus_kPa  = max(p_rec_kPa);
    p_minus_kPa = min(p_rec_kPa);
    p_pp_kPa    = p_plus_kPa - p_minus_kPa;
    
    % ----------------------------- optional bp filtering for envelope
    f0 = 10.5e6;
    bpFracBW = 0.35;  % like your gate option
    p_bp = bandpass(p_rec_kPa, [(1-bpFracBW)*f0 (1+bpFracBW)*f0], fs);
    env_p = abs(hilbert(p_bp));
    % -----------------------------------------------------------------
    
    % Uncomment if you're not using the filter
    % env_p = abs(hilbert(p_rec_kPa)); 
    env_peak_kPa = max(env_p);

    fprintf("\n--- Time-domain pressure metrics (reconstructed) ---\n");
    fprintf("p+   = %.4g kPa\n", p_plus_kPa);
    fprintf("p-   = %.4g kPa\n", p_minus_kPa);
    fprintf("p_pp = %.4g kPa\n", p_pp_kPa);
    fprintf("Envelope peak = %.4g kPa\n", env_peak_kPa);

    % Plot reconstructed pressure waveform + envelope
    figure('Color','w');
    plot(t_rec*1e6, p_rec_kPa, 'b', 'LineWidth', 1.2); hold on;
    plot(t_rec*1e6, env_p, 'k--', 'LineWidth', 1.0);
    yline(p_plus_kPa,  ':r', 'p+');
    yline(p_minus_kPa, ':r', 'p-');
    grid on;
    xlabel('Time (\mus)');
    ylabel('Pressure (kPa)');
    title('Reconstructed time-domain pressure waveform p(t) + envelope');
    legend('p(t) reconstructed','Envelope |hilbert(p(t))|','Location','best');
end
% ========================================================================

% ========================== Pressure Spectrum ============================
Vmag_mV = Vmag_mV(:);
Sinterp_mV_per_MPa = Sinterp_mV_per_MPa(:);

Pmag_MPa = Vmag_mV ./ Sinterp_mV_per_MPa;
Pmag_kPa = Pmag_MPa * 1000;

% ======================== Pressure per Volt ==============================
Vdrive_Vpk = 10;  % set/measure this at the transducer terminals if possible
Kspec_kPa_per_V = Pmag_kPa / Vdrive_Vpk;

% Extract value at f0 (optional scalar report)
f0_MHz = 10.50;
f0 = f0_MHz * 1e6;
[~, idx0] = min(abs(f_Hz - f0));
p0_kPa = Pmag_kPa(idx0);
K0_kPa_per_V = Kspec_kPa_per_V(idx0);

fprintf("At %.3f MHz: p = %.4g kPa, Vdrive = %.4g Vpk -> K = %.4g kPa/V\n", ...
    f0_MHz, p0_kPa, Vdrive_Vpk, K0_kPa_per_V);

% ============================ Figures/Plots ==============================
% Plot voltage magnitude spectrum
figure('Color','w');
plot(f_MHz, Vmag_mV, 'LineWidth', 1.2);
grid on;
xlabel("Frequency (MHz)");
ylabel("|V(f)| (mV peak)");
title("Voltage magnitude spectrum (gated toneburst)", 'Interpreter','none');
if ~isempty(plotMaxMHz), xlim([0 plotMaxMHz]); end

% Plot pressure spectrum
figure('Color','w');
plot(f_MHz, Pmag_kPa, 'LineWidth', 1.2);
grid on;
xlabel("Frequency (MHz)");
ylabel("|P(f)| (kPa peak)");
title("Pressure magnitude spectrum (gated toneburst)", 'Interpreter','none');
if ~isempty(plotMaxMHz), xlim([0 plotMaxMHz]); end

% Plot kPa/V spectrum
figure('Color','w');
plot(f_MHz, Kspec_kPa_per_V, 'LineWidth', 1.2);
grid on;
xlabel("Frequency (MHz)");
ylabel("kPa/V (peak)");
title("Pressure per Volt spectrum (gated toneburst)", 'Interpreter','none');
if ~isempty(plotMaxMHz), xlim([0 plotMaxMHz]); end

% Plot interpolated sensitivity curve
figure('Color','w');
plot(f_MHz, Sinterp_mV_per_MPa, 'LineWidth', 1.2);
grid on;
xlabel("Frequency (MHz)");
ylabel("mV/Mpa (peak)");
title("Mario hydrophone (intrepolated) Sensitivity Curve)", 'Interpreter','none');
if ~isempty(plotMaxMHz), xlim([0 plotMaxMHz]); end

% ===================== Output Summary ===================================
fprintf("\nDone.\n");
fprintf("Gate samples: Ng = %d\n", Ng);
fprintf("FFT length  : nfft = %d\n", nfft);
fprintf("Freq step   : df = %.3f kHz\n", (fs/nfft)/1e3);
fprintf("Nyquist     : %.2f MHz\n", (fs/2)/1e6);

end


%% =======================================================================
% ************************* Parse .txt helper ****************************
% ========================================================================
% Reads a file that contains a header/metadata block followed by numeric rows.
% Expected numeric columns: time(s), voltage, (optional third column) 
function [t, y] = read_time_series_with_header(filePath)

txt = fileread(filePath);
lines = splitlines(string(txt));

% Find the first numeric line
startIdx = 0;
for i = 1:numel(lines)
    s = strtrim(lines(i));
    if s == "" ; continue; end
    if ~isempty(regexp(s, '^[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?\s+', 'once'))
        startIdx = i; break;
    end
end
if startIdx == 0
    error("No numeric data found in %s", filePath);
end

dataBlock = strjoin(lines(startIdx:end), newline);

% Read up to 3 columns. If 3rd exists, it is ignored.
C = textscan(dataBlock, '%f %f %f', ...
    'Delimiter', {' ', '\t'}, ...
    'MultipleDelimsAsOne', true, ...
    'CollectOutput', true);

A = C{1};
if isempty(A) || size(A,2) < 2
    error("Expected at least 2 columns (time, voltage) in %s", filePath);
end

t = A(:,1);
y = A(:,2);
end



%% =======================================================================
% ************************* main  ****************************
% ========================================================================

% Array 20 - @ freq = 10MHz; Beam Max; depth 55mm, 37us
% txtFile = "C:\Users\cioca100\Desktop\TX water tank\Alex and Maarten\TX measurements water-tank\array20\Waveform_BeamMax_depth_55mm_37us\waveform_depth_55mm_37us.txt";
% Array 20 - @ freq = 10MHz; depth 55mm, 37us
% txtFile = "C:\Users\cioca100\Desktop\TX water tank\Alex and Maarten\TX measurements water-tank\array20\Waveform_depth_55mm_37us\waveform_depth_55mm_37us.txt";
% Array 20 - @ freq = 10.5MHz; Beam Max; depth 55mm, 37us
% txtFile = "C:\Users\cioca100\Desktop\TX water tank\Alex and Maarten\TX measurements water-tank\array20\Waveform_BeamMax_depth_55mm_37us_10_5MHz\waveform_depth_55mm_37us_10_5MHz.txt";

% Array 27 - @ freq = 10.5MHz; Beam Max; depth 55mm, 37us (300kb)
% txtFile = "C:\Users\cioca100\Desktop\TX water tank\Alex and Maarten\TX measurements water-tank\array27_connected_ground\freqsweep_10p4MHz\10p50MHz.txt"; 
% Array 27 - @ refreq = 10.5MHz; Beam Max; depth 55mm, 37us (149kb - redone)
txtFile = "C:\Users\cioca100\Desktop\TX water tank\Alex and Maarten\TX measurements water-tank\array27_connected_ground\Waveform_depth 55mm_37us_10_5MHz\waveform_depth_55mm_37us_10_5MHz.txt";

% Array 21 - @ freq = 10MHz; Beam Max; depth 55mm, 37us
% txtFile = "C:\Users\cioca100\Desktop\TX water tank\Alex and Maarten\TX measurements water-tank\array21\Waveform\wavefrom_depth_55mm_37us.txt";
% Array 21 - @ freq = 10.5MHz; Beam Max; depth 55mm, 37us
% txtFile = "C:\Users\cioca100\Desktop\TX water tank\freqsweep_manual\10p5MHz.txt";

% Signal sent to drive the PMUT array
% txtFile = "C:\Users\cioca100\Desktop\TX water tank\Signal_Sent_to_Drive_PMUT\5Vamplitude_10MHz_verasonics_signal.txt";

% Sanity Check Z-axis Linear Scan (5th point) - theory: 2.59kPa/V
% txtFile = "C:\Users\cioca100\Desktop\TX water tank\Alex and Maarten\TX measurements water-tank\array27_connected_ground\Z_scan (10_5MHz - Linear) - axial profile\Scan_ZZ004.txt"

% Mario (Calibration) Sensitivity Curve
calCsv  = "C:\Users\cioca100\Desktop\TX water tank\Alex and Maarten\TX measurements water-tank\array27_connected_ground\Hydrophone_Mario_sensitivity_curve.csv";   % 30-point calibration curve

compute_pressure(txtFile, calCsv)