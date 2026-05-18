"""
Comprehensive matched-filter figure for thesis appendix, using the real
RF trace from rf_trace.mat.

Layout (4 rows x 2 cols):
  Row 1: (a) reference square-wave toneburst, time domain
         (b) reference toneburst, frequency domain
  Row 2: (c) raw RF (real echo + noise), time domain
         (d) raw RF, frequency domain
  Row 3: (e) explanation: cross-correlation as a sliding inner product
         (f) operation in equations
  Row 4: (g) matched-filter output, time domain
         (h) matched-filter output, frequency domain

Same robustness logic as the rest of the appendix scripts:
  - blanks the saturating transmit transient at the start of the trace
  - searches for the actual target echo inside a configurable time
    window so the matched-filter peak detection doesn't lock onto the
    transmit-coupling spike
"""
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.signal import hilbert
from scipy.io import loadmat
import textwrap

# ============ USER CONFIG =====================================================
RF_FILE = 'rf_trace.mat'
OUT_DIR = 'out'

BLANK_END_US      = 5.0       # zero out everything before this time
ECHO_SEARCH_LO_US = 40.0      # echo must be inside this window
ECHO_SEARCH_HI_US = 50.0
# ==============================================================================

os.makedirs(OUT_DIR, exist_ok=True)

mpl.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['DejaVu Serif', 'Times New Roman', 'Times'],
    'font.size': 10, 'axes.labelsize': 10, 'axes.titlesize': 10,
    'xtick.labelsize': 8.5, 'ytick.labelsize': 8.5,
    'axes.linewidth': 0.8,
    'xtick.major.width': 0.8, 'ytick.major.width': 0.8,
    'xtick.major.size': 3,    'ytick.major.size': 3,
    'xtick.direction': 'in',  'ytick.direction': 'in',
    'mathtext.fontset': 'cm',
    'axes.spines.top': False, 'axes.spines.right': False,
    'pdf.fonttype': 42, 'ps.fonttype': 42,
})

C_RAW   = '#1f4e79'
C_PULSE = '#2f6f4f'
C_MF    = '#a8323b'
C_FILL_RAW = '#cdd5dc'
C_FILL_MF  = '#e6cdd1'
C_FILL_PUL = '#cfe0d4'
C_GREY  = '#5F5E5A'
C_ACCENT = '#0e5a73'

# ------------------------------------------------------------------
# Hardware parameters
# ------------------------------------------------------------------
fc = 10.5e6
fs = 42e6

# ------------------------------------------------------------------
# Reference square-wave toneburst (8 half-cycles = 4 full periods at fc)
# This is the actual TX drive waveform; the matched filter correlates
# against this.
# ------------------------------------------------------------------
T_pulse = 4 / fc
t_pulse = np.arange(0, T_pulse + 1/fs, 1/fs)
duty = 0.67
ref_pulse_sq = np.where(np.mod(t_pulse * fc, 1) < duty, 1.0, -1.0)
ref_pulse_sq_n = ref_pulse_sq / np.linalg.norm(ref_pulse_sq)

# ------------------------------------------------------------------
# Load real RF trace + blank the transmit transient
# ------------------------------------------------------------------
rf_full = loadmat(RF_FILE)['trace'].squeeze().astype(float)
N = len(rf_full)
t = np.arange(N) / fs
t_us = t * 1e6
print(f"Loaded RF trace: {N} samples = {N/fs*1e6:.1f} us")

i_blank_end = int(BLANK_END_US * 1e-6 * fs)
rf = rf_full.copy()
rf[:i_blank_end] = 0.0
print(f"Blanked samples 0..{i_blank_end} (t < {BLANK_END_US} us)")

# Matched filter
mf_kernel = ref_pulse_sq_n[::-1]
rf_mf = np.convolve(rf, mf_kernel, mode='same')

# Per-trace normalisation (works because the transmit transient is gone)
rf_disp    = rf    / np.max(np.abs(rf))
rf_mf_disp = rf_mf / np.max(np.abs(rf_mf))
env_mf     = np.abs(hilbert(rf_mf_disp))

# Echo location: argmax of MF envelope inside the search window
i_search_lo = int(ECHO_SEARCH_LO_US * 1e-6 * fs)
i_search_hi = int(ECHO_SEARCH_HI_US * 1e-6 * fs)
i_echo = i_search_lo + int(np.argmax(env_mf[i_search_lo:i_search_hi]))
t_echo_us = i_echo / fs * 1e6
print(f"Echo search window: {ECHO_SEARCH_LO_US}-{ECHO_SEARCH_HI_US} us")
print(f"Detected echo: i = {i_echo}, t = {t_echo_us:.2f} us, "
      f"env = {env_mf[i_echo]:.3f}")

# Display window: +/- 9 us around the echo
window_us = 9.0
t_lo = max(0.0, t_echo_us - window_us)
t_hi = min(t_us[-1], t_echo_us + window_us)
window = (t_us >= t_lo) & (t_us <= t_hi)

# ------------------------------------------------------------------
# Spectra
# ------------------------------------------------------------------
def single_sided(x, fs):
    X = np.fft.fft(x)
    f = np.fft.fftfreq(len(x), 1/fs)
    pos = f >= 0
    Xm = np.abs(X[pos])
    return f[pos] / 1e6, Xm / np.max(Xm)

# Reference pulse spectrum (zero-padded for smoother frequency resolution)
ref_padded = np.concatenate([ref_pulse_sq_n, np.zeros(2048 - len(ref_pulse_sq_n))])
f_pulse, X_pulse = single_sided(ref_padded, fs)

# RF and MF spectra: take a segment around the detected echo
n_pre, n_post = 200, 1500
i_seg_lo = max(0, i_echo - n_pre)
i_seg_hi = min(len(rf_mf), i_echo + n_post)
rf_seg    = rf   [i_seg_lo:i_seg_hi]
rf_mf_seg = rf_mf[i_seg_lo:i_seg_hi]
f_rf, X_rf      = single_sided(rf_seg,    fs)
f_mf, X_mf_spec = single_sided(rf_mf_seg, fs)

# ==================================================================
# FIGURE
# ==================================================================
fig = plt.figure(figsize=(8.4, 11.4))
gs = fig.add_gridspec(4, 2, hspace=0.85, wspace=0.30,
                      left=0.08, right=0.97, top=0.965, bottom=0.045,
                      width_ratios=[1.6, 1])

# def panel_label(ax, tag, title, subtitle=None, wrap=None):
#     ax.text(0.0, 1.20, f'({tag})  {title}',
#             transform=ax.transAxes,
#             fontsize=10.5, weight='bold', color='#2C2C2A')
#     if subtitle:
#         if wrap:
#             subtitle = textwrap.fill(subtitle, width=wrap)
#         ax.text(0.0, 1.05, subtitle,
#                 transform=ax.transAxes,
#                 fontsize=8.5, color=C_GREY, style='italic',
#                 verticalalignment='bottom')   # <-- key bit

def panel_label(ax, tag, title, subtitle=None, wrap=None):
    if subtitle and wrap:
        subtitle = textwrap.fill(subtitle, width=wrap)

    n_extra = subtitle.count('\n') if subtitle else 0
    title_y = 1.20 + 0.13 * n_extra   # bump per wrapped line

    ax.text(0.0, title_y, f'({tag})  {title}',
            transform=ax.transAxes,
            fontsize=10.5, weight='bold', color='#2C2C2A')
    if subtitle:
        ax.text(0.0, 1.05, subtitle,
                transform=ax.transAxes,
                fontsize=8.5, color=C_GREY, style='italic',
                verticalalignment='bottom')
        
# ------------------------------------------------------------------
# Row 1: reference pulse (a) time, (b) frequency
# ------------------------------------------------------------------
ax_a = fig.add_subplot(gs[0, 0])
fs_fine = 50 * fc
t_fine = np.arange(0, T_pulse + 1/fs_fine, 1/fs_fine)
ref_fine = np.where(np.mod(t_fine * fc, 1) < duty, 1.0, -1.0)
ax_a.fill_between(t_fine * 1e6, 0, ref_fine, color=C_FILL_PUL, lw=0)
ax_a.plot(t_fine * 1e6, ref_fine, color=C_PULSE, lw=1.2,
          drawstyle='steps-post')
ax_a.plot(t_pulse * 1e6, ref_pulse_sq, 'o', color=C_PULSE, ms=3.8,
          mfc='white', mew=1.0)
ax_a.axhline(0, color=C_GREY, lw=0.4)
ax_a.set_xlim(0, T_pulse * 1e6)
ax_a.set_ylim(-1.4, 1.4)
ax_a.set_xlabel(r'Time ($\mu$s)', labelpad=2)
ax_a.set_ylabel('Drive level', labelpad=4)
panel_label(ax_a, 'a', 'Reference pulse, time domain',
            r'Square wave at $f_c=10.5$ MHz, 67% duty, 8 half-cycles.')

ax_b = fig.add_subplot(gs[0, 1])
ax_b.fill_between(f_pulse, 0, X_pulse, color=C_PULSE, alpha=0.18, lw=0)
ax_b.plot(f_pulse, X_pulse, color=C_PULSE, lw=1.2)
ax_b.axvline(fc/1e6, color=C_ACCENT, lw=0.7, ls='--')
ax_b.text(fc/1e6 + 0.4, 0.92, r'$f_c$', fontsize=9, color=C_ACCENT)
for k, lbl in [(3, '3'), (5, '5')]:
    fk = k * fc / 1e6
    if fk < 21:
        ax_b.axvline(fk, color='#888', lw=0.5, ls=':')
        ax_b.text(fk + 0.3, 0.5, f'{lbl}$f_c$', fontsize=8, color='#666')
ax_b.set_xlim(0, 21)
ax_b.set_ylim(0, 1.08)
ax_b.set_xlabel('Frequency (MHz)', labelpad=2)
ax_b.set_ylabel('Norm. magnitude', labelpad=4)
# panel_label(ax_b, 'b', 'Reference pulse, frequency',
#             'Fundamental at $f_c$ plus odd harmonics (square-wave signature).')
panel_label(ax_b, 'b', 'Reference pulse, frequency',
            'Fundamental at $f_c$ plus odd harmonics (square-wave signature).',
            wrap=45)

# ------------------------------------------------------------------
# Row 2: raw RF (c) time, (d) frequency
# ------------------------------------------------------------------
ax_c = fig.add_subplot(gs[1, 0])
ax_c.plot(t_us[window], rf_disp[window], color=C_RAW, lw=0.6)
ax_c.axhline(0, color=C_GREY, lw=0.4)
# Annotate the echo position (single echo only -- real data, no synthetic
# late echo).
ax_c.annotate('echo', xy=(t_echo_us, 0.55),
              xytext=(t_echo_us, 1.05),
              ha='center', fontsize=8, color='#444',
              arrowprops=dict(arrowstyle='-', color='#888', lw=0.5))
ax_c.set_xlim(t_lo, t_hi)
ax_c.set_ylim(-1.2, 1.4)
ax_c.set_xlabel(r'Time ($\mu$s)', labelpad=2)
ax_c.set_ylabel('Norm. amp.', labelpad=4)
panel_label(ax_c, 'c', 'Raw RF, single channel',
            'Hydrophone voltage. Echo lasts ~3 $\\mu$s due to ringdown.')

ax_d = fig.add_subplot(gs[1, 1])
ax_d.fill_between(f_rf, 0, X_rf, color=C_RAW, alpha=0.18, lw=0)
ax_d.plot(f_rf, X_rf, color=C_RAW, lw=1.2)
ax_d.axvline(fc/1e6, color=C_ACCENT, lw=0.7, ls='--')
ax_d.text(fc/1e6 + 0.4, 0.92, r'$f_c$', fontsize=9, color=C_ACCENT)
ax_d.set_xlim(0, 21)
ax_d.set_ylim(0, 1.08)
ax_d.set_xlabel('Frequency (MHz)', labelpad=2)
ax_d.set_ylabel('Norm. magnitude', labelpad=4)
# panel_label(ax_d, 'd', 'Raw RF, frequency',
#             'Echo energy concentrated at $f_c$ (transducer is bandlimited). Noise spread across the band.')
panel_label(ax_d, 'd', 'Raw RF, frequency',
            'Echo energy concentrated at $f_c$ (transducer is bandlimited). '
            'Noise spread across the band.',
            wrap=45)

# ------------------------------------------------------------------
# Row 3: cross-correlation explanation
# ------------------------------------------------------------------
ax_e = fig.add_subplot(gs[2, 0])

def overlay(ax, idx_start, color, alpha=0.7, label=None):
    n = len(ref_pulse_sq_n)
    t_kern_us = (np.arange(n) + idx_start) / fs * 1e6
    k_disp = ref_pulse_sq_n / np.max(np.abs(ref_pulse_sq_n)) * 0.5
    ax.plot(t_kern_us, k_disp, color=color, lw=1.0, alpha=alpha,
            drawstyle='steps-post', label=label)

# Local zoom for the explanation panel: tighter than the main display
# so the kernel overlay is visible alongside the echo.
e_lo = max(0.0, t_echo_us - 4.0)
e_hi = min(t_us[-1], t_echo_us + 4.0)
e_mask = (t_us >= e_lo) & (t_us <= e_hi)
ax_e.plot(t_us[e_mask], rf_disp[e_mask], color=C_RAW, lw=0.5, alpha=0.7)

# Place the kernel at the detected echo (high correlation) and a few
# microseconds earlier in the noise (low correlation).
idx_aligned    = i_echo
idx_misaligned = max(0, i_echo - 130)   # ~3 us before
overlay(ax_e, idx_aligned,    C_PULSE, alpha=0.95,
        label='kernel aligned with echo')
overlay(ax_e, idx_misaligned, '#888',  alpha=0.85,
        label='kernel in noise region')
ax_e.axhline(0, color=C_GREY, lw=0.4)
ax_e.set_xlim(e_lo, e_hi)
ax_e.set_ylim(-1.0, 1.0)
ax_e.set_xlabel(r'Time ($\mu$s)', labelpad=2)
ax_e.set_ylabel('Amplitude', labelpad=4)
ax_e.annotate('high\n correlation',
              xy=(idx_aligned/fs*1e6 + 0.2, -0.55),
              xytext=(idx_aligned/fs*1e6 + 0.5, -0.85),
              fontsize=8, color=C_PULSE, weight='bold',
              arrowprops=dict(arrowstyle='-', color=C_PULSE, lw=0.6))
ax_e.annotate('low\n correlation',
              xy=(idx_misaligned/fs*1e6 + 0.2, -0.55),
              xytext=(idx_misaligned/fs*1e6 + 0.5, -0.85),
              fontsize=8, color='#666', weight='bold',
              arrowprops=dict(arrowstyle='-', color='#666', lw=0.6))
ax_e.legend(loc='upper right', fontsize=7.5, frameon=False)
# panel_label(ax_e, 'e', 'Cross-correlation as a sliding inner product',
#             'At each lag, multiply the kernel pointwise with the RF and sum. Peaks where the kernel matches.')

panel_label(ax_e, 'e', 'Cross-correlation as a sliding inner product',
            'At each lag, multiply the kernel pointwise with the RF and sum. '
            'Peaks where the kernel matches.',
            wrap=72)

# Right panel: equations
ax_f = fig.add_subplot(gs[2, 1])
ax_f.set_xlim(0, 1)
ax_f.set_ylim(0, 1)
ax_f.axis('off')
ax_f.text(0.5, 0.92, 'Per-sample matched-filter output:',
          ha='center', fontsize=9.5, weight='bold', color='#2C2C2A')
ax_f.text(0.5, 0.74,
          r'$y[n] = \sum_{k} \, h[k] \, \cdot \, x[n+k]$',
          ha='center', fontsize=12, color='#2C2C2A')
ax_f.text(0.5, 0.60,
          r'$h[k]$ = reference pulse (kernel)',
          ha='center', fontsize=9, color=C_PULSE)
ax_f.text(0.5, 0.50,
          r'$x[n]$ = raw RF trace',
          ha='center', fontsize=9, color=C_RAW)
ax_f.text(0.5, 0.40,
          r'$y[n]$ = correlation at lag $n$',
          ha='center', fontsize=9, color=C_MF)
ax_f.text(0.5, 0.22,
          'Implemented as convolution\nwith the time-reversed kernel:',
          ha='center', fontsize=8.5, color=C_GREY, style='italic')
ax_f.text(0.5, 0.08,
          r'$y[n] = (x * \tilde{h})[n], \quad \tilde{h}[k] = h[-k]$',
          ha='center', fontsize=10.5, color='#2C2C2A')
panel_label(ax_f, 'f', 'The operation, in equations', None)

# ------------------------------------------------------------------
# Row 4: matched-filter output (g) time, (h) frequency
# ------------------------------------------------------------------
ax_g = fig.add_subplot(gs[3, 0])
ax_g.fill_between(t_us[window],  env_mf[window], -env_mf[window],
                  color=C_FILL_MF, alpha=0.6, lw=0)
ax_g.plot(t_us[window], rf_mf_disp[window], color=C_MF, lw=0.7)
ax_g.axhline(0, color=C_GREY, lw=0.4)
# Mark the compressed-echo peak
ax_g.axvline(t_echo_us, color='#444', lw=0.5, ls=':')
ax_g.text(t_echo_us + 0.15, 1.05, 'compressed\nmain echo',
          fontsize=8, color='#444')
ax_g.set_xlim(t_lo, t_hi)
ax_g.set_ylim(-1.2, 1.4)
ax_g.set_xlabel(r'Time ($\mu$s)', labelpad=2)
ax_g.set_ylabel('Norm. amp.', labelpad=4)
panel_label(ax_g, 'g', 'Matched-filter output, time domain',
            'Echo concentrated into a tighter envelope; broadband noise floor suppressed.')

ax_h = fig.add_subplot(gs[3, 1])
ax_h.fill_between(f_mf, 0, X_mf_spec, color=C_MF, alpha=0.18, lw=0)
ax_h.plot(f_mf, X_mf_spec, color=C_MF, lw=1.2)
ax_h.axvline(fc/1e6, color=C_ACCENT, lw=0.7, ls='--')
ax_h.text(fc/1e6 + 0.4, 0.92, r'$f_c$', fontsize=9, color=C_ACCENT)
ax_h.set_xlim(0, 21)
ax_h.set_ylim(0, 1.08)
ax_h.set_xlabel('Frequency (MHz)', labelpad=2)
ax_h.set_ylabel('Norm. magnitude', labelpad=4)
# panel_label(ax_h, 'h', 'Matched-filter output, frequency',
#             'MF $=$ multiplication in freq.: $X(f)\\cdot H^*(f)$. Out-of-band noise suppressed.')

# panel_label(ax_h, 'h', 'Matched-filter output, frequency',
#             'MF $=$ multiplication in freq.: $X(f)\\cdot H^*(f)$. '
#             'Out-of-band noise suppressed.',
#             wrap=45)

panel_label(ax_h, 'h', 'Matched-filter output, frequency',
            'MF $=$ multiplication in freq.:\n'
            r'$X(f)\cdot H^*(f)$. Out-of-band noise suppressed.')



fig.savefig(os.path.join(OUT_DIR, 'figure_a2_matched_filter_v3.pdf'),
            bbox_inches='tight')
fig.savefig(os.path.join(OUT_DIR, 'figure_a2_matched_filter_v3.png'),
            dpi=300, bbox_inches='tight')
plt.close(fig)
print('Done.')
