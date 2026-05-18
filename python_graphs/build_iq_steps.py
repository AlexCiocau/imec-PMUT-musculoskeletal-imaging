# """
# IQ demodulation step-by-step figure using the REAL synthetic RF data
# from build_a2_a3.py (so the same toneburst echo, same noise, same
# matched filtering). Renders four panels stacked vertically:
#   Step 0: matched-filtered RF spectrum (real-valued, symmetric)
#   Step 1: complex local oscillator at -f_c (single delta)
#   Step 2: product spectrum (echo at DC + image at -2 f_c)
#   Step 3: lowpass-filtered baseband (only DC bump remains)
# Plus a fifth panel: time-domain I(t) and Q(t) for one specific time
# sample, with the (I, Q) point shown in the complex plane.
# """
# import os
# import numpy as np
# import matplotlib as mpl
# import matplotlib.pyplot as plt
# from scipy.signal import butter, filtfilt
# from matplotlib.patches import FancyArrowPatch

# mpl.rcParams.update({
#     'font.family': 'serif',
#     'font.serif': ['DejaVu Serif', 'Times New Roman', 'Times'],
#     'font.size': 10, 'axes.labelsize': 10, 'axes.titlesize': 10,
#     'xtick.labelsize': 9, 'ytick.labelsize': 9,
#     'axes.linewidth': 0.8,
#     'xtick.major.width': 0.8, 'ytick.major.width': 0.8,
#     'xtick.major.size': 3,    'ytick.major.size': 3,
#     'xtick.direction': 'in',  'ytick.direction': 'in',
#     'mathtext.fontset': 'cm',
#     'axes.spines.top': False, 'axes.spines.right': False,
#     'pdf.fonttype': 42, 'ps.fonttype': 42,
# })

# C_PURPLE  = '#534AB7'
# C_PURPLE_D= '#26215C'
# C_TEAL    = '#0F6E56'
# C_TEAL_D  = '#04342C'
# C_GREY    = '#5F5E5A'
# C_FILL    = '#EEEDFE'
# C_FILL_T  = '#E1F5EE'
# C_I       = '#185FA5'   # blue for I (real)
# C_Q       = '#A32D2D'   # red  for Q (imag)
# C_ACCENT  = '#0e5a73'

# # ============================================================
# # 1. RECREATE THE SYNTHETIC RF (must match build_a2_a3.py)
# # ============================================================
# fc = 10.5e6
# fs = 42e6
# T_total = 70e-6
# N = int(T_total * fs)
# t = np.arange(N) / fs

# T_pulse = 8 / (2 * fc)
# t_pulse = np.arange(0, T_pulse + 1/fs, 1/fs)
# duty = 0.67
# ref_pulse = (np.mod(t_pulse * fc, 1) < duty).astype(float) * np.sin(2*np.pi*fc*t_pulse)
# ref_pulse /= np.linalg.norm(ref_pulse)

# c_sound = 1480.0
# depth_mm = 33.8
# t_echo = 2 * (depth_mm * 1e-3) / c_sound
# i_echo = int(t_echo * fs)
# echo_decay = 1.5e-6
# n_echo_len = int(8e-6 * fs)
# t_echo_local = np.arange(n_echo_len) / fs

# # Echo with a small phase offset (in practice the echo phase is set by
# # the round-trip travel time and typically isn't 0 mod 2*pi at the
# # demodulator -- this gives a meaningful I/Q decomposition rather than
# # something that's almost purely real).
# phi0 = np.deg2rad(35)
# echo = np.exp(-t_echo_local / echo_decay) * np.sin(2*np.pi*fc*t_echo_local + phi0) * 0.55

# rng = np.random.default_rng(7)
# rf = np.zeros(N)
# rf[i_echo:i_echo + n_echo_len] += echo
# i_echo2 = i_echo + int(7e-6 * fs)
# echo2 = np.exp(-t_echo_local / echo_decay) * np.sin(2*np.pi*fc*t_echo_local) * 0.18
# rf[i_echo2:i_echo2 + n_echo_len] += echo2
# rf += rng.normal(0, 0.06, size=N)

# # matched filter
# rf_mf = np.convolve(rf, ref_pulse[::-1], mode='same')

# # segment around the echo (matches Figure A.3)
# seg = rf_mf[i_echo - 200 : i_echo + 1500]
# seg = seg / np.max(np.abs(seg))
# n_seg = len(seg)
# t_seg = np.arange(n_seg) / fs

# # IQ steps
# demod_carrier = np.exp(-1j * 2 * np.pi * fc * t_seg)
# product = seg * demod_carrier              # Step 2: before LPF
# b_lp, a_lp = butter(3, (fc * 0.8) / (fs/2), btype='low')
# iq = filtfilt(b_lp, a_lp, product)         # Step 3: final IQ baseband

# # ============================================================
# # 2. SPECTRA (two-sided, frequency in MHz, peak-normalised)
# # ============================================================
# def two_sided(x, fs):
#     X = np.fft.fftshift(np.fft.fft(x))
#     f = np.fft.fftshift(np.fft.fftfreq(len(x), 1/fs)) / 1e6
#     return f, np.abs(X) / np.max(np.abs(X))

# f_rf,   X_rf   = two_sided(seg,            fs)
# f_prod, X_prod = two_sided(product,        fs)
# f_iq,   X_iq   = two_sided(iq,             fs)

# # Local oscillator: a single complex exponential at -f_c, drawn as a
# # Kronecker-delta-like spectrum (a single nonzero bin at -f_c)
# lo = demod_carrier
# f_lo, X_lo_full = two_sided(lo, fs)
# # Smooth the LO display so it reads as a clear narrow spike rather than
# # numerical noise across the band
# X_lo = np.zeros_like(X_lo_full)
# i_neg_fc = int(np.argmin(np.abs(f_lo - (-fc/1e6))))
# # place a single spike at the bin closest to -fc
# X_lo[i_neg_fc] = 1.0

# # ============================================================
# # 3. PICK ONE TIME-DOMAIN SAMPLE TO DEMO I AND Q
# # ============================================================
# # Find a sample inside the echo region with non-trivial I AND Q (so the
# # complex-plane visualisation actually shows a tilted vector). Constraint:
# # |I| > 0.15 |IQ_max|, |Q| > 0.15 |IQ_max|, in the echo region.
# iq_max = np.max(np.abs(iq))
# i_search_lo = 200          # ~ start of segment with echo
# i_search_hi = 800          # well inside the echo
# search_slice = slice(i_search_lo, i_search_hi)
# abs_iq = np.abs(iq[search_slice])
# I_arr = np.real(iq[search_slice])
# Q_arr = np.imag(iq[search_slice])
# # score: prefer high magnitude AND high min(|I|, |Q|)
# score = np.minimum(np.abs(I_arr), np.abs(Q_arr)) * abs_iq
# i_demo = i_search_lo + int(np.argmax(score))
# I_val = float(np.real(iq[i_demo]))
# Q_val = float(np.imag(iq[i_demo]))
# t_demo_us = i_demo / fs * 1e6
# mag_val = np.hypot(I_val, Q_val)
# phase_val = np.degrees(np.arctan2(Q_val, I_val))

# print(f"Demo sample: i = {i_demo}, t = {t_demo_us:.3f} us within segment")
# print(f"  I = {I_val:+.4f},  Q = {Q_val:+.4f}")
# print(f"  |IQ| = {mag_val:.4f},  phase = {phase_val:+.1f} deg")

# # ============================================================
# # 4. BUILD FIGURE
# # ============================================================
# fig = plt.figure(figsize=(7.2, 11.2))
# gs = fig.add_gridspec(5, 2, height_ratios=[1, 1, 1, 1, 1.4],
#                       hspace=0.65, wspace=0.30,
#                       left=0.085, right=0.97, top=0.975, bottom=0.045)

# def setup_spectrum_axes(ax, xlim=(-25, 25)):
#     ax.axhline(0, color=C_GREY, lw=0.5)
#     ax.axvline(0, color=C_GREY, lw=0.4)
#     ax.set_xlim(*xlim)
#     ax.set_ylim(0, 1.1)
#     ax.set_yticks([])
#     ax.spines['left'].set_visible(False)
#     # mark f_c, -f_c
#     for f0, lbl in [(-fc/1e6, r'$-f_c$'), (fc/1e6, r'$+f_c$')]:
#         ax.axvline(f0, color=C_GREY, lw=0.3, ls=':')
#         ax.text(f0, -0.07, lbl, ha='center', va='top',
#                 fontsize=8.5, color=C_GREY,
#                 transform=ax.get_xaxis_transform())

# def panel_title(ax, step, title, subtitle):
#     ax.text(0.0, 1.18, f'Step {step} — {title}',
#             transform=ax.transAxes,
#             fontsize=11, weight='bold', color='#2C2C2A')
#     ax.text(0.0, 1.04, subtitle,
#             transform=ax.transAxes,
#             fontsize=9, color=C_GREY, style='italic')

# # ---- Step 0: RF spectrum (matched-filtered) -------------------
# ax0 = fig.add_subplot(gs[0, :])
# ax0.fill_between(f_rf, 0, X_rf, color=C_PURPLE, alpha=0.25, lw=0)
# ax0.plot(f_rf, X_rf, color=C_PURPLE_D, lw=0.9)
# setup_spectrum_axes(ax0)
# ax0.set_xlabel('Frequency (MHz)', labelpad=2)
# panel_title(ax0, 0, 'Received RF echo (after matched filter)',
#             'Real-valued signal — its spectrum is symmetric: a peak at $+f_c$ and its mirror at $-f_c$.')

# # ---- Step 1: Local oscillator -------------------------------
# ax1 = fig.add_subplot(gs[1, :])
# # Draw the LO as a stem at -fc
# df_bin = f_lo[1] - f_lo[0]
# spike_w = 0.25
# ax1.add_patch(plt.Rectangle((-fc/1e6 - spike_w/2, 0), spike_w, 1.0,
#                             color=C_TEAL_D, alpha=0.85))
# # A small label arrow
# ax1.annotate(r'$\delta(f + f_c)$',
#              xy=(-fc/1e6, 1.0), xytext=(-fc/1e6 + 4, 0.85),
#              fontsize=9, color=C_TEAL_D,
#              arrowprops=dict(arrowstyle='-', color=C_TEAL_D, lw=0.6))
# setup_spectrum_axes(ax1)
# ax1.set_xlabel('Frequency (MHz)', labelpad=2)
# panel_title(ax1, 1, 'Complex local oscillator: $e^{-j 2\\pi f_c t}$',
#             'A complex exponential is a single spike in the spectrum, located at $-f_c$.')

# # ---- Step 2: Product spectrum -------------------------------
# ax2 = fig.add_subplot(gs[2, :])
# ax2.fill_between(f_prod, 0, X_prod, color=C_PURPLE, alpha=0.25, lw=0)
# ax2.plot(f_prod, X_prod, color=C_PURPLE_D, lw=0.9)
# setup_spectrum_axes(ax2)
# ax2.set_xlabel('Frequency (MHz)', labelpad=2)
# # Annotate the two surviving copies
# ax2.annotate('echo at DC',
#              xy=(0, 0.85), xytext=(7, 1.0),
#              fontsize=9, color=C_PURPLE_D, weight='bold',
#              arrowprops=dict(arrowstyle='-', color=C_PURPLE_D, lw=0.6))
# ax2.annotate(r'image at $-2 f_c$',
#              xy=(-2*fc/1e6, 0.85), xytext=(-22, 1.0),
#              fontsize=9, color=C_PURPLE_D, weight='bold',
#              arrowprops=dict(arrowstyle='-', color=C_PURPLE_D, lw=0.6))
# # Show the LPF cutoff (only positive side, label it once below)
# fcut = fc * 0.8 / 1e6
# ax2.axvline( fcut, color=C_ACCENT, lw=0.7, ls='--', alpha=0.8)
# ax2.axvline(-fcut, color=C_ACCENT, lw=0.7, ls='--', alpha=0.8)
# ax2.text(fcut + 0.5, 0.52, 'LPF\ncutoff', fontsize=8, color=C_ACCENT,
#          ha='left', va='center')
# panel_title(ax2, 2, 'Product: signal × LO',
#             'Multiplication shifts the spectrum left by $f_c$. Result: echo at DC, image at $-2 f_c$.')

# # ---- Step 3: After LPF (final IQ) ---------------------------
# ax3 = fig.add_subplot(gs[3, :])
# ax3.fill_between(f_iq, 0, X_iq, color=C_PURPLE, alpha=0.25, lw=0)
# ax3.plot(f_iq, X_iq, color=C_PURPLE_D, lw=0.9)
# setup_spectrum_axes(ax3)
# ax3.set_xlabel('Frequency (MHz)', labelpad=2)
# ax3.annotate('echo, at baseband',
#              xy=(0, 0.95), xytext=(5, 0.95),
#              fontsize=9, color=C_PURPLE_D, weight='bold',
#              arrowprops=dict(arrowstyle='-', color=C_PURPLE_D, lw=0.6))
# panel_title(ax3, 3, 'IQ baseband signal (final result)',
#             r'Only the echo at DC remains. Spectrum is asymmetric — that is allowed because $iq(t) = I(t) + j\,Q(t)$ is complex.')

# # ---- Step 4: I(t), Q(t) for one sample ---------------------
# ax4a = fig.add_subplot(gs[4, 0])

# # Plot I(t) and Q(t) over a small region around the demo sample
# t_seg_us = t_seg * 1e6
# win_lo, win_hi = t_demo_us - 1.0, t_demo_us + 1.0
# m_win = (t_seg_us >= win_lo) & (t_seg_us <= win_hi)
# ax4a.plot(t_seg_us[m_win], np.real(iq)[m_win],
#           color=C_I, lw=1.2, label=r'$I(t) = \mathrm{Re}\{iq(t)\}$')
# ax4a.plot(t_seg_us[m_win], np.imag(iq)[m_win],
#           color=C_Q, lw=1.2, label=r'$Q(t) = \mathrm{Im}\{iq(t)\}$')
# # Mark the demo sample with vertical guide and dots
# ax4a.axvline(t_demo_us, color='#444', lw=0.6, ls='--', alpha=0.7)
# ax4a.plot(t_demo_us, I_val, 'o', color=C_I, ms=7,
#           mfc='white', mew=1.6, zorder=5)
# ax4a.plot(t_demo_us, Q_val, 'o', color=C_Q, ms=7,
#           mfc='white', mew=1.6, zorder=5)
# # annotate values
# ax4a.annotate(f'$I = {I_val:+.3f}$', xy=(t_demo_us, I_val),
#               xytext=(t_demo_us + 0.18, I_val - 0.10),
#               fontsize=8.5, color=C_I,
#               arrowprops=dict(arrowstyle='-', color=C_I, lw=0.5))
# ax4a.annotate(f'$Q = {Q_val:+.3f}$', xy=(t_demo_us, Q_val),
#               xytext=(t_demo_us - 0.55, Q_val - 0.18),
#               fontsize=8.5, color=C_Q,
#               arrowprops=dict(arrowstyle='-', color=C_Q, lw=0.5))
# ax4a.axhline(0, color=C_GREY, lw=0.4)
# ax4a.set_xlabel(r'Time within segment ($\mu$s)', labelpad=2)
# ax4a.set_ylabel('IQ amplitude', labelpad=4)
# ax4a.set_xlim(win_lo, win_hi)
# ax4a.legend(loc='lower right', fontsize=8, frameon=False)
# ax4a.text(0.0, 1.18, 'Step 4a — Time-domain I and Q',
#           transform=ax4a.transAxes,
#           fontsize=11, weight='bold', color='#2C2C2A')
# ax4a.text(0.0, 1.04, r'At each $t$, real and imaginary parts of $iq(t)$.',
#           transform=ax4a.transAxes,
#           fontsize=9, color=C_GREY, style='italic')

# # ---- Step 4b: Complex plane visualization -----------------
# ax4b = fig.add_subplot(gs[4, 1])
# # Draw I/Q axes
# lim = max(0.6, 1.25 * mag_val)
# ax4b.axhline(0, color=C_GREY, lw=0.6)
# ax4b.axvline(0, color=C_GREY, lw=0.6)
# ax4b.set_xlim(-lim, lim)
# ax4b.set_ylim(-lim, lim)
# ax4b.set_aspect('equal')
# ax4b.set_xticks([-0.4, 0, 0.4])
# ax4b.set_yticks([-0.4, 0, 0.4])
# ax4b.set_xlabel(r'$I$  (real)', labelpad=2, color=C_I)
# ax4b.set_ylabel(r'$Q$  (imag)', labelpad=4, color=C_Q)

# # Show the trajectory of iq(t) within the same window as a thin trace
# ax4b.plot(np.real(iq)[m_win], np.imag(iq)[m_win],
#           color='#bbbbbb', lw=0.8, alpha=0.7, zorder=1)
# # The demo point as a vector from origin
# ax4b.plot([0, I_val], [0, Q_val], color='#222', lw=1.0, zorder=3)
# ax4b.plot(I_val, Q_val, 'o', color='#222', ms=8,
#           mfc='white', mew=1.6, zorder=4)
# # Annotate the point
# ax4b.annotate(f'$({I_val:+.3f},\\, {Q_val:+.3f})$',
#               xy=(I_val, Q_val),
#               xytext=(I_val * 1.15 + 0.06, Q_val * 1.15 + 0.04),
#               fontsize=8.5, color='#222',
#               arrowprops=dict(arrowstyle='-', color='#222', lw=0.5))
# # Magnitude/phase note
# ax4b.text(0.02, 0.98,
#           f'$|iq| = {mag_val:.3f}$\n$\\angle iq = {phase_val:+.1f}^\\circ$',
#           transform=ax4b.transAxes,
#           fontsize=8.5, va='top', ha='left',
#           bbox=dict(boxstyle='round,pad=0.3', fc='white', ec=C_GREY, lw=0.5))
# ax4b.text(0.0, 1.18, 'Step 4b — In the complex plane',
#           transform=ax4b.transAxes,
#           fontsize=11, weight='bold', color='#2C2C2A')
# ax4b.text(0.0, 1.04, r'Same sample as a vector $iq = I + jQ$.',
#           transform=ax4b.transAxes,
#           fontsize=9, color=C_GREY, style='italic')

# fig.savefig('/mnt/user-data/outputs/figure_iq_steps.pdf', bbox_inches='tight')
# fig.savefig('/mnt/user-data/outputs/figure_iq_steps.png', dpi=300, bbox_inches='tight')
# plt.close(fig)
# print('Done.')



"""
IQ demodulation step-by-step figure using a real RF trace.

Same robustness logic as build_appendix_a2_a3.py: blanks the
saturating transmit transient at t<BLANK_END_US, then searches for
the actual target echo inside ECHO_SEARCH_LO_US..ECHO_SEARCH_HI_US.
"""
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt, hilbert
from scipy.io import loadmat

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
    'xtick.labelsize': 9, 'ytick.labelsize': 9,
    'axes.linewidth': 0.8,
    'xtick.major.width': 0.8, 'ytick.major.width': 0.8,
    'xtick.major.size': 3,    'ytick.major.size': 3,
    'xtick.direction': 'in',  'ytick.direction': 'in',
    'mathtext.fontset': 'cm',
    'axes.spines.top': False, 'axes.spines.right': False,
    'pdf.fonttype': 42, 'ps.fonttype': 42,
})

C_PURPLE  = '#534AB7'
C_PURPLE_D= '#26215C'
C_TEAL    = '#0F6E56'
C_TEAL_D  = '#04342C'
C_GREY    = '#5F5E5A'
C_FILL    = '#EEEDFE'
C_FILL_T  = '#E1F5EE'
C_I       = '#185FA5'
C_Q       = '#A32D2D'
C_ACCENT  = '#0e5a73'

fc = 10.5e6
fs = 42e6

# ------------------------------------------------------------------
# Load + blank transmit
# ------------------------------------------------------------------
rf = loadmat(RF_FILE)['trace'].squeeze().astype(float)
N = len(rf)
print(f"Loaded RF trace: {N} samples = {N/fs*1e6:.1f} us")

i_blank_end = int(BLANK_END_US * 1e-6 * fs)
rf_blank = rf.copy()
rf_blank[:i_blank_end] = 0.0
print(f"Blanked samples 0..{i_blank_end} (t < {BLANK_END_US} us)")

# ------------------------------------------------------------------
# Reference + matched filter
# ------------------------------------------------------------------
T_pulse = 8 / (2 * fc)
t_pulse = np.arange(0, T_pulse + 1/fs, 1/fs)
duty = 0.67
ref_pulse = (np.mod(t_pulse * fc, 1) < duty).astype(float) * np.sin(2*np.pi*fc*t_pulse)
ref_pulse /= np.linalg.norm(ref_pulse)

mf_kernel = ref_pulse[::-1]
rf_mf = np.convolve(rf_blank, mf_kernel, mode='same')
env_mf = np.abs(hilbert(rf_mf / np.max(np.abs(rf_mf))))

# ------------------------------------------------------------------
# Find echo inside the search window
# ------------------------------------------------------------------
i_search_lo = int(ECHO_SEARCH_LO_US * 1e-6 * fs)
i_search_hi = int(ECHO_SEARCH_HI_US * 1e-6 * fs)
i_echo = i_search_lo + int(np.argmax(env_mf[i_search_lo:i_search_hi]))
print(f"Detected echo at sample {i_echo} (t = {i_echo/fs*1e6:.2f} us)")

# Segment around the echo
n_pre, n_post = 200, 1500
i_lo = max(0, i_echo - n_pre)
i_hi = min(len(rf_mf), i_echo + n_post)
seg = rf_mf[i_lo:i_hi].astype(float)
seg_max = np.max(np.abs(seg))
if seg_max > 0:
    seg = seg / seg_max
n_seg = len(seg)
t_seg = np.arange(n_seg) / fs
print(f"IQ segment: {n_seg} samples")

# IQ steps
demod_carrier = np.exp(-1j * 2 * np.pi * fc * t_seg)
product = seg * demod_carrier
b_lp, a_lp = butter(3, (fc * 0.8) / (fs/2), btype='low')
iq = filtfilt(b_lp, a_lp, product)

# ------------------------------------------------------------------
# Spectra (two-sided, MHz)
# ------------------------------------------------------------------
def two_sided(x, fs):
    X = np.fft.fftshift(np.fft.fft(x))
    f = np.fft.fftshift(np.fft.fftfreq(len(x), 1/fs)) / 1e6
    return f, np.abs(X) / np.max(np.abs(X))

f_rf,   X_rf   = two_sided(seg,     fs)
f_prod, X_prod = two_sided(product, fs)
f_iq,   X_iq   = two_sided(iq,      fs)

# ------------------------------------------------------------------
# Pick a sample with non-trivial I AND Q in the strong-signal region
# ------------------------------------------------------------------
abs_iq = np.abs(iq)
threshold = 0.5 * abs_iq.max()
strong_mask = abs_iq > threshold
score = np.where(strong_mask,
                 np.minimum(np.abs(np.real(iq)), np.abs(np.imag(iq))) * abs_iq,
                 0)
i_demo = int(np.argmax(score))
I_val = float(np.real(iq[i_demo]))
Q_val = float(np.imag(iq[i_demo]))
t_demo_us = i_demo / fs * 1e6
mag_val = np.hypot(I_val, Q_val)
phase_val = np.degrees(np.arctan2(Q_val, I_val))
print(f"Demo sample: i = {i_demo}, t in segment = {t_demo_us:.3f} us")
print(f"  I = {I_val:+.4f},  Q = {Q_val:+.4f}")
print(f"  |IQ| = {mag_val:.4f},  phase = {phase_val:+.1f} deg")

# ------------------------------------------------------------------
# Build figure
# ------------------------------------------------------------------
fig = plt.figure(figsize=(7.2, 11.2))
gs = fig.add_gridspec(5, 2, height_ratios=[1, 1, 1, 1, 1.4],
                      hspace=0.65, wspace=0.30,
                      left=0.085, right=0.97, top=0.975, bottom=0.045)

def setup_spectrum_axes(ax, xlim=(-25, 25)):
    ax.axhline(0, color=C_GREY, lw=0.5)
    ax.axvline(0, color=C_GREY, lw=0.4)
    ax.set_xlim(*xlim)
    ax.set_ylim(0, 1.1)
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    for f0, lbl in [(-fc/1e6, r'$-f_c$'), (fc/1e6, r'$+f_c$')]:
        ax.axvline(f0, color=C_GREY, lw=0.3, ls=':')
        ax.text(f0, -0.07, lbl, ha='center', va='top',
                fontsize=8.5, color=C_GREY,
                transform=ax.get_xaxis_transform())

def panel_title(ax, step, title, subtitle):
    ax.text(0.0, 1.18, f'Step {step} — {title}',
            transform=ax.transAxes,
            fontsize=11, weight='bold', color='#2C2C2A')
    ax.text(0.0, 1.04, subtitle,
            transform=ax.transAxes,
            fontsize=9, color=C_GREY, style='italic')

# ---- Step 0: RF spectrum ----------------------------------------
ax0 = fig.add_subplot(gs[0, :])
ax0.fill_between(f_rf, 0, X_rf, color=C_PURPLE, alpha=0.25, lw=0)
ax0.plot(f_rf, X_rf, color=C_PURPLE_D, lw=0.9)
setup_spectrum_axes(ax0)
ax0.set_xlabel('Frequency (MHz)', labelpad=2)
panel_title(ax0, 0, 'Received RF echo (after matched filter)',
            'Real-valued signal — its spectrum is symmetric: a peak at $+f_c$ and its mirror at $-f_c$.')

# ---- Step 1: Local oscillator -----------------------------------
ax1 = fig.add_subplot(gs[1, :])
spike_w = 0.25
ax1.add_patch(plt.Rectangle((-fc/1e6 - spike_w/2, 0), spike_w, 1.0,
                            color=C_TEAL_D, alpha=0.85))
ax1.annotate(r'$\delta(f + f_c)$',
             xy=(-fc/1e6, 1.0), xytext=(-fc/1e6 + 4, 0.85),
             fontsize=9, color=C_TEAL_D,
             arrowprops=dict(arrowstyle='-', color=C_TEAL_D, lw=0.6))
setup_spectrum_axes(ax1)
ax1.set_xlabel('Frequency (MHz)', labelpad=2)
panel_title(ax1, 1, r'Complex local oscillator: $e^{-j 2\pi f_c t}$',
            'A complex exponential is a single spike in the spectrum, located at $-f_c$.')

# ---- Step 2: Product spectrum -----------------------------------
ax2 = fig.add_subplot(gs[2, :])
ax2.fill_between(f_prod, 0, X_prod, color=C_PURPLE, alpha=0.25, lw=0)
ax2.plot(f_prod, X_prod, color=C_PURPLE_D, lw=0.9)
setup_spectrum_axes(ax2)
ax2.set_xlabel('Frequency (MHz)', labelpad=2)
ax2.annotate('echo at DC',
             xy=(0, 0.85), xytext=(7, 1.0),
             fontsize=9, color=C_PURPLE_D, weight='bold',
             arrowprops=dict(arrowstyle='-', color=C_PURPLE_D, lw=0.6))
ax2.annotate(r'image at $-2 f_c$',
             xy=(-2*fc/1e6, 0.85), xytext=(-22, 1.0),
             fontsize=9, color=C_PURPLE_D, weight='bold',
             arrowprops=dict(arrowstyle='-', color=C_PURPLE_D, lw=0.6))
fcut = fc * 0.8 / 1e6
ax2.axvline( fcut, color=C_ACCENT, lw=0.7, ls='--', alpha=0.8)
ax2.axvline(-fcut, color=C_ACCENT, lw=0.7, ls='--', alpha=0.8)
ax2.text(fcut + 0.5, 0.52, 'LPF\ncutoff', fontsize=8, color=C_ACCENT,
         ha='left', va='center')
panel_title(ax2, 2, 'Product: signal × LO',
            'Multiplication shifts the spectrum left by $f_c$. Result: echo at DC, image at $-2 f_c$.')

# ---- Step 3: IQ baseband ----------------------------------------
ax3 = fig.add_subplot(gs[3, :])
ax3.fill_between(f_iq, 0, X_iq, color=C_PURPLE, alpha=0.25, lw=0)
ax3.plot(f_iq, X_iq, color=C_PURPLE_D, lw=0.9)
setup_spectrum_axes(ax3)
ax3.set_xlabel('Frequency (MHz)', labelpad=2)
ax3.annotate('echo, at baseband',
             xy=(0, 0.95), xytext=(5, 0.95),
             fontsize=9, color=C_PURPLE_D, weight='bold',
             arrowprops=dict(arrowstyle='-', color=C_PURPLE_D, lw=0.6))
panel_title(ax3, 3, 'IQ baseband signal (final result)',
            r'Only the echo at DC remains. Spectrum is asymmetric — that is allowed because $iq(t) = I(t) + j\,Q(t)$ is complex.')

# ---- Step 4a: Time-domain I and Q -------------------------------
ax4a = fig.add_subplot(gs[4, 0])
t_seg_us = t_seg * 1e6
win_lo, win_hi = t_demo_us - 1.0, t_demo_us + 1.0
m_win = (t_seg_us >= win_lo) & (t_seg_us <= win_hi)
ax4a.plot(t_seg_us[m_win], np.real(iq)[m_win],
          color=C_I, lw=1.2, label=r'$I(t) = \mathrm{Re}\{iq(t)\}$')
ax4a.plot(t_seg_us[m_win], np.imag(iq)[m_win],
          color=C_Q, lw=1.2, label=r'$Q(t) = \mathrm{Im}\{iq(t)\}$')
ax4a.axvline(t_demo_us, color='#444', lw=0.6, ls='--', alpha=0.7)
ax4a.plot(t_demo_us, I_val, 'o', color=C_I, ms=7,
          mfc='white', mew=1.6, zorder=5)
ax4a.plot(t_demo_us, Q_val, 'o', color=C_Q, ms=7,
          mfc='white', mew=1.6, zorder=5)
ax4a.annotate(f'$I = {I_val:+.3f}$', xy=(t_demo_us, I_val),
              xytext=(t_demo_us + 0.18, I_val - 0.10),
              fontsize=8.5, color=C_I,
              arrowprops=dict(arrowstyle='-', color=C_I, lw=0.5))
ax4a.annotate(f'$Q = {Q_val:+.3f}$', xy=(t_demo_us, Q_val),
              xytext=(t_demo_us - 0.55, Q_val - 0.18),
              fontsize=8.5, color=C_Q,
              arrowprops=dict(arrowstyle='-', color=C_Q, lw=0.5))
ax4a.axhline(0, color=C_GREY, lw=0.4)
ax4a.set_xlabel(r'Time within segment ($\mu$s)', labelpad=2)
ax4a.set_ylabel('IQ amplitude', labelpad=4)
ax4a.set_xlim(win_lo, win_hi)
ax4a.legend(loc='lower right', fontsize=8, frameon=False)
ax4a.text(0.0, 1.18, 'Step 4a — Time-domain I and Q',
          transform=ax4a.transAxes,
          fontsize=11, weight='bold', color='#2C2C2A')
ax4a.text(0.0, 1.04, r'At each $t$, real and imaginary parts of $iq(t)$.',
          transform=ax4a.transAxes,
          fontsize=9, color=C_GREY, style='italic')

# ---- Step 4b: Complex plane -------------------------------------
ax4b = fig.add_subplot(gs[4, 1])
lim = max(0.6, 1.25 * mag_val)
ax4b.axhline(0, color=C_GREY, lw=0.6)
ax4b.axvline(0, color=C_GREY, lw=0.6)
ax4b.set_xlim(-lim, lim)
ax4b.set_ylim(-lim, lim)
ax4b.set_aspect('equal')
ax4b.set_xticks([-0.4, 0, 0.4])
ax4b.set_yticks([-0.4, 0, 0.4])
ax4b.set_xlabel(r'$I$  (real)', labelpad=2, color=C_I)
ax4b.set_ylabel(r'$Q$  (imag)', labelpad=4, color=C_Q)
ax4b.plot(np.real(iq)[m_win], np.imag(iq)[m_win],
          color='#bbbbbb', lw=0.8, alpha=0.7, zorder=1)
ax4b.plot([0, I_val], [0, Q_val], color='#222', lw=1.0, zorder=3)
ax4b.plot(I_val, Q_val, 'o', color='#222', ms=8,
          mfc='white', mew=1.6, zorder=4)
ax4b.annotate(f'$({I_val:+.3f},\\, {Q_val:+.3f})$',
              xy=(I_val, Q_val),
              xytext=(I_val * 1.15 + 0.06, Q_val * 1.15 + 0.04),
              fontsize=8.5, color='#222',
              arrowprops=dict(arrowstyle='-', color='#222', lw=0.5))
ax4b.text(0.02, 0.98,
          f'$|iq| = {mag_val:.3f}$\n$\\angle iq = {phase_val:+.1f}^\\circ$',
          transform=ax4b.transAxes,
          fontsize=8.5, va='top', ha='left',
          bbox=dict(boxstyle='round,pad=0.3', fc='white', ec=C_GREY, lw=0.5))
ax4b.text(0.0, 1.18, 'Step 4b — In the complex plane',
          transform=ax4b.transAxes,
          fontsize=11, weight='bold', color='#2C2C2A')
ax4b.text(0.0, 1.04, r'Same sample as a vector $iq = I + jQ$.',
          transform=ax4b.transAxes,
          fontsize=9, color=C_GREY, style='italic')

fig.savefig(os.path.join(OUT_DIR, 'figure_iq_steps.pdf'), bbox_inches='tight')
fig.savefig(os.path.join(OUT_DIR, 'figure_iq_steps.png'), dpi=300, bbox_inches='tight')
plt.close(fig)
print('Wrote figure_iq_steps.{pdf,png}')