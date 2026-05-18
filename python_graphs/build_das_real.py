# """
# Stylised reproduction of the user's PowerPoint sketch, populated with
# matched-filtered RF from a real nylon-string acquisition.

# Layout:
#   - Probe along the top (drawn schematically), with seven highlighted
#     receive elements that we picked.
#   - Below each highlighted element, a vertical channel trace showing
#     the matched-filtered RF in a tight time window around the echo.
#   - A reflector dot below the centre channel.
#   - A few dashed horizontal lines indicating "phase markers" running
#     across all channels (matching the sketch).
# """
# import os
# import numpy as np
# import matplotlib as mpl
# import matplotlib.pyplot as plt
# from scipy.signal import hilbert

# mpl.rcParams.update({
#     'font.family': 'serif',
#     'font.serif': ['DejaVu Serif', 'Times New Roman', 'Times'],
#     'font.size': 10, 'axes.labelsize': 10,
#     'pdf.fonttype': 42, 'ps.fonttype': 42,
# })

# C_TRACE  = '#1f4e79'
# C_PIXEL  = '#a8323b'
# C_PROBE  = '#d97706'    # orange like in the sketch
# C_ELEM   = '#facc15'    # element bumps (yellow-ish)
# C_GUIDE  = '#a8323b'    # phase markers (red)
# C_GREY   = '#6b6b66'

# # ------------------------------------------------------------------
# # Load matched-filtered data
# # ------------------------------------------------------------------
# data = np.load('C:\Users\cioca100\Desktop\PythonGraph\NylonString_Rcv_35.mat')
# Frame_MF = data['Frame_MF']
# fs       = float(data['fs'])
# fc       = float(data['fc'])

# # Channels to display (1-indexed for readability; convert to 0-indexed)
# # 7 channels symmetric around channel 52 (the acoustic centre of the
# # target, which is offset from the probe centre because of the actual
# # target placement during the experiment).
# picks = [40, 44, 48, 52, 56, 60, 64]
# center_idx_in_picks = 3   # picks[3] == 52, the centre channel
# N_PICK = len(picks)

# # Time window around the echo
# T_CENTER_US = 44.0          # echo nominally at this time
# T_HALF_US   = 6           # +/- this many us
# t_lo_us = T_CENTER_US - T_HALF_US
# t_hi_us = T_CENTER_US + T_HALF_US
# i_lo = int(t_lo_us * 1e-6 * fs)
# i_hi = int(t_hi_us * 1e-6 * fs)
# n_samples_disp = i_hi - i_lo
# t_us = np.arange(i_lo, i_hi) / fs * 1e6
# print(f"Display window: {t_lo_us:.2f}-{t_hi_us:.2f} us "
#       f"({n_samples_disp} samples)")

# # Extract per-channel traces for the chosen window
# traces = np.array([Frame_MF[i_lo:i_hi, c-1] for c in picks])  # shape (N_PICK, n_samples)

# # Per-channel envelope (Hilbert) for marking the "echo arrival" dot
# envelopes = np.array([np.abs(hilbert(tr)) for tr in traces])
# peak_t_us = np.array([t_us[np.argmax(e)] for e in envelopes])

# # Normalise each trace to a common scale (use the global max so relative
# # amplitudes are preserved; the centre channel is genuinely the strongest)
# global_max = np.max(np.abs(traces))
# traces_norm = traces / global_max
# env_norm    = envelopes / global_max

# print(f"Normalised peak times per channel:")
# for i, c in enumerate(picks):
#     print(f"  ch {c}: t_peak = {peak_t_us[i]:.3f} us, "
#           f"|peak| = {np.max(env_norm[i]):.3f}")

# # ==================================================================
# # FIGURE
# # ==================================================================
# # We draw everything on ONE axes object, not multiple subplots, so the
# # probe-element <-> trace alignment is exact.

# fig, ax = plt.subplots(figsize=(9.2, 9.6))
# ax.set_aspect('auto')

# # Layout coordinates (figure-fraction-ish, but we use axes data units)
# # x-axis: spans the picked channels horizontally (use 1..7 then map to
# # their physical positions). Use the channel index (1..N_PICK) directly
# # as x because we want them evenly spaced like in the sketch.
# # y-axis: 0 = top of figure (probe sits here), increasing downward
# # represents time/depth in the channel traces.

# X_CH = np.arange(1, N_PICK + 1, dtype=float)  # x-position of each channel
# DX   = 1.0  # x-spacing between adjacent picked channels (in axes units)

# # Probe extent: spans the picked channels with a small margin on each side
# X_PROBE_L = X_CH[0]  - 0.5
# X_PROBE_R = X_CH[-1] + 0.5

# # Vertical layout
# Y_PROBE_TOP = 0.0
# Y_PROBE_BOT = -0.6     # probe height
# Y_TRACE_TOP = -1.4     # top of trace region
# Y_TRACE_BOT = -10.0    # bottom of trace region
# TRACE_AMP_DATA = 0.40  # max horizontal swing of trace, in axes-x units

# # ------------------------------------------------------------------
# # 1. Draw the probe (orange bar)
# # ------------------------------------------------------------------
# ax.fill_between(
#     [X_PROBE_L, X_PROBE_R], Y_PROBE_BOT, Y_PROBE_TOP,
#     color=C_PROBE, lw=0, zorder=2
# )
# # A subtle highlight band on top
# ax.fill_between(
#     [X_PROBE_L, X_PROBE_R], Y_PROBE_TOP - 0.08, Y_PROBE_TOP,
#     color='#ec9c2b', lw=0, zorder=3
# )
# # Probe label
# ax.text(
#     X_PROBE_L - 0.15, (Y_PROBE_TOP + Y_PROBE_BOT) / 2,
#     'PMUT\nelements', ha='right', va='center',
#     fontsize=11, weight='bold', color='#3f1d00'
# )

# # Channel-number labels above the probe
# for i, c in enumerate(picks):
#     ax.text(
#         X_CH[i], Y_PROBE_TOP + 0.18, f'ch {c}',
#         ha='center', va='bottom', fontsize=9, color='#3f1d00'
#     )

# # ------------------------------------------------------------------
# # 2. Element bumps under the probe (yellow semicircles, like in the sketch)
# # ------------------------------------------------------------------
# # Use ellipses so we control horizontal vs vertical scale independently
# from matplotlib.patches import Ellipse
# for x in X_CH:
#     e = Ellipse(
#         xy=(x, Y_PROBE_BOT - 0.15),
#         width=DX * 0.7, height=0.30,
#         facecolor=C_ELEM, edgecolor='#b4860c', lw=0.6, zorder=4
#     )
#     ax.add_patch(e)

# # ------------------------------------------------------------------
# # 3. Phase markers: dashed horizontal red lines (matching the sketch)
# # ------------------------------------------------------------------
# n_markers = 5
# y_marker_lo = Y_TRACE_TOP - 1.0   # below the top of the trace region
# y_marker_hi = Y_TRACE_BOT + 0.8
# y_markers = np.linspace(y_marker_lo, y_marker_hi, n_markers)
# for y in y_markers:
#     ax.plot([X_CH[0] - 0.4, X_CH[-1] + 0.4], [y, y],
#             ls='--', color=C_GUIDE, lw=0.9, alpha=0.55, zorder=2)

# # Annotate the phase markers (just one labeled, top-right)
# ax.annotate(
#     'phase reference\nmarkers', xy=(X_CH[-1] + 0.4, y_markers[0]),
#     xytext=(X_CH[-1] + 1.0, y_markers[0] - 0.6),
#     fontsize=8.5, color=C_GUIDE, style='italic', ha='left', va='center',
#     arrowprops=dict(arrowstyle='-', color=C_GUIDE, lw=0.5)
# )

# # ------------------------------------------------------------------
# # 4. Channel traces: vertical, centred under each element
# # ------------------------------------------------------------------
# # Each trace is plotted as a curve x = X_CH[i] + signal * scale, y = time-axis
# # We map t_us to y in the trace region: t_lo_us -> Y_TRACE_TOP, t_hi_us -> Y_TRACE_BOT
# y_of_t = lambda tt: Y_TRACE_TOP + (tt - t_lo_us) / (t_hi_us - t_lo_us) * \
#                                      (Y_TRACE_BOT - Y_TRACE_TOP)
# y_trace = y_of_t(t_us)

# # Time-of-flight curvature line: connect the per-channel peak times
# y_peak_per_ch = y_of_t(peak_t_us)

# # Plot each trace
# for i in range(N_PICK):
#     x_curve = X_CH[i] + traces_norm[i] * TRACE_AMP_DATA
#     ax.plot(x_curve, y_trace, color=C_TRACE, lw=0.9, zorder=5)

#     # Soft envelope fill behind the trace (gives the trace volume like the
#     # sketch's curves, where the bumps are visible)
#     env_x_lo = X_CH[i] - env_norm[i] * TRACE_AMP_DATA
#     env_x_hi = X_CH[i] + env_norm[i] * TRACE_AMP_DATA
#     ax.fill_betweenx(y_trace, env_x_lo, env_x_hi,
#                      color=C_TRACE, alpha=0.08, lw=0, zorder=4)

# # Draw the curvature line connecting peak times (dashed grey)
# ax.plot(X_CH, y_peak_per_ch, '--', color=C_PIXEL, lw=1.4,
#         alpha=0.85, zorder=6)
# # Markers on the curve at each peak
# for i in range(N_PICK):
#     ax.plot(X_CH[i], y_peak_per_ch[i], 'o',
#             color=C_PIXEL, ms=7, mfc='white', mew=2.0, zorder=7)

# # Annotate the curve
# ax.annotate(
#     'echo arrival times\n(later as channel moves\noff-axis)',
#     xy=(X_CH[-2], y_peak_per_ch[-2]),
#     xytext=(X_CH[-1] + 0.3, y_peak_per_ch[-1] + 1.5),
#     fontsize=9, color=C_PIXEL, style='italic', ha='left', va='top',
#     arrowprops=dict(arrowstyle='->', color=C_PIXEL, lw=0.7,
#                     connectionstyle='arc3,rad=0.2')
# )

# # ------------------------------------------------------------------
# # 5. Reflector dot below the centre channel
# # ------------------------------------------------------------------
# y_target = Y_TRACE_BOT - 0.6
# ax.plot(X_CH[center_idx_in_picks], y_target, 'o',
#         color=C_PIXEL, ms=14, zorder=10)
# ax.text(
#     X_CH[center_idx_in_picks], y_target + 0.5,
#     'nylon target', ha='center', va='top',
#     fontsize=10, weight='bold', color=C_PIXEL
# )

# # ------------------------------------------------------------------
# # 6. Time-axis annotation on the left
# # ------------------------------------------------------------------
# # Show a time scale next to the leftmost channel
# y_t_lo = y_of_t(t_lo_us)
# y_t_hi = y_of_t(t_hi_us)
# x_axis = X_PROBE_L - 0.05
# ax.annotate(
#     '', xy=(x_axis, y_t_hi), xytext=(x_axis, y_t_lo),
#     arrowprops=dict(arrowstyle='->', color=C_GREY, lw=0.8)
# )
# # (round-trip time label removed - axis is self-explanatory with the us markers)
# # Tick marks at known times
# for tt_us in [t_lo_us, T_CENTER_US, t_hi_us]:
#     yy = y_of_t(tt_us)
#     ax.plot([x_axis - 0.04, x_axis + 0.04], [yy, yy],
#             color=C_GREY, lw=0.7)
#     ax.text(x_axis - 0.07, yy, f'{tt_us:.2f} $\\mu$s',
#             ha='right', va='center', fontsize=7.5, color=C_GREY)

# # ------------------------------------------------------------------
# # Cosmetics
# # ------------------------------------------------------------------
# ax.set_xlim(X_PROBE_L - 1.4, X_PROBE_R + 1.7)
# ax.set_ylim(y_target - 0.7, Y_PROBE_TOP + 0.65)
# ax.set_xticks([])
# ax.set_yticks([])
# for sp in ax.spines.values():
#     sp.set_visible(False)

# ax.set_title(
#     'Per-channel echo from a nylon-string target (matched-filtered RF)',
#     fontsize=12, weight='bold', color='#2C2C2A', pad=22, loc='left'
# )
# fig.text(
#     0.04, 0.945,
#     'Each channel sees the same echo, but at a slightly different time '
#     'due to the different path\nlengths to the target. DAS shifts each '
#     'channel back into alignment.',
#     fontsize=9, color=C_GREY, style='italic', ha='left'
# )

# os.makedirs('/mnt/user-data/outputs', exist_ok=True)
# fig.savefig('/mnt/user-data/outputs/figure_a4_das_real.pdf',
#             bbox_inches='tight')
# fig.savefig('/mnt/user-data/outputs/figure_a4_das_real.png',
#             dpi=300, bbox_inches='tight')
# plt.close(fig)
# print("Done.")









"""
Stylised reproduction of the user's PowerPoint sketch, populated with
matched-filtered RF from a real nylon-string acquisition.

This script is self-contained: it loads NylonString_Rcv_35.mat,
applies the same averaging + channel reorder + matched filter as the
main MATLAB pipeline, then plots the result.
"""
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.signal import hilbert
from scipy.io import loadmat
from matplotlib.patches import Ellipse

# ============ USER CONFIG =====================================================
MAT_FILE = 'NylonString_Rcv_35.mat'   # path to the input .mat file
OUT_DIR  = 'out'                       # output folder

# Time window (in microseconds) shown in each channel trace.
# - T_CENTER_US: where the echo nominally sits in the matched-filtered RF
# - T_HALF_US:   how many microseconds before/after to show
# For example, T_CENTER_US=44, T_HALF_US=8 gives the window 36-52 us.
T_CENTER_US = 44.0
T_HALF_US   = 8

# Channels to display (1-indexed, 1..64). The figure picks 7 channels
# spread across the array; channel 52 is the acoustic centre of the
# target in this dataset (target is laterally offset from probe centre).
PICKS = [4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64]
CENTER_IDX_IN_PICKS = 3   # which one of PICKS is "the central channel"

# Matched-filter pulse parameters (must match the acquisition)
fs = 42e6
fc = 10.5e6
NUM_FRAMES_AVG = 30
# ==============================================================================

os.makedirs(OUT_DIR, exist_ok=True)

mpl.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['DejaVu Serif', 'Times New Roman', 'Times'],
    'font.size': 10, 'axes.labelsize': 10,
    'pdf.fonttype': 42, 'ps.fonttype': 42,
})

C_TRACE  = '#1f4e79'
C_PIXEL  = '#a8323b'
C_PROBE  = '#d97706'    # orange like in the sketch
C_ELEM   = '#facc15'    # yellow-ish element bumps
C_GUIDE  = '#a8323b'    # phase markers (red)
C_GREY   = '#6b6b66'

# ==============================================================================
# 1. LOAD + PROCESS THE RF DATA  (mirrors the main MATLAB pipeline)
# ==============================================================================
print(f"Loading {MAT_FILE}...")
m = loadmat(MAT_FILE)
RF = m['RcvData'][0, 0]   # (n_samples, 128, n_frames) int16
n_samples_total, n_total_ch, n_frames = RF.shape
print(f"  RF shape: {RF.shape}, dtype: {RF.dtype}")

# 1a. Take the receive aperture: MATLAB columns 65..128 -> Python 64..128
Frame_unavg = RF[:, 64:128, :].astype(np.float64)

# 1b. Frame averaging (use the same number as the main script, capped at
# the available frame count)
nfa = min(NUM_FRAMES_AVG, n_frames)
Frame = Frame_unavg[:, :, :nfa].mean(axis=2)   # (n_samples_total, 64)
print(f"  Averaged {nfa} frames")

# 1c. Channel reorder (interleaved -> physical sequential).
# MATLAB code (1-indexed):
#   for k = 1..32:
#     reorder_idx(2*k - 1) = 32 + k     (odd physical -> col 33,34,..,64)
#     reorder_idx(2*k)     = k          (even physical -> col 1,2,..,32)
# Python (0-indexed) equivalents:
reorder_idx = np.zeros(64, dtype=int)
for k in range(1, 33):
    reorder_idx[2*k - 2] = 32 + k - 1
    reorder_idx[2*k - 1] = k - 1
Frame = Frame[:, reorder_idx]
print("  Channel reorder applied")

# 1d. Matched filter: time-reversed transmit pulse
T_pulse = 8 / (2 * fc)
t_pulse = np.arange(0, T_pulse + 1/fs, 1/fs)
duty = 0.67
ref = (np.mod(t_pulse * fc, 1) < duty).astype(float) * np.sin(2*np.pi*fc*t_pulse)
ref /= np.linalg.norm(ref)
mf_kernel = ref[::-1]

n_samples, n_ch = Frame.shape
Frame_MF = np.zeros_like(Frame)
for c in range(n_ch):
    Frame_MF[:, c] = np.convolve(Frame[:, c], mf_kernel, mode='same')
print(f"  Matched filter applied to {n_ch} channels")

# ==============================================================================
# 2. EXTRACT THE PER-CHANNEL TRACES FOR THE CHOSEN WINDOW
# ==============================================================================
N_PICK = len(PICKS)
t_lo_us = T_CENTER_US - T_HALF_US
t_hi_us = T_CENTER_US + T_HALF_US
i_lo = int(t_lo_us * 1e-6 * fs)
i_hi = int(t_hi_us * 1e-6 * fs)
i_lo = max(0, i_lo)
i_hi = min(n_samples, i_hi)
n_samples_disp = i_hi - i_lo
t_us = np.arange(i_lo, i_hi) / fs * 1e6
print(f"\nDisplay window: {t_lo_us:.2f}-{t_hi_us:.2f} us "
      f"({n_samples_disp} samples)")

traces = np.array([Frame_MF[i_lo:i_hi, c-1] for c in PICKS])
envelopes = np.array([np.abs(hilbert(tr)) for tr in traces])
peak_t_us = np.array([t_us[np.argmax(e)] for e in envelopes])

global_max = np.max(np.abs(traces))
traces_norm = traces / global_max
env_norm    = envelopes / global_max

print("Per-channel echo arrival times:")
for i, c in enumerate(PICKS):
    print(f"  ch {c:2d}: t_peak = {peak_t_us[i]:.3f} us, "
          f"|peak| = {np.max(env_norm[i]):.3f}")

# ==============================================================================
# 3. BUILD THE FIGURE
# ==============================================================================
fig, ax = plt.subplots(figsize=(9.2, 9.6))
ax.set_aspect('auto')

# Layout coordinates (data units; x = channel index, y = time)
X_CH = np.arange(1, N_PICK + 1, dtype=float)
DX = 1.0
X_PROBE_L = X_CH[0]  - 0.5
X_PROBE_R = X_CH[-1] + 0.5

Y_PROBE_TOP = 0.0
Y_PROBE_BOT = -0.6
Y_TRACE_TOP = -1.4
Y_TRACE_BOT = -10.0
TRACE_AMP_DATA = 0.40

# ----- Probe -----
ax.fill_between(
    [X_PROBE_L, X_PROBE_R], Y_PROBE_BOT, Y_PROBE_TOP,
    color=C_PROBE, lw=0, zorder=2
)
ax.fill_between(
    [X_PROBE_L, X_PROBE_R], Y_PROBE_TOP - 0.08, Y_PROBE_TOP,
    color='#ec9c2b', lw=0, zorder=3
)
ax.text(
    X_PROBE_L - 0.15, (Y_PROBE_TOP + Y_PROBE_BOT) / 2,
    'PMUT\nelements', ha='right', va='center',
    fontsize=11, weight='bold', color='#3f1d00'
)

for i, c in enumerate(PICKS):
    ax.text(
        X_CH[i], Y_PROBE_TOP + 0.18, f'ch {c}',
        ha='center', va='bottom', fontsize=9, color='#3f1d00'
    )

# ----- Element ellipses -----
for x in X_CH:
    e = Ellipse(
        xy=(x, Y_PROBE_BOT - 0.15),
        width=DX * 0.7, height=0.30,
        facecolor=C_ELEM, edgecolor='#b4860c', lw=0.6, zorder=4
    )
    ax.add_patch(e)

# ----- Phase reference markers -----
n_markers = 5
y_marker_lo = Y_TRACE_TOP - 1.0
y_marker_hi = Y_TRACE_BOT + 0.8
y_markers = np.linspace(y_marker_lo, y_marker_hi, n_markers)
for y in y_markers:
    ax.plot([X_CH[0] - 0.4, X_CH[-1] + 0.4], [y, y],
            ls='--', color=C_GUIDE, lw=0.9, alpha=0.55, zorder=2)

ax.annotate(
    'phase reference\nmarkers', xy=(X_CH[-1] + 0.4, y_markers[0]),
    xytext=(X_CH[-1] + 1.0, y_markers[0] - 0.6),
    fontsize=8.5, color=C_GUIDE, style='italic', ha='left', va='center',
    arrowprops=dict(arrowstyle='-', color=C_GUIDE, lw=0.5)
)

# ----- Channel traces (vertical) -----
y_of_t = lambda tt: Y_TRACE_TOP + (tt - t_lo_us) / (t_hi_us - t_lo_us) * \
                                     (Y_TRACE_BOT - Y_TRACE_TOP)
y_trace = y_of_t(t_us)
y_peak_per_ch = y_of_t(peak_t_us)

for i in range(N_PICK):
    x_curve = X_CH[i] + traces_norm[i] * TRACE_AMP_DATA
    ax.plot(x_curve, y_trace, color=C_TRACE, lw=0.9, zorder=5)
    env_x_lo = X_CH[i] - env_norm[i] * TRACE_AMP_DATA
    env_x_hi = X_CH[i] + env_norm[i] * TRACE_AMP_DATA
    ax.fill_betweenx(y_trace, env_x_lo, env_x_hi,
                     color=C_TRACE, alpha=0.08, lw=0, zorder=4)

# Curvature line connecting per-channel peak times
ax.plot(X_CH, y_peak_per_ch, '--', color=C_PIXEL, lw=1.4,
        alpha=0.85, zorder=6)
for i in range(N_PICK):
    ax.plot(X_CH[i], y_peak_per_ch[i], 'o',
            color=C_PIXEL, ms=7, mfc='white', mew=2.0, zorder=7)

ax.annotate(
    'echo arrival times\n(later as channel moves\noff-axis)',
    xy=(X_CH[-2], y_peak_per_ch[-2]),
    xytext=(X_CH[-1] + 0.3, y_peak_per_ch[-1] + 1.5),
    fontsize=9, color=C_PIXEL, style='italic', ha='left', va='top',
    arrowprops=dict(arrowstyle='->', color=C_PIXEL, lw=0.7,
                    connectionstyle='arc3,rad=0.2')
)

# ----- Reflector dot below the centre channel -----
y_target = Y_TRACE_BOT - 0.6
ax.plot(X_CH[CENTER_IDX_IN_PICKS], y_target, 'o',
        color=C_PIXEL, ms=14, zorder=10)
ax.text(
    X_CH[CENTER_IDX_IN_PICKS], y_target + 0.5,
    'nylon target', ha='center', va='top',
    fontsize=10, weight='bold', color=C_PIXEL
)

# ----- Time axis on the left -----
y_t_lo = y_of_t(t_lo_us)
y_t_hi = y_of_t(t_hi_us)
x_axis = X_PROBE_L - 0.05
ax.annotate(
    '', xy=(x_axis, y_t_hi), xytext=(x_axis, y_t_lo),
    arrowprops=dict(arrowstyle='->', color=C_GREY, lw=0.8)
)
# Tick marks: at t_lo, t_center, t_hi
for tt_us in [t_lo_us, T_CENTER_US, t_hi_us]:
    yy = y_of_t(tt_us)
    ax.plot([x_axis - 0.04, x_axis + 0.04], [yy, yy],
            color=C_GREY, lw=0.7)
    ax.text(x_axis - 0.07, yy, f'{tt_us:.2f} $\\mu$s',
            ha='right', va='center', fontsize=7.5, color=C_GREY)

# ----- Cosmetics -----
ax.set_xlim(X_PROBE_L - 1.4, X_PROBE_R + 1.7)
ax.set_ylim(y_target - 0.7, Y_PROBE_TOP + 0.65)
ax.set_xticks([])
ax.set_yticks([])
for sp in ax.spines.values():
    sp.set_visible(False)

ax.set_title(
    'Per-channel echo from a nylon-string target (matched-filtered RF)',
    fontsize=12, weight='bold', color='#2C2C2A', pad=22, loc='left'
)
fig.text(
    0.04, 0.945,
    'Each channel sees the same echo, but at a slightly different time '
    'due to the different path\nlengths to the target. DAS shifts each '
    'channel back into alignment.',
    fontsize=9, color=C_GREY, style='italic', ha='left'
)

# ----- Save -----
fig.savefig(os.path.join(OUT_DIR, 'figure_a4_das_real.pdf'),
            bbox_inches='tight')
fig.savefig(os.path.join(OUT_DIR, 'figure_a4_das_real.png'),
            dpi=300, bbox_inches='tight')
plt.close(fig)
print(f"\nWrote figure_a4_das_real.{{pdf,png}} to {OUT_DIR}/")