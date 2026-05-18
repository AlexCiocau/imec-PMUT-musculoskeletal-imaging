"""
Figure: "Delay-and-Sum reconstruction" diagrammatic explanation,
in the style of the user's PowerPoint sketch.

Layout (left-to-right):
  LEFT:   Image grid showing pixel states (done/current/todo) with a
          red distance arrow from the probe to the current pixel.
  MIDDLE: Caption + horizontal connecting arrow.
  RIGHT:  Probe + N channels with synthetic RF traces hanging below.
          A red dot on each trace marks the sample picked at the
          per-channel time-of-flight from the current pixel.
          Below: convergence lines from each picked dot into a Sigma
          symbol -> the resulting pixel value.

Synthetic data: each channel's trace is a Gaussian-modulated cosine
centred at tau_i, the round-trip TOF from the chosen pixel to that
channel. This is illustrative; the figure shows the *operation*,
not a real acquisition.
"""
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Ellipse

mpl.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['DejaVu Serif', 'Times New Roman', 'Times'],
    'font.size': 10, 'axes.labelsize': 10,
    'pdf.fonttype': 42, 'ps.fonttype': 42,
})

# Palette
C_DONE    = '#2563eb'
C_CURR    = '#dc2626'
C_CURR_F  = '#fee2e2'
C_TODO    = '#d4d4d8'
C_TODO_E  = '#a1a1aa'
C_TRACE   = '#1e40af'
C_PROBE   = '#d97706'
C_ELEM    = '#facc15'
C_PICK    = '#dc2626'
C_GREY    = '#6b6b66'

# ==================================================================
# IMAGE GRID
# ==================================================================
N_COLS = 6
N_ROWS = 8
CURR_COL, CURR_ROW = 1, 5
state = np.full((N_ROWS, N_COLS), 2, dtype=int)
for c in range(CURR_COL):
    state[:, c] = 0
for r in range(CURR_ROW):
    state[r, CURR_COL] = 0
state[CURR_ROW, CURR_COL] = 1

# ==================================================================
# CHANNELS + SYNTHETIC TRACES
# ==================================================================
N_CH = 7
x_pix = 1.5
z_pix = 2.5
ch_x = np.arange(N_CH, dtype=float)
PSEUDO_C = 1.5
tau = (z_pix + np.sqrt((ch_x - x_pix)**2 + z_pix**2)) / PSEUDO_C
t0 = tau.min() - 1.5
tau = tau - t0

n_samples = 220
t_samples = np.linspace(0, 12, n_samples)
rng = np.random.default_rng(3)
def make_trace(t_peak):
    sigma = 0.45
    fc = 1.6
    sig = np.exp(-((t_samples - t_peak)**2) / (2 * sigma**2)) \
          * np.cos(2*np.pi*fc*(t_samples - t_peak))
    return sig + rng.normal(0, 0.04, size=n_samples)

traces = np.array([make_trace(tau[i]) for i in range(N_CH)])
traces /= np.max(np.abs(traces))
picked_idx = [int(np.argmin(np.abs(t_samples - tau[i])))
              for i in range(N_CH)]
picked_t = [t_samples[picked_idx[i]] for i in range(N_CH)]

# ==================================================================
# FIGURE
# ==================================================================
fig = plt.figure(figsize=(13.5, 8.6))
gs = fig.add_gridspec(
    1, 3, width_ratios=[1.0, 0.30, 1.4],
    wspace=0.05,
    left=0.04, right=0.99, top=0.93, bottom=0.05,
)

# ------------------------------------------------------------------
# LEFT: image grid
# ------------------------------------------------------------------
ax_grid = fig.add_subplot(gs[0, 0])
ax_grid.set_aspect('equal')

GRID_X0, GRID_Y0 = 0.0, 0.0
CELL = 1.0
GRID_W = N_COLS * CELL
GRID_H = N_ROWS * CELL

ax_grid.annotate(
    '', xy=(GRID_X0 + GRID_W, GRID_Y0 - 0.3),
        xytext=(GRID_X0,         GRID_Y0 - 0.3),
    arrowprops=dict(arrowstyle='->', color=C_DONE, lw=1.4),
)
ax_grid.text(
    GRID_X0 + GRID_W / 2, GRID_Y0 - 0.7,
    'Reconstructed image\nwidth (x-axis)',
    ha='center', va='bottom', fontsize=11, color=C_DONE, weight='bold',
)
ax_grid.annotate(
    '', xy=(GRID_X0 - 0.3, GRID_Y0 + GRID_H),
        xytext=(GRID_X0 - 0.3, GRID_Y0),
    arrowprops=dict(arrowstyle='->', color=C_DONE, lw=1.4),
)
ax_grid.text(
    GRID_X0 - 0.6, GRID_Y0 + GRID_H / 2,
    'Depth (z-axis)', rotation=90,
    ha='right', va='center', fontsize=11, color=C_DONE, weight='bold',
)

for r in range(N_ROWS):
    for c in range(N_COLS):
        x = GRID_X0 + c * CELL
        y = GRID_Y0 + r * CELL
        s = state[r, c]
        if s == 0:
            face, edge, lw, ls = C_DONE, C_DONE, 0.8, '-'
        elif s == 1:
            face, edge, lw, ls = 'white', C_CURR, 1.8, '--'
        else:
            face, edge, lw, ls = C_TODO, C_TODO_E, 0.6, '-'
        ax_grid.add_patch(Rectangle(
            (x, y), CELL * 0.95, CELL * 0.95,
            facecolor=face, edgecolor=edge, lw=lw, linestyle=ls,
        ))

band_x = GRID_X0 + CURR_COL * CELL
band_w = CELL * 0.95
band_h = (CURR_ROW + 0.95) * CELL
ax_grid.add_patch(Rectangle(
    (band_x, GRID_Y0), band_w, band_h,
    facecolor=C_CURR_F, edgecolor='none', alpha=0.45, zorder=0,
))

arrow_x = band_x + band_w / 2
arrow_y_top = GRID_Y0 + 0.05
arrow_y_bot = GRID_Y0 + (CURR_ROW + 0.45) * CELL
ax_grid.annotate(
    '', xy=(arrow_x, arrow_y_bot),
        xytext=(arrow_x, arrow_y_top),
    arrowprops=dict(arrowstyle='->', color=C_CURR, lw=1.8,
                    mutation_scale=18),
)
ax_grid.text(
    arrow_x + 0.10, (arrow_y_top + arrow_y_bot) / 2,
    'Distance to reconstructed pixel',
    rotation=-90, ha='left', va='center',
    fontsize=10, weight='bold', color=C_CURR,
)

LEG_Y = GRID_H + 1.0
LEG_X = GRID_X0
items = [
    (C_DONE,  '-',  'Reconstructed pixels'),
    (C_TODO,  '-',  'Yet to be reconstructed'),
    ('white', '--', 'Current pixel'),
]
for k, (face, ls, lbl) in enumerate(items):
    yy = LEG_Y + k * 0.7
    edge = (C_DONE if face == C_DONE
            else C_TODO_E if face == C_TODO else C_CURR)
    lw = 1.8 if face == 'white' else 0.8
    ax_grid.add_patch(Rectangle(
        (LEG_X, yy), 0.55, 0.5,
        facecolor=face, edgecolor=edge, lw=lw, linestyle=ls,
    ))
    ax_grid.text(LEG_X + 0.75, yy + 0.25, lbl,
                 va='center', ha='left', fontsize=10)

ax_grid.set_xlim(GRID_X0 - 1.4, GRID_X0 + GRID_W + 0.4)
ax_grid.set_ylim(GRID_Y0 + GRID_H + 2.6, GRID_Y0 - 1.4)
ax_grid.set_xticks([])
ax_grid.set_yticks([])
for sp in ax_grid.spines.values():
    sp.set_visible(False)

# ------------------------------------------------------------------
# MIDDLE
# ------------------------------------------------------------------
ax_mid = fig.add_subplot(gs[0, 1])
ax_mid.axis('off')
ax_mid.set_xlim(0, 1)
ax_mid.set_ylim(0, 1)
ax_mid.annotate(
    '', xy=(0.95, 0.5), xytext=(0.05, 0.5),
    arrowprops=dict(arrowstyle='->', color=C_GREY, lw=2.2,
                    mutation_scale=22),
)
ax_mid.text(
    0.5, 0.70,
    'For each channel,\ncompute the time of flight\nto the pixel and pick the\nsample at that time.',
    ha='center', va='center', fontsize=10, color='#2C2C2A',
)

# ------------------------------------------------------------------
# RIGHT
# ------------------------------------------------------------------
ax_ch = fig.add_subplot(gs[0, 2])
X_CH_VALS = np.arange(N_CH, dtype=float)
X_CH_L = -0.5
X_CH_R = (N_CH - 1) + 0.5
Y_PROBE_TOP = 0.0
Y_PROBE_BOT = -0.55
Y_TRACE_TOP = -1.3
Y_TRACE_BOT = -8.5
TRACE_AMP = 0.38

ax_ch.add_patch(Rectangle(
    (X_CH_L, Y_PROBE_BOT), X_CH_R - X_CH_L,
    Y_PROBE_TOP - Y_PROBE_BOT,
    facecolor=C_PROBE, edgecolor='none', zorder=2,
))
ax_ch.add_patch(Rectangle(
    (X_CH_L, Y_PROBE_TOP - 0.08), X_CH_R - X_CH_L, 0.08,
    facecolor='#ec9c2b', edgecolor='none', zorder=3,
))
ax_ch.text(
    X_CH_L - 0.20, (Y_PROBE_TOP + Y_PROBE_BOT) / 2,
    'PMUT\ntransducer\ncells',
    ha='right', va='center', fontsize=10, weight='bold', color='#3f1d00',
)
ax_ch.text(X_CH_VALS[0],  Y_PROBE_TOP + 0.15, '1',
           ha='center', va='bottom', fontsize=10, color='#3f1d00')
ax_ch.text(X_CH_VALS[-1], Y_PROBE_TOP + 0.15, 'N',
           ha='center', va='bottom', fontsize=10, color='#3f1d00')
ax_ch.text(
    X_CH_R + 0.05, Y_PROBE_TOP + 0.65,
    'RF data per channel',
    ha='right', va='center', fontsize=10, weight='bold', color='#2C2C2A',
)

for x in X_CH_VALS:
    ax_ch.add_patch(Ellipse(
        xy=(x, Y_PROBE_BOT - 0.12),
        width=0.65, height=0.26,
        facecolor=C_ELEM, edgecolor='#b4860c', lw=0.5, zorder=4,
    ))

y_of_t = lambda tt: Y_TRACE_TOP + (tt - t_samples[0]) / \
                                  (t_samples[-1] - t_samples[0]) * \
                                  (Y_TRACE_BOT - Y_TRACE_TOP)
y_trace = y_of_t(t_samples)

for i in range(N_CH):
    x_curve = X_CH_VALS[i] + traces[i] * TRACE_AMP
    ax_ch.plot(x_curve, y_trace, color=C_TRACE, lw=0.9, zorder=5)
    env = np.abs(traces[i])
    env_x_lo = X_CH_VALS[i] - env * TRACE_AMP
    env_x_hi = X_CH_VALS[i] + env * TRACE_AMP
    ax_ch.fill_betweenx(y_trace, env_x_lo, env_x_hi,
                        color=C_TRACE, alpha=0.06, lw=0, zorder=4)

xs_dots = [X_CH_VALS[i] + traces[i, picked_idx[i]] * TRACE_AMP
           for i in range(N_CH)]
ys_dots = [y_of_t(picked_t[i]) for i in range(N_CH)]
for i in range(N_CH):
    ax_ch.plot(xs_dots[i], ys_dots[i], 'o',
               color=C_PICK, ms=10, zorder=10,
               mec='white', mew=1.5)

ax_ch.annotate(
    'sample picked at\n$\\tau_i = t_{TX}+t_{RX,i}$',
    xy=(xs_dots[-1], ys_dots[-1]),
    xytext=(X_CH_VALS[-1] + 1.0, ys_dots[-1] - 1.6),
    fontsize=9.5, color=C_PICK, weight='bold',
    ha='left', va='center',
    arrowprops=dict(arrowstyle='->', color=C_PICK, lw=0.7,
                    connectionstyle='arc3,rad=-0.25'),
)

Y_OUT = Y_TRACE_BOT - 1.5
X_OUT = (X_CH_L + X_CH_R) / 2
sigma_box_w, sigma_box_h = 1.0, 0.7
ax_ch.add_patch(Rectangle(
    (X_OUT - sigma_box_w/2, Y_OUT - sigma_box_h/2),
    sigma_box_w, sigma_box_h,
    facecolor=C_CURR_F, edgecolor=C_CURR, lw=1.2, zorder=6,
))
ax_ch.text(X_OUT, Y_OUT, r'$\sum$',
           ha='center', va='center',
           fontsize=18, color=C_CURR, weight='bold', zorder=7)
for i in range(N_CH):
    ax_ch.plot(
        [xs_dots[i], X_OUT],
        [ys_dots[i], Y_OUT + sigma_box_h / 2],
        ls='-', color=C_PICK, lw=0.5, alpha=0.30, zorder=5,
    )
ax_ch.text(
    X_OUT + sigma_box_w / 2 + 0.25, Y_OUT,
    r'$\to\;\;$pixel value $= \sum_{i=1}^{N} s_i(\tau_i)$',
    ha='left', va='center',
    fontsize=12, color=C_CURR, weight='bold',
)

ax_ch.set_xlim(X_CH_L - 1.4, X_CH_R + 2.5)
ax_ch.set_ylim(Y_OUT - 1.3, Y_PROBE_TOP + 1.0)
ax_ch.set_xticks([])
ax_ch.set_yticks([])
for sp in ax_ch.spines.values():
    sp.set_visible(False)

fig.suptitle(
    'Delay-and-sum reconstruction: per-pixel sample selection and summation',
    fontsize=12.5, weight='bold', color='#2C2C2A', y=0.97,
)

os.makedirs('out', exist_ok=True)
fig.savefig('out/figure_a4_das_reconstruction.pdf', bbox_inches='tight')
fig.savefig('out/figure_a4_das_reconstruction.png',
            dpi=300, bbox_inches='tight')
plt.close(fig)
print('Done.')
