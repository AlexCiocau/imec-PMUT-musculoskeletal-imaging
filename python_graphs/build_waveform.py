"""
waveform_figure.py
------------------
Reproduces the RAW waveforms figure from waveform.fig.

Requirements:
    pip install scipy numpy matplotlib

Usage:
    python waveform_figure.py
"""

import scipy.io as sio
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os

C_RAW   = '#1f4e79'
C_PULSE = '#2f6f4f'
C_MF    = '#a8323b'
C_FILL_RAW = '#cdd5dc'
C_FILL_MF  = '#e6cdd1'
C_FILL_PUL = '#cfe0d4'
C_GREY  = '#5F5E5A'
C_ACCENT = '#0e5a73'
C_DARKGREY  = "#20201E"

# ── 1. Load data from the .fig file ──────────────────────────────────────────
FIG_FILE = 'waveform.fig'   # <-- change path if needed

data = sio.loadmat(FIG_FILE, squeeze_me=True, struct_as_record=False)
hgs  = data['hgS_070000']
axes = hgs.children          # 3 axes objects

def get_line(ax):
    """Return the first lineseries child of an axes as a dict."""
    for ch in np.atleast_1d(ax.children):
        if hasattr(ch, 'type') and 'lineseries' in str(ch.type):
            p = ch.properties
            return {
                'x':     p.XData,
                'y':     p.YData,
                'color': np.atleast_1d(p.Color).astype(float),
                'lw':    float(p.LineWidth),
            }
    return None

l0 = get_line(axes[0])   # TX Driving signal  – black
l1 = get_line(axes[1])   # RX 2-Layer PCB     – blue
l2 = get_line(axes[2])   # RX 4-Layer PCB     – red


# ── 2. Style ──────────────────────────────────────────────────────────────────
mpl.rcParams.update({
    'font.family':       'serif',
    'font.serif':        ['DejaVu Serif', 'Times New Roman', 'Times'],
    'font.size':          10,
    'axes.labelsize':     10,
    'axes.titlesize':     10,
    'xtick.labelsize':    8.5,
    'ytick.labelsize':    8.5,
    'axes.linewidth':     0.8,
    'xtick.major.width':  0.8,
    'ytick.major.width':  0.8,
    'xtick.major.size':   3,
    'ytick.major.size':   3,
    'xtick.direction':   'in',
    'ytick.direction':   'in',
    'mathtext.fontset':  'cm',
    'pdf.fonttype':       42,    # embeds fonts properly in PDF
    'ps.fonttype':        42,
    'axes.spines.top':    True,  # keep all four spines (box=on)
    'axes.spines.right':  True,
})


# ── 3. Layout ─────────────────────────────────────────────────────────────────
fig, axs = plt.subplots(
    3, 1,
    figsize=(8, 7.5),                        # width, height in inches
    gridspec_kw={'hspace': 0.48},            # vertical spacing between subplots
)
# fig.suptitle('RAW waveforms', fontsize=11, weight='bold', color='#2C2C2A', y=0.98)

C_ZERO = '#5F5E5A'    # colour of the y=0 reference line
GRID_KW = dict(color='#cccccc', linewidth=0.5, linestyle='--', zorder=0)

TITLE_KW  = dict(fontsize=10, weight='bold', color='#2C2C2A', pad=4)
YLABEL_KW = dict(labelpad=4)


# # ── 4a. Subplot 1 – TX Driving signal ────────────────────────────────────────
# ax = axs[0]
# ax.plot(l0['x'], l0['y'] / 1000,
#         color=C_DARKGREY, lw=1.1, zorder=3)

# ax.set_title('Oscilloscope measurement of TX Driving signal', **TITLE_KW)
# ax.set_xlim(-0.2, 1.0)
# ax.set_ylim(-14, 14)
# ax.set_xticks([-0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0])
# ax.set_yticks([-10, -5, 0, 5, 10])
# ax.set_ylabel('Amplitude (V)', **YLABEL_KW)
# ax.grid(True, **GRID_KW)
# ax.set_axisbelow(True)
# ax.axhline(0, color=C_ZERO, lw=0.4, zorder=1)


# # ── 4b. Subplot 2 – RX 2-Layer PCB ───────────────────────────────────────────
# ax = axs[1]
# ax.plot(l1['x'], l1['y'],
#         color=C_RAW, lw=1.1, zorder=3)

# ax.set_title('Hydrophone reading of pulse transmitted using 2-Layer PCB', **TITLE_KW)
# ax.set_xlim(37, 39.5)
# ax.set_ylim(-2, 2)
# ax.set_xticks([37, 37.5, 38, 38.5, 39, 39.5])
# ax.set_yticks([-2, -1, 0, 1, 2])
# ax.set_ylabel('Amplitude (mV)', **YLABEL_KW)
# ax.grid(True, **GRID_KW)
# ax.set_axisbelow(True)
# ax.axhline(0, color=C_ZERO, lw=0.4, zorder=1)


# # ── 4c. Subplot 3 – RX 4-Layer PCB ───────────────────────────────────────────
# ax = axs[2]
# ax.plot(l2['x'], l2['y'],
#         color=C_MF, lw=1.1, zorder=3)

# ax.set_title('Hydrophone reading of pulse transmitted using 4-Layer PCB', **TITLE_KW)
# ax.set_xlim(37, 39.5)
# ax.set_ylim(-2, 2)
# ax.set_xticks([37, 37.5, 38, 38.5, 39, 39.5])
# ax.set_yticks([-2, -1, 0, 1, 2])
# ax.set_ylabel('Amplitude (mV)', **YLABEL_KW)
# ax.set_xlabel(r'Time ($\mu$s)', labelpad=3)
# ax.grid(True, **GRID_KW)
# ax.set_axisbelow(True)
# ax.axhline(0, color=C_ZERO, lw=0.4, zorder=1)

# ------------------------------------------------------------------------------

# ── 3.5 Helper Function for Vpp Markers ──────────────────────────────────────
def add_vpp_markers(ax, y_data, unit='V'):
    v_plus = np.max(y_data)
    v_minus = np.min(y_data)
    v_pp = v_plus - v_minus

    # Define a distinct marker color (e.g., a nice muted orange so it stands out)
    m_color = C_PULSE

    # Draw the dashed horizontal lines
    ax.axhline(v_plus, color=m_color, linestyle='--', lw=1.2, alpha=0.7, zorder=2)
    ax.axhline(v_minus, color=m_color, linestyle='--', lw=1.2, alpha=0.7, zorder=2)

    # Add V+ and V- labels just outside the right edge of the plot
    # Using get_yaxis_transform() maps X to 0-1 (axes coords) and Y to data coords
    ax.text(1.01, v_plus, 'V+', transform=ax.get_yaxis_transform(),
            va='center', ha='left', color=m_color, fontsize=9, weight='bold', clip_on=False)
    ax.text(1.01, v_minus, 'V-', transform=ax.get_yaxis_transform(),
            va='center', ha='left', color=m_color, fontsize=9, weight='bold', clip_on=False)

    # Add the overall Vpp text box inside the plot, in the top right corner
    ax.text(0.98, 0.92, f'Vpp = {v_pp:.2f} {unit}', 
            transform=ax.transAxes, ha='right', va='top', 
            color=m_color, fontsize=9, weight='bold',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor=m_color, alpha=0.9))


# ── 4a. Subplot 1 – TX Driving signal ────────────────────────────────────────
ax = axs[0]

# Scale down by 1000 like we fixed earlier
y0_scaled = l0['y'] / 1000.0 

ax.plot(l0['x'], y0_scaled, color=C_DARKGREY, lw=1.1, zorder=3)

ax.set_title('Oscilloscope measurement of TX Driving signal', **TITLE_KW)
ax.set_xlim(-0.2, 1.0)
ax.set_ylim(-14, 14)
ax.set_xticks([-0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax.set_yticks([-10, -5, 0, 5, 10])
ax.set_ylabel('Amplitude (V)', **YLABEL_KW)
ax.grid(True, **GRID_KW)
ax.set_axisbelow(True)
ax.axhline(0, color=C_ZERO, lw=0.4, zorder=1)

# Add our Vpp markers
# add_vpp_markers(ax, y0_scaled, unit='V')


# ── 4b. Subplot 2 – RX 2-Layer PCB ───────────────────────────────────────────
ax = axs[1]
ax.plot(l1['x'], l1['y'], color=C_RAW, lw=1.1, zorder=3)

ax.set_title('Hydrophone reading of pulse transmitted using 2-Layer PCB', **TITLE_KW)
ax.set_xlim(37, 39.5)
ax.set_ylim(-2, 2)
ax.set_xticks([37, 37.5, 38, 38.5, 39, 39.5])
ax.set_yticks([-2, -1, 0, 1, 2])
ax.set_ylabel('Amplitude (mV)', **YLABEL_KW)
ax.grid(True, **GRID_KW)
ax.set_axisbelow(True)
ax.axhline(0, color=C_ZERO, lw=0.4, zorder=1)

# Add our Vpp markers
add_vpp_markers(ax, l1['y'], unit='mV')


# ── 4c. Subplot 3 – RX 4-Layer PCB ───────────────────────────────────────────
ax = axs[2]
ax.plot(l2['x'], l2['y'], color=C_MF, lw=1.1, zorder=3)

ax.set_title('Hydrophone reading of pulse transmitted using 4-Layer PCB', **TITLE_KW)
ax.set_xlim(37, 39.5)
ax.set_ylim(-2, 2)
ax.set_xticks([37, 37.5, 38, 38.5, 39, 39.5])
ax.set_yticks([-2, -1, 0, 1, 2])
ax.set_ylabel('Amplitude (mV)', **YLABEL_KW)
ax.set_xlabel(r'Time ($\mu$s)', labelpad=3)
ax.grid(True, **GRID_KW)
ax.set_axisbelow(True)
ax.axhline(0, color=C_ZERO, lw=0.4, zorder=1)

# Add our Vpp markers
add_vpp_markers(ax, l2['y'], unit='mV')


# ── 5. Save ───────────────────────────────────────────────────────────────────
OUT_DIR = 'out'
os.makedirs(OUT_DIR, exist_ok=True)

fig.savefig(os.path.join(OUT_DIR, 'waveform_v2.pdf'), bbox_inches='tight')
fig.savefig(os.path.join(OUT_DIR, 'waveform_v2.png'), dpi=300, bbox_inches='tight')
# plt.show()
print('Saved to', OUT_DIR)



