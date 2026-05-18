#!/usr/bin/env python3
"""
figure_directivity_portrait.py  –  Two A4-portrait PNG files.

All tunable numbers live in the USER-ADJUSTABLE PARAMETERS block.

Key polar transAxes geometry (theta_zero='N', thetamin=-90..90):
  Arc tip   transAxes y = 0.75   (r=1, theta=0 deg)
  Flat edge transAxes y = 0.25   (r=1, theta=+/-90 deg)
  Upper empty bounding-box space : y = 0.75 .. 1.00
  Lower empty bounding-box space : y = 0.00 .. 0.25

Title and subtitle are placed with ax.text() inside the upper empty space
(between y=0.75 and y=1.0) to keep them tight against the arc.
"""
import textwrap
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# ═══════════════════════════════════════════════════════════════════════════
#  USER-ADJUSTABLE PARAMETERS
# ═══════════════════════════════════════════════════════════════════════════

# ── Figure (A4 portrait) ────────────────────────────────────────────────────
FIG_W, FIG_H   = 8.27, 11.69   # inches

# ── Polar axes bounding box (figure-fraction units) ──────────────────────────
POLAR_LEFT     = 0.09
POLAR_W        = 0.82
POLAR_H_P1     = 0.52   # box height – PAGE 1 polar
POLAR_H_P2     = 0.44   # box height – PAGE 2 polars (same for both)

# ── Hanning weight plot (page 1 only) ────────────────────────────────────────
HAN_LEFT       = 0.11
HAN_W          = 0.78
HAN_H          = 0.20   # initial was 0.13

# ── Vertical gaps (figure-height fractions) ──────────────────────────────────
MARGIN_BOT     = 0.00
GAP_P1         = 0.04   # polar(a) bottom → Hanning top
GAP_P2         = 0.00   # polar(d) top → polar(c) bottom (page 2)
                         # title+subtitle live inside the polar box now,
                         # so this only needs to provide visual separation.

# ── Title typography ──────────────────────────────────────────────────────────
# Title is placed via ax.text() at (0.5, TITLE_Y_Px) in transAxes coords,
# inside the upper empty bounding-box space (arc tip = y=0.75, box top = y=1.0).
# Increase TITLE_Y_Px to push the title toward the bounding-box top.
# Use 'center' or 'left' for TITLE_HA.
TITLE_FS       = 11
TITLE_HA       = 'center'
TITLE_Y_P1     = 0.88   # transAxes y for title bottom – PAGE 1 polar
TITLE_Y_P2     = 0.93   # transAxes y for title bottom – PAGE 2 polars
                         # (slightly higher because box is smaller → subtitle is
                         #  proportionally taller; needs more clearance)
TITLE_PAD_B    = 6       # pt – Hanning (b) axes uses ax.set_title normally

# ── Subtitle (grey italic description) ────────────────────────────────────────
# Placed via ax.text() at (SUB_X, SUB_Y) in transAxes coords.
# Arc tip = transAxes y=0.75. SUB_Y slightly above 0.75 puts the subtitle
# just above the arc with minimal gap.
# Increase SUB_Y to add more space between arc and subtitle.
SUB_FS         = 8.5
SUB_HA         = 'center'
SUB_X          = 0.5
SUB_Y          = 0.78   # transAxes y for subtitle bottom  (arc tip = 0.75)
SUB_WRAP       = 88
SUB_LINESPACING = 1.4

# ── Legend inside polar plots ─────────────────────────────────────────────────
# Flat edge = transAxes y=0.25.  LEGEND_Y just below 0.25 places the legend
# top in the lower empty space, directly beneath the arc – close, no gap.
LEGEND_Y       = 0.20
LEGEND_FS      = 9.5

# ── Tick font sizes ───────────────────────────────────────────────────────────
R_FS           = 8.5
TH_FS          = 10.5

# ═══════════════════════════════════════════════════════════════════════════
#  PHYSICS
# ═══════════════════════════════════════════════════════════════════════════
c      = 1540.0;  fc = 10.5e6;  lam = c / fc
N      = 64;      p  = 75e-6
D_half = (N - 1) / 2 * p
x_el   = (np.arange(N) - (N - 1) / 2) * p
FLOOR  = -60.0;  Z_COMP = 5;  F_VALS = [2, 4, 8]

def array_factor(w, th_deg):
    k  = 2 * np.pi / lam
    ph = np.exp(1j * k * np.outer(x_el, np.sin(np.deg2rad(th_deg))))
    v  = np.abs(w @ ph);  return v / v.max()

def to_dB(v): return np.maximum(20 * np.log10(np.abs(v) + 1e-30), FLOOR)
def to_r(db): return np.clip((db - FLOOR) / (-FLOOR), 0, 1)

def hanning_fnum(z_mm, F):
    Rz = z_mm * 1e-3 / (2.0 * F);  w = np.zeros(N)
    ins = np.abs(x_el) <= Rz
    if ins.any():
        Re = np.max(np.abs(x_el[ins]))
        w[ins] = 0.5 * (1 + np.cos(np.pi * x_el[ins] / Re))
    return w

W_UNI = np.ones(N)
W_HAN = 0.5 * (1 + np.cos(np.pi * x_el / D_half))
W_F   = {F: hanning_fnum(Z_COMP, F) for F in F_VALS}
N_ACT = {F: int((W_F[F] > 0).sum()) for F in F_VALS}

TH     = np.linspace(-90, 90, 60000);  TH_rad = np.deg2rad(TH)
r_uni  = to_r(to_dB(array_factor(W_UNI, TH)))
r_han  = to_r(to_dB(array_factor(W_HAN, TH)))
r_f    = {F: to_r(to_dB(array_factor(W_F[F], TH))) for F in F_VALS}

C_UNI = '#2B2D42';  C_HAN = '#1D9E75';  C_REF = '#555555'
F_COL = {2: '#3A86C8', 4: '#BA7517', 8: '#C0392B'}

matplotlib.rcParams.update({
    'font.family': 'serif', 'font.size': 10,
    'axes.labelsize': 10, 'xtick.labelsize': 9.5, 'ytick.labelsize': 9.5,
    'axes.linewidth': 0.7,
    'axes.spines.top': False, 'axes.spines.right': False,
    'figure.facecolor': 'white',
    'grid.alpha': 0.30, 'grid.linewidth': 0.40, 'axes.grid': True,
})

# ═══════════════════════════════════════════════════════════════════════════
#  DRAWING HELPERS
# ═══════════════════════════════════════════════════════════════════════════
def setup_polar(ax):
    ax.set_theta_zero_location('N');  ax.set_theta_direction(-1)
    ax.set_thetamin(-90);             ax.set_thetamax(90)
    ax.set_rlim(0, 1)
    ax.set_yticks([0.25, 0.5, 0.75, 1.0])
    ax.set_yticklabels(['-45', '-30', '-15', '0 dB'], fontsize=R_FS, color='#555')
    ax.set_thetagrids([-90, -60, -30, 30, 60, 90])
    ax.tick_params(axis='x', labelsize=TH_FS)
    ax.grid(True, alpha=0.28, linewidth=0.40)


def polar_legend(ax, **kw):
    """Legend top at LEGEND_Y (< 0.25 = below arc flat edge)."""
    ax.legend(
        bbox_to_anchor=(0.5, LEGEND_Y),
        loc='upper center',
        bbox_transform=ax.transAxes,
        framealpha=0.93, edgecolor='#ddd', handlelength=2.0,
        fontsize=LEGEND_FS, **kw,
    )


def polar_title(ax, text, title_y):
    """Bold title placed inside the upper empty bounding-box space."""
    tx = 0.5 if TITLE_HA == 'center' else 0.0
    ax.text(tx, title_y, text,
            transform=ax.transAxes,
            fontsize=TITLE_FS, fontweight='bold',
            ha=TITLE_HA, va='bottom')


def add_subtitle(ax, text):
    """Grey italic description placed just above the arc tip (transAxes y=SUB_Y)."""
    wrapped = '\n'.join(textwrap.wrap(text, width=SUB_WRAP))
    ax.text(SUB_X, SUB_Y, wrapped,
            transform=ax.transAxes,
            fontsize=SUB_FS, style='italic', color='#555555',
            va='bottom', ha=SUB_HA, linespacing=SUB_LINESPACING)


# ═══════════════════════════════════════════════════════════════════════════
#  PANEL DRAWERS
# ═══════════════════════════════════════════════════════════════════════════
def draw_a(ax, title_y=TITLE_Y_P1, subtitle=True):
    setup_polar(ax)
    ax.fill(TH_rad, r_uni, color=C_UNI, alpha=0.08)
    ax.plot(TH_rad, r_uni, color=C_UNI, lw=0.75,
            label='Uniform  (64 ch)  --  first sidelobe approx. -13 dB')
    polar_legend(ax)
    polar_title(ax, '(a)  Receive beam pattern  --  uniform aperture', title_y)
    if subtitle:
        add_subtitle(ax,
            f'64 channels, pitch = 75 um, f0 = 10.5 MHz, '
            f'wavelength = {lam*1e6:.0f} um. '
            'Uniform (rectangular) window produces a sinc-shaped beam; '
            'the first sidelobe at approximately -13 dB causes bright targets '
            'to bleed into neighbouring image pixels.')


def draw_b(ax):
    ch = np.arange(N) + 1
    ax.fill_between(ch, W_UNI, step='mid', alpha=0.07, color=C_UNI)
    ax.step(ch, W_UNI, where='mid', color=C_UNI, lw=1.2, ls='--',
            label='Uniform  (rectangular)')
    ax.fill_between(ch, W_HAN, step='mid', alpha=0.16, color=C_HAN)
    ax.step(ch, W_HAN, where='mid', color=C_HAN, lw=1.6,
            label='Hanning taper')
    ax.set_xlabel('Channel index  $i$', fontsize=10)
    ax.set_ylabel('Weight  $w_i$', fontsize=10)
    ax.set_xlim(0.5, N + 0.5);  ax.set_ylim(0, 1.25)
    ax.text(0.50, 0.56,
            r'$w_i = \frac{1}{2}\left[1+\cos\left(\frac{\pi x_i}{R}\right)\right]$',
            transform=ax.transAxes, ha='center', fontsize=11, color=C_HAN,
            bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='none', alpha=0.9))
    ax.legend(fontsize=10, loc='upper right')
    ax.set_title('(b)  Hanning aperture weight function',
                 fontsize=TITLE_FS, fontweight='bold',
                 pad=TITLE_PAD_B, loc=TITLE_HA)


def draw_c(ax, title_y=TITLE_Y_P1, subtitle=True):
    setup_polar(ax)
    ax.fill(TH_rad, r_uni, color='#CCCCCC', alpha=0.04)
    ax.plot(TH_rad, r_uni, color='#BBBBBB', lw=0.65,
            label='Uniform  (reference)')
    ax.fill(TH_rad, r_han, color=C_HAN, alpha=0.11)
    ax.plot(TH_rad, r_han, color=C_HAN, lw=0.75,
            label='Hanning  (approx. -31 dB sidelobes)')
    polar_legend(ax, ncol=2)
    polar_title(ax, '(c)  Receive beam pattern  --  after Hanning apodisation', title_y)
    if subtitle:
        add_subtitle(ax,
            'Applying the Hanning taper suppresses sidelobes from approximately '
            '-13 dB to approximately -31 dB, at the cost of a slightly wider '
            'main lobe. The uniform reference is shown in light grey for comparison.')


def draw_d(ax, title_y=TITLE_Y_P1, subtitle=True):
    setup_polar(ax)
    ax.plot(TH_rad, r_han, color=C_REF, lw=0.60,
            label='Hanning, full aperture  (64 ch)')
    for F in F_VALS:
        lw  = 0.90 if F == 4 else 0.65
        lbl = f'F# = {F}  ({N_ACT[F]} ch)' + ('  [used]' if F == 4 else '')
        ax.fill(TH_rad, r_f[F], color=F_COL[F], alpha=0.07)
        ax.plot(TH_rad, r_f[F], color=F_COL[F], lw=lw, label=lbl)
    polar_legend(ax, ncol=2)
    polar_title(ax, f'(d)  Effect of F#  (Hanning, z = {Z_COMP} mm)', title_y)
    if subtitle:
        add_subtitle(ax,
            f'Hanning weights applied within the F#-limited aperture at '
            f'z = {Z_COMP} mm. '
            'A smaller F# activates more channels, narrowing the receive beam. '
            f'Reconstruction uses F# = 4 [used]; '
            f'full aperture reached at z ~ {4*N*p*1e3:.0f} mm. '
            'Full Hanning reference (no F# restriction) shown in grey.')


# ═══════════════════════════════════════════════════════════════════════════
#  PAGE 1  –  (a) uniform polar  +  (b) Hanning weight function
# ═══════════════════════════════════════════════════════════════════════════
han_rect    = [HAN_LEFT,   MARGIN_BOT,                   HAN_W,   HAN_H]
polar1_rect = [POLAR_LEFT, MARGIN_BOT + HAN_H + GAP_P1,  POLAR_W, POLAR_H_P1]

fig1 = plt.figure(figsize=(FIG_W, FIG_H), facecolor='white')
ax1b = fig1.add_axes(han_rect)
ax1a = fig1.add_axes(polar1_rect, projection='polar')
draw_a(ax1a, title_y=TITLE_Y_P1, subtitle=True)
draw_b(ax1b)
fig1.savefig('./out/figure_directivity_page1.png',
             dpi=200, bbox_inches='tight', facecolor='white')
print('Saved page 1')
plt.close(fig1)


# ═══════════════════════════════════════════════════════════════════════════
#  PAGE 2  –  (c) Hanning polar  +  (d) F# polar
# ═══════════════════════════════════════════════════════════════════════════
polar_d_rect = [POLAR_LEFT, MARGIN_BOT,                         POLAR_W, POLAR_H_P2]
polar_c_rect = [POLAR_LEFT, MARGIN_BOT + POLAR_H_P2 + GAP_P2,  POLAR_W, POLAR_H_P2]

fig2 = plt.figure(figsize=(FIG_W, FIG_H), facecolor='white')
ax2d = fig2.add_axes(polar_d_rect, projection='polar')
ax2c = fig2.add_axes(polar_c_rect, projection='polar')
draw_d(ax2d, title_y=TITLE_Y_P2, subtitle=True)
draw_c(ax2c, title_y=TITLE_Y_P2, subtitle=True)
fig2.savefig('./out/figure_directivity_page2.png',
             dpi=200, bbox_inches='tight', facecolor='white')
print('Saved page 2')
plt.close(fig2)
