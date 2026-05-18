# """
# resolution_figure.py  (v2)
# --------------------------
# Shows the two closest distinguishable peaks along each profile slice,
# with their individual –6 dB bands highlighted and a valley marker
# indicating whether they are resolved.

# Key changes from v1
# -------------------
#   - Peak detection via scipy.signal.find_peaks
#   - Finds the two CLOSEST peaks above a minimum prominence
#   - Zooms each profile panel around those two peaks
#   - Each peak gets its own coloured –6 dB band
#   - Valley between peaks is annotated (resolved / unresolved)

# MATLAB export block (add to end of your reconstruction script)
# --------------------------------------------------------------
#   target_z_mm =  30.0;   % depth of the horizontal slice (mm)
#   target_x_mm =   2.0;   % lateral position of the vertical slice (mm)

#   [~, iz] = min(abs(z_axis*1000 - target_z_mm));
#   [~, ix] = min(abs(x_axis*1000 - target_x_mm));

#   export = struct();
#   export.B_Mode_dB       = B_Mode_dB;
#   export.x_axis_mm       = x_axis * 1000;
#   export.z_axis_mm       = z_axis * 1000;
#   export.lateral_profile = B_Mode_dB(iz, :);
#   export.axial_profile   = B_Mode_dB(:, ix);
#   export.target_iz       = iz;
#   export.target_ix       = ix;
#   export.target_z_mm     = target_z_mm;
#   export.target_x_mm     = target_x_mm;
#   save('resolution_export.mat', '-struct', 'export');

# Usage
# -----
#   python resolution_figure.py                      # uses resolution_export.mat
#   python resolution_figure.py path/to/export.mat

# Figure size
# -----------
#   Edit FIGSIZE below.  A few useful presets:
#     Half A4 page  : (13.0, 6.2)   <- current default
#     Two-thirds A4 : (13.0, 7.8)
#     Full column   : ( 8.5, 4.2)
#   Width is set first, height second.

# Requirements:  pip install scipy numpy matplotlib
# """

# import sys, os
# import numpy as np
# import scipy.io as sio
# import scipy.signal as sig
# import matplotlib as mpl
# import matplotlib.pyplot as plt
# import matplotlib.gridspec as gridspec
# from   matplotlib.ticker import MultipleLocator, AutoMinorLocator

# # ─────────────────────────────────────────────────────────────────────────────
# # 0.  User configuration  <- all the knobs in one place
# # ─────────────────────────────────────────────────────────────────────────────
# MAT_FILE  = sys.argv[1] if len(sys.argv) > 1 else 'resolution_export.mat'
# OUT_DIR   = 'out'

# # Figure size (width, height) in inches
# FIGSIZE   = (13.0, 6.2)    # <- EDIT THIS to resize the whole figure

# # B-mode display
# DB_RANGE  = (-30, 0)       # dynamic range shown on B-mode
# CMAP      = 'inferno'         # 'gray' | 'hot' | 'inferno'

# # Resolution analysis
# DB_THRESH    = -6          # dB threshold for FWHM  (-6 dB = half-power)

# # Peak detection — tune these if peaks are missed or noise is picked up:
# #   PEAK_MIN_PROMINENCE : how many dB a peak must rise above its surroundings
# #   PEAK_MIN_DISTANCE   : minimum separation between peaks in samples
# PEAK_MIN_PROMINENCE = 6
# PEAK_MIN_DISTANCE   = 1

# # Padding around the two-peak zoom window (mm)
# ZOOM_PAD_MM = 3.0          # increase if annotation labels are clipped

# # Phantom caption on the B-mode panel (set to None to hide)
# PHANTOM_LABEL = 'CIRS 040GSE — resolution group'

# # ─────────────────────────────────────────────────────────────────────────────
# # 1.  Style
# # ─────────────────────────────────────────────────────────────────────────────
# mpl.rcParams.update({
#     'font.family':      'serif',
#     'font.serif':       ['DejaVu Serif', 'Times New Roman', 'Times'],
#     'font.size':         9.5,
#     'axes.labelsize':    9.5,
#     'axes.titlesize':    10,
#     'xtick.labelsize':   8.5,
#     'ytick.labelsize':   8.5,
#     'axes.linewidth':    0.7,
#     'xtick.major.width': 0.7,
#     'ytick.major.width': 0.7,
#     'xtick.major.size':  3.5,
#     'ytick.major.size':  3.5,
#     'xtick.direction':  'in',
#     'ytick.direction':  'in',
#     'mathtext.fontset': 'cm',
#     'pdf.fonttype':      42,
#     'ps.fonttype':       42,
# })

# C_CROSS      = '#E05A2B'
# C_VALLEY     = '#444444'
# C_PROF       = '#1A1A1A'
# PEAK_COLORS  = ['#1158A0', '#C8580A']   # blue, orange — one per peak
# FILL_COLORS  = ['#C5DEF0', '#FAE0C5']

# # ─────────────────────────────────────────────────────────────────────────────
# # 2.  Load data
# # ─────────────────────────────────────────────────────────────────────────────
# print(f'Loading: {MAT_FILE}')
# d = sio.loadmat(MAT_FILE, squeeze_me=True, struct_as_record=False)

# B_dB     = np.asarray(d['B_Mode_dB'],       dtype=float)
# x_mm     = np.asarray(d['x_axis_mm'],       dtype=float).ravel()
# z_mm     = np.asarray(d['z_axis_mm'],       dtype=float).ravel()
# lat_prof = np.asarray(d['lateral_profile'], dtype=float).ravel()
# ax_prof  = np.asarray(d['axial_profile'],   dtype=float).ravel()
# iz       = int(d['target_iz']) - 1      # MATLAB 1-based -> Python 0-based
# ix       = int(d['target_ix']) - 1
# tz_mm    = float(d['target_z_mm'])
# tx_mm    = float(d['target_x_mm'])

# print(f'  B-mode: {B_dB.shape}   target slice: x={tx_mm:.1f} mm, z={tz_mm:.1f} mm')

# # ─────────────────────────────────────────────────────────────────────────────
# # 3.  Peak analysis helpers
# # ─────────────────────────────────────────────────────────────────────────────

# def find_two_closest_peaks(axis, profile):
#     """
#     Find all prominent peaks, then return the indices of the two closest
#     together — the hardest pair to resolve, i.e. the resolution limit.
#     Returns (None, None) if fewer than 2 peaks are found.
#     """
#     peaks, _ = sig.find_peaks(
#         profile,
#         prominence=PEAK_MIN_PROMINENCE,
#         distance=PEAK_MIN_DISTANCE,
#     )
#     if len(peaks) < 2:
#         print(f'  WARNING: only {len(peaks)} peak(s) found. '
#               f'Try lowering PEAK_MIN_PROMINENCE (currently {PEAK_MIN_PROMINENCE}).')
#         return None, None

#     # Find closest pair by axis distance
#     best_gap, best_pair = np.inf, (peaks[0], peaks[1])
#     for i in range(len(peaks)):
#         for j in range(i + 1, len(peaks)):
#             gap = abs(axis[peaks[i]] - axis[peaks[j]])
#             if gap < best_gap:
#                 best_gap  = gap
#                 best_pair = (peaks[i], peaks[j])

#     return best_pair


# def fwhm_of_peak(axis, profile, peak_idx):
#     """
#     Compute the –6 dB full-width of a single peak by linear interpolation
#     at the threshold crossings.
#     Returns (fwhm_mm, left_mm, right_mm, level_dB).
#     """
#     level = profile[peak_idx] + DB_THRESH

#     # Left crossing
#     left_seg   = profile[:peak_idx + 1]
#     below_left = np.where(left_seg <= level)[0]
#     if len(below_left) == 0:
#         left = axis[0]
#     else:
#         j = below_left[-1]
#         if j + 1 <= peak_idx:
#             dy   = profile[j + 1] - profile[j]
#             t    = (level - profile[j]) / dy if dy != 0 else 0.5
#             left = axis[j] + t * (axis[j + 1] - axis[j])
#         else:
#             left = axis[j]

#     # Right crossing
#     right_seg   = profile[peak_idx:]
#     below_right = np.where(right_seg <= level)[0]
#     if len(below_right) == 0:
#         right = axis[-1]
#     else:
#         j = below_right[0] + peak_idx
#         if j - 1 >= 0:
#             dy    = profile[j] - profile[j - 1]
#             t     = (level - profile[j - 1]) / dy if dy != 0 else 0.5
#             right = axis[j - 1] + t * (axis[j] - axis[j - 1])
#         else:
#             right = axis[j]

#     return right - left, left, right, level


# def valley_between(axis, profile, peak_a, peak_b):
#     """Return (position_mm, amplitude_dB) of the minimum between two peaks."""
#     lo, hi = sorted([peak_a, peak_b])
#     vi     = lo + int(np.argmin(profile[lo:hi + 1]))
#     return axis[vi], profile[vi]


# # ─────────────────────────────────────────────────────────────────────────────
# # 4.  Run analysis
# # ─────────────────────────────────────────────────────────────────────────────
# print('\n-- Lateral profile --')
# lat_p1, lat_p2 = find_two_closest_peaks(x_mm, lat_prof)

# print('\n-- Axial profile --')
# ax_p1, ax_p2 = find_two_closest_peaks(z_mm, ax_prof)

# # ─────────────────────────────────────────────────────────────────────────────
# # 5.  Build figure
# # ─────────────────────────────────────────────────────────────────────────────
# fig = plt.figure(figsize=FIGSIZE)
# gs  = gridspec.GridSpec(
#     2, 2,
#     width_ratios=[1.15, 1],
#     height_ratios=[1, 1],
#     hspace=0.52, wspace=0.38,
#     left=0.07, right=0.96, top=0.93, bottom=0.10,
# )
# ax_bm  = fig.add_subplot(gs[:, 0])
# ax_lat = fig.add_subplot(gs[0, 1])
# ax_ax  = fig.add_subplot(gs[1, 1])


# # ── 5a. B-mode ────────────────────────────────────────────────────────────────
# im = ax_bm.imshow(
#     B_dB,
#     aspect='auto',
#     extent=[x_mm[0], x_mm[-1], z_mm[-1], z_mm[0]],
#     vmin=DB_RANGE[0], vmax=DB_RANGE[1],
#     cmap=CMAP, interpolation='bilinear',
# )
# cbar = fig.colorbar(im, ax=ax_bm, fraction=0.046, pad=0.03, aspect=38)
# cbar.set_label('Normalised amplitude (dB)', fontsize=8.5)
# cbar.ax.tick_params(labelsize=8)
# cbar.set_ticks(range(DB_RANGE[0], DB_RANGE[1] + 1, 10))

# # Slice crosshairs (colour-coded to match the two panels)
# ax_bm.axhline(tz_mm, color=C_CROSS,          lw=0.9, ls='--', alpha=0.85,
#               label=f'Lateral slice  z = {tz_mm:.0f} mm')
# ax_bm.axvline(tx_mm, color=PEAK_COLORS[1],   lw=0.9, ls=':',  alpha=0.85,
#               label=f'Axial slice  x = {tx_mm:.1f} mm')
# ax_bm.plot(tx_mm, tz_mm, 'o',
#            color='none', mec=C_CROSS, mew=1.3, ms=10, zorder=5)

# ax_bm.set_xlabel('Lateral (mm)', labelpad=3)
# ax_bm.set_ylabel('Axial / Depth (mm)', labelpad=3)
# ax_bm.set_title('B-mode reconstruction — DAS + CF', pad=5)
# ax_bm.legend(fontsize=7.5, framealpha=0.45, loc='lower right',
#              facecolor='black', labelcolor='white', edgecolor='none')

# if PHANTOM_LABEL:
#     ax_bm.text(0.02, 0.98, PHANTOM_LABEL, transform=ax_bm.transAxes,
#                fontsize=7.5, color='white', va='top', ha='left',
#                bbox=dict(boxstyle='round,pad=0.25', fc='#1A1A1A',
#                          alpha=0.55, ec='none'))


# # ── 5b/5c. Profile panels ─────────────────────────────────────────────────────

# def draw_profile_panel(ax, axis, profile, p1_idx, p2_idx,
#                        xlabel, title, is_axial=False):
#     """
#     Draw profile with two –6 dB bands, valley annotation, FWHM arrows,
#     zoomed to the two-peak region.
#     """
#     if p1_idx is None:
#         ax.plot(axis, profile, color=C_PROF, lw=1.4)
#         ax.set_ylim(DB_RANGE)
#         ax.set_xlabel(xlabel, labelpad=2)
#         ax.set_title(title, pad=4)
#         return

#     # Ensure left-to-right ordering
#     if axis[p1_idx] > axis[p2_idx]:
#         p1_idx, p2_idx = p2_idx, p1_idx

#     fwhm1, l1, r1, lev1 = fwhm_of_peak(axis, profile, p1_idx)
#     fwhm2, l2, r2, lev2 = fwhm_of_peak(axis, profile, p2_idx)
#     val_pos, val_amp     = valley_between(axis, profile, p1_idx, p2_idx)

#     is_resolved = (val_amp < lev1) and (val_amp < lev2)
#     res_str     = 'Resolved' if is_resolved else 'Unresolved'
#     res_col     = '#1C7A3A' if is_resolved else '#AA2222'

#     # Zoom window
#     x_lo = min(l1, l2) - ZOOM_PAD_MM
#     x_hi = max(r1, r2) + ZOOM_PAD_MM

#     # ── Shaded -6 dB bands ────────────────────────────────────────────────────
#     for (l, r, lev, pc, fc) in [
#         (l1, r1, lev1, PEAK_COLORS[0], FILL_COLORS[0]),
#         (l2, r2, lev2, PEAK_COLORS[1], FILL_COLORS[1]),
#     ]:
#         band_mask = (axis >= l) & (axis <= r) & (profile >= lev)
#         ax.fill_between(axis, DB_RANGE[0], profile,
#                         where=band_mask,
#                         color=fc, alpha=0.80, zorder=2,
#                         label=f'−{abs(DB_THRESH)} dB band')
#         # Dashed vertical lines at the two FWHM edges
#         ax.vlines([l, r], DB_RANGE[0], lev,
#                   color=pc, lw=0.85, ls=':', alpha=0.65, zorder=3)

#     # ── Profile curve ─────────────────────────────────────────────────────────
#     ax.plot(axis, profile, color=C_PROF, lw=1.5, zorder=5)

#     # ── –6 dB horizontal lines ────────────────────────────────────────────────
#     for lev, pc in [(lev1, PEAK_COLORS[0]), (lev2, PEAK_COLORS[1])]:
#         ax.hlines(lev, x_lo, x_hi, color=pc, lw=0.85, ls='--',
#                   alpha=0.60, zorder=3)

#     # ── FWHM double-headed arrows (staggered vertically so they don't overlap) ─
#     dy_span = 2 - DB_RANGE[0]
#     for (l, r, fwhm, pc), frac in zip(
#         [(l1, r1, fwhm1, PEAK_COLORS[0]),
#          (l2, r2, fwhm2, PEAK_COLORS[1])],
#         [0.07, 0.15],
#     ):
#         y_arr = DB_RANGE[0] + dy_span * frac
#         ax.annotate('', xy=(r, y_arr), xytext=(l, y_arr),
#                     arrowprops=dict(arrowstyle='<->', color=pc, lw=1.3,
#                                    shrinkA=0, shrinkB=0),
#                     zorder=6)
#         ax.text((l + r) / 2, y_arr - dy_span * 0.025,
#                 f'{fwhm * 1000:.0f} µm',
#                 ha='center', va='top', fontsize=8.5, color=pc,
#                 weight='bold', zorder=7)

#     # ── Valley marker ─────────────────────────────────────────────────────────
#     ax.plot(val_pos, val_amp, 'v', color=C_VALLEY, ms=6.5, zorder=6,
#             clip_on=False)
#     ax.text(val_pos, val_amp - dy_span * 0.03,
#             f'{val_amp:.1f} dB  —  {res_str}',
#             ha='center', va='top', fontsize=8, color=res_col,
#             weight='bold', zorder=7,
#             bbox=dict(boxstyle='round,pad=0.25', fc='white',
#                       alpha=0.80, ec=res_col, lw=0.6))

#     # ── Axis formatting ───────────────────────────────────────────────────────
#     ax.set_xlim(x_lo, x_hi)
#     ax.set_ylim(DB_RANGE[0], 2)
#     ax.set_xlabel(xlabel, labelpad=2)
#     ax.set_ylabel('Amplitude (dB)', labelpad=2)
#     ax.set_title(title, pad=4)
#     ax.xaxis.set_minor_locator(AutoMinorLocator(2))
#     ax.yaxis.set_minor_locator(MultipleLocator(5))
#     ax.tick_params(which='minor', length=1.8, width=0.5)

#     # dB labels outside the right spine
#     for lev, pc in [(lev1, PEAK_COLORS[0]), (lev2, PEAK_COLORS[1])]:
#         ax.text(x_hi + 0.01 * (x_hi - x_lo), lev,
#                 f'{DB_THRESH} dB',
#                 ha='left', va='center', fontsize=7.5, color=pc,
#                 clip_on=False)

#     # Console summary
#     tag = 'Axial' if is_axial else 'Lateral'
#     sep = abs(axis[p2_idx] - axis[p1_idx])
#     print(f'  {tag} FWHM peak-1 = {fwhm1 * 1000:.0f} µm')
#     print(f'  {tag} FWHM peak-2 = {fwhm2 * 1000:.0f} µm')
#     print(f'  {tag} separation  = {sep  * 1000:.0f} µm')
#     print(f'  {tag} valley      = {val_amp:.1f} dB  ->  {res_str}')


# print('\n-- Drawing lateral panel --')
# draw_profile_panel(
#     ax_lat, x_mm, lat_prof, lat_p1, lat_p2,
#     xlabel='Lateral (mm)',
#     title=f'Lateral profile  at  z = {tz_mm:.1f} mm',
# )

# print('\n-- Drawing axial panel --')
# draw_profile_panel(
#     ax_ax, z_mm, ax_prof, ax_p1, ax_p2,
#     xlabel='Depth (mm)',
#     title=f'Axial profile  at  x = {tx_mm:.1f} mm',
#     is_axial=True,
# )

# # Panel labels (a), (b), (c)
# for ax, lbl in [(ax_bm, '(a)'), (ax_lat, '(b)'), (ax_ax, '(c)')]:
#     ax.text(-0.09, 1.02, lbl, transform=ax.transAxes,
#             fontsize=11, weight='bold', color='#1A1A1A', va='bottom')

# # ─────────────────────────────────────────────────────────────────────────────
# # 6.  Save
# # ─────────────────────────────────────────────────────────────────────────────
# os.makedirs(OUT_DIR, exist_ok=True)
# stem = os.path.join(OUT_DIR, 'resolution_figure_v2')
# fig.savefig(stem + '.pdf', bbox_inches='tight')
# fig.savefig(stem + '.png', dpi=300, bbox_inches='tight')
# plt.show()
# print(f'\nSaved -> {stem}.pdf  /  {stem}.png')







#  --------------------------------------------------------------------------------
# """
# resolution_figure.py  (v2)
# --------------------------
# Shows the two closest distinguishable peaks along each profile slice,
# with their individual –6 dB bands highlighted and a valley marker
# indicating whether they are resolved.

# Key changes from v1
# -------------------
#   - Peak detection via scipy.signal.find_peaks
#   - Finds the two CLOSEST peaks above a minimum prominence
#   - Zooms each profile panel around those two peaks
#   - Each peak gets its own coloured –6 dB band
#   - Valley between peaks is annotated (resolved / unresolved)

# MATLAB export block (add to end of your reconstruction script)
# --------------------------------------------------------------
#   target_z_mm =  30.0;   % depth of the horizontal slice (mm)
#   target_x_mm =   2.0;   % lateral position of the vertical slice (mm)

#   [~, iz] = min(abs(z_axis*1000 - target_z_mm));
#   [~, ix] = min(abs(x_axis*1000 - target_x_mm));

#   export = struct();
#   export.B_Mode_dB       = B_Mode_dB;
#   export.x_axis_mm       = x_axis * 1000;
#   export.z_axis_mm       = z_axis * 1000;
#   export.lateral_profile = B_Mode_dB(iz, :);
#   export.axial_profile   = B_Mode_dB(:, ix);
#   export.target_iz       = iz;
#   export.target_ix       = ix;
#   export.target_z_mm     = target_z_mm;
#   export.target_x_mm     = target_x_mm;
#   save('resolution_export.mat', '-struct', 'export');

# Usage
# -----
#   python resolution_figure.py                      # uses resolution_export.mat
#   python resolution_figure.py path/to/export.mat

# Figure size
# -----------
#   Edit FIGSIZE below.  A few useful presets:
#     Half A4 page  : (13.0, 6.8)   <- current default
#     Two-thirds A4 : (13.0, 8.5)
#     Full column   : ( 8.5, 5.0)
#   Width is set first, height second.

# Requirements:  pip install scipy numpy matplotlib
# """

# import sys, os
# import numpy as np
# import scipy.io as sio
# import scipy.signal as sig
# import matplotlib as mpl
# import matplotlib.pyplot as plt
# import matplotlib.gridspec as gridspec
# from   matplotlib.ticker import MultipleLocator, AutoMinorLocator

# # ─────────────────────────────────────────────────────────────────────────────
# # 0.  User configuration  <- all the knobs in one place
# # ─────────────────────────────────────────────────────────────────────────────
# MAT_FILE  = sys.argv[1] if len(sys.argv) > 1 else 'resolution_export.mat'
# OUT_DIR   = 'out'

# # Figure size (width, height) in inches
# # FIGSIZE   = (13.0, 6.8)    # <- EDIT THIS to resize the whole figure
# FIGSIZE   = (17.0, 6.8)

# # B-mode display
# DB_RANGE  = (-50, 0)       # dynamic range shown on B-mode
# CMAP      = 'inferno'         # 'gray' | 'hot' | 'inferno'

# # Resolution analysis
# DB_THRESH    = -6          # dB threshold for FWHM  (-6 dB = half-power)

# # ── Peak detection ────────────────────────────────────────────────────────────
# # THIS IS THE MOST IMPORTANT PARAMETER TO SET:
# #   PEAK_MIN_HEIGHT : absolute dB floor — any peak below this is ignored.
# #                     Look at your B-mode, read off the dB value of the real
# #                     wire targets from the colourbar, then set this ~10 dB
# #                     below that so -90 dB noise is excluded.
# #                     Example: targets at ~-5 dB   ->  set to -20
# #                              targets at ~-65 dB  ->  set to -75
# PEAK_MIN_HEIGHT     = -40  # <- SET THIS FIRST based on your B-mode colourbar

# # Secondary controls (usually fine as-is):
# #   PEAK_MIN_PROMINENCE : peak must also rise this many dB above surroundings
# #   PEAK_MIN_DISTANCE   : minimum separation in samples (avoids double-picks)
# PEAK_MIN_PROMINENCE = 3
# PEAK_MIN_DISTANCE   = 1

# # Padding around the two-peak zoom window (mm)
# ZOOM_PAD_MM = 6.0          # increase if annotation labels are clipped

# # Phantom caption on the B-mode panel (set to None to hide)
# PHANTOM_LABEL = 'CIRS 040GSE — resolution group'

# # ─────────────────────────────────────────────────────────────────────────────
# # 1.  Style
# # ─────────────────────────────────────────────────────────────────────────────
# mpl.rcParams.update({
#     'font.family':      'serif',
#     'font.serif':       ['DejaVu Serif', 'Times New Roman', 'Times'],
#     'font.size':         9.5,
#     'axes.labelsize':    9.5,
#     'axes.titlesize':    10,
#     'xtick.labelsize':   8.5,
#     'ytick.labelsize':   8.5,
#     'axes.linewidth':    0.7,
#     'xtick.major.width': 0.7,
#     'ytick.major.width': 0.7,
#     'xtick.major.size':  3.5,
#     'ytick.major.size':  3.5,
#     'xtick.direction':  'in',
#     'ytick.direction':  'in',
#     'mathtext.fontset': 'cm',
#     'pdf.fonttype':      42,
#     'ps.fonttype':       42,
# })

# C_CROSS      = '#E05A2B'
# C_VALLEY     = '#444444'
# C_PROF       = '#1A1A1A'
# PEAK_COLORS  = ['#1158A0', '#C8580A']   # blue, orange — one per peak
# FILL_COLORS  = ['#C5DEF0', '#FAE0C5']

# # ─────────────────────────────────────────────────────────────────────────────
# # 2.  Load data
# # ─────────────────────────────────────────────────────────────────────────────
# print(f'Loading: {MAT_FILE}')
# d = sio.loadmat(MAT_FILE, squeeze_me=True, struct_as_record=False)

# B_dB     = np.asarray(d['B_Mode_dB'],       dtype=float)
# x_mm     = np.asarray(d['x_axis_mm'],       dtype=float).ravel()
# z_mm     = np.asarray(d['z_axis_mm'],       dtype=float).ravel()
# lat_prof = np.asarray(d['lateral_profile'], dtype=float).ravel()
# ax_prof  = np.asarray(d['axial_profile'],   dtype=float).ravel()
# iz       = int(d['target_iz']) - 1      # MATLAB 1-based -> Python 0-based
# ix       = int(d['target_ix']) - 1
# # tz_mm    = float(d['target_z_mm'])
# # tx_mm    = float(d['target_x_mm'])
# tz_mm = 29.75
# tx_mm = 1.9

# print(f'  B-mode: {B_dB.shape}   target slice: x={tx_mm:.1f} mm, z={tz_mm:.1f} mm')

# # ─────────────────────────────────────────────────────────────────────────────
# # 3.  Peak analysis helpers
# # ─────────────────────────────────────────────────────────────────────────────

# def find_two_closest_peaks(axis, profile):
#     """
#     Find all prominent peaks, then return the indices of the two closest
#     together — the hardest pair to resolve, i.e. the resolution limit.
#     Returns (None, None) if fewer than 2 peaks are found.
#     """
#     peaks, _ = sig.find_peaks(
#         profile,
#         height=PEAK_MIN_HEIGHT,         # absolute floor — excludes -90 dB noise
#         prominence=PEAK_MIN_PROMINENCE,
#         distance=PEAK_MIN_DISTANCE,
#     )
#     if len(peaks) < 2:
#         found_str = ', '.join(f'{profile[p]:.1f}' for p in peaks) or 'none'
#         print(f'  WARNING: only {len(peaks)} peak(s) above PEAK_MIN_HEIGHT='
#               f'{PEAK_MIN_HEIGHT} dB (found: {found_str} dB). '
#               f'Lower PEAK_MIN_HEIGHT until real wire targets appear.')
#         return None, None

#     # Find closest pair by axis distance
#     best_gap, best_pair = np.inf, (peaks[0], peaks[1])
#     for i in range(len(peaks)):
#         for j in range(i + 1, len(peaks)):
#             gap = abs(axis[peaks[i]] - axis[peaks[j]])
#             if gap < best_gap:
#                 best_gap  = gap
#                 best_pair = (peaks[i], peaks[j])

#     return best_pair


# def fwhm_of_peak(axis, profile, peak_idx):
#     """
#     Compute the –6 dB full-width of a single peak by linear interpolation
#     at the threshold crossings.
#     Returns (fwhm_mm, left_mm, right_mm, level_dB).
#     """
#     level = profile[peak_idx] + DB_THRESH

#     # Left crossing
#     left_seg   = profile[:peak_idx + 1]
#     below_left = np.where(left_seg <= level)[0]
#     if len(below_left) == 0:
#         left = axis[0]
#     else:
#         j = below_left[-1]
#         if j + 1 <= peak_idx:
#             dy   = profile[j + 1] - profile[j]
#             t    = (level - profile[j]) / dy if dy != 0 else 0.5
#             left = axis[j] + t * (axis[j + 1] - axis[j])
#         else:
#             left = axis[j]

#     # Right crossing
#     right_seg   = profile[peak_idx:]
#     below_right = np.where(right_seg <= level)[0]
#     if len(below_right) == 0:
#         right = axis[-1]
#     else:
#         j = below_right[0] + peak_idx
#         if j - 1 >= 0:
#             dy    = profile[j] - profile[j - 1]
#             t     = (level - profile[j - 1]) / dy if dy != 0 else 0.5
#             right = axis[j - 1] + t * (axis[j] - axis[j - 1])
#         else:
#             right = axis[j]

#     return right - left, left, right, level


# def valley_between(axis, profile, peak_a, peak_b):
#     """Return (position_mm, amplitude_dB) of the minimum between two peaks."""
#     lo, hi = sorted([peak_a, peak_b])
#     vi     = lo + int(np.argmin(profile[lo:hi + 1]))
#     return axis[vi], profile[vi]


# # ─────────────────────────────────────────────────────────────────────────────
# # 4.  Run analysis
# # ─────────────────────────────────────────────────────────────────────────────
# print('\n-- Lateral profile --')
# lat_p1, lat_p2 = find_two_closest_peaks(x_mm, lat_prof)

# print('\n-- Axial profile --')
# ax_p1, ax_p2 = find_two_closest_peaks(z_mm, ax_prof)

# # ─────────────────────────────────────────────────────────────────────────────
# # 5.  Build figure
# # ─────────────────────────────────────────────────────────────────────────────
# fig = plt.figure(figsize=FIGSIZE)
# gs  = gridspec.GridSpec(
#     2, 2,
#     width_ratios=[1.0, 1.35],   # right column (profiles) slightly wider than left; 1.15 original
#     height_ratios=[1, 1],
#     hspace=0.55, wspace=0.42,
#     left=0.07, right=0.97, top=0.93, bottom=0.10,
# )
# ax_bm  = fig.add_subplot(gs[:, 0])
# ax_lat = fig.add_subplot(gs[0, 1])
# ax_ax  = fig.add_subplot(gs[1, 1])


# # ── 5a. B-mode ────────────────────────────────────────────────────────────────
# im = ax_bm.imshow(
#     B_dB,
#     aspect='auto',
#     extent=[x_mm[0], x_mm[-1], z_mm[-1], z_mm[0]],
#     vmin=DB_RANGE[0], vmax=DB_RANGE[1],
#     cmap=CMAP, interpolation='bilinear',
# )
# cbar = fig.colorbar(im, ax=ax_bm, fraction=0.046, pad=0.03, aspect=38)
# cbar.set_label('Normalised amplitude (dB)', fontsize=8.5)
# cbar.ax.tick_params(labelsize=8)
# cbar.set_ticks(range(DB_RANGE[0], DB_RANGE[1] + 1, 10))

# # Slice crosshairs (colour-coded to match the two panels)
# ax_bm.axhline(tz_mm, color=C_CROSS,          lw=0.9, ls='--', alpha=0.85,
#               label=f'Lateral slice  z = {tz_mm:.0f} mm')
# ax_bm.axvline(tx_mm, color=PEAK_COLORS[1],   lw=0.9, ls=':',  alpha=0.85,
#               label=f'Axial slice  x = {tx_mm:.1f} mm')
# ax_bm.plot(tx_mm, tz_mm, 'o',
#            color='none', mec=C_CROSS, mew=1.3, ms=10, zorder=5)

# ax_bm.set_xlabel('Lateral (mm)', labelpad=3)
# ax_bm.set_ylabel('Axial / Depth (mm)', labelpad=3)
# ax_bm.set_title('B-mode reconstruction — DAS + CF', pad=5)
# ax_bm.legend(fontsize=7.5, framealpha=0.45, loc='lower right',
#              facecolor='black', labelcolor='white', edgecolor='none')

# if PHANTOM_LABEL:
#     ax_bm.text(0.02, 0.98, PHANTOM_LABEL, transform=ax_bm.transAxes,
#                fontsize=7.5, color='white', va='top', ha='left',
#                bbox=dict(boxstyle='round,pad=0.25', fc='#1A1A1A',
#                          alpha=0.55, ec='none'))


# # ── 5b/5c. Profile panels ─────────────────────────────────────────────────────

# def draw_profile_panel(ax, axis, profile, p1_idx, p2_idx,
#                        xlabel, title, is_axial=False):
#     """
#     Draw profile with two –6 dB bands, valley annotation, FWHM arrows,
#     zoomed to the two-peak region.
#     """
#     if p1_idx is None:
#         ax.plot(axis, profile, color=C_PROF, lw=1.4)
#         ax.set_ylim(DB_RANGE)
#         ax.set_xlabel(xlabel, labelpad=2)
#         ax.set_title(title, pad=4)
#         return

#     # Ensure left-to-right ordering
#     if axis[p1_idx] > axis[p2_idx]:
#         p1_idx, p2_idx = p2_idx, p1_idx

#     fwhm1, l1, r1, lev1 = fwhm_of_peak(axis, profile, p1_idx)
#     fwhm2, l2, r2, lev2 = fwhm_of_peak(axis, profile, p2_idx)
#     val_pos, val_amp     = valley_between(axis, profile, p1_idx, p2_idx)

#     is_resolved = (val_amp < lev1) and (val_amp < lev2)
#     res_str     = 'Resolved' if is_resolved else 'Unresolved'
#     res_col     = '#1C7A3A' if is_resolved else '#AA2222'

#     # Zoom window
#     x_lo = min(l1, l2) - ZOOM_PAD_MM
#     x_hi = max(r1, r2) + ZOOM_PAD_MM

#     # ── Shaded -6 dB bands ────────────────────────────────────────────────────
#     for (l, r, lev, pc, fc) in [
#         (l1, r1, lev1, PEAK_COLORS[0], FILL_COLORS[0]),
#         (l2, r2, lev2, PEAK_COLORS[1], FILL_COLORS[1]),
#     ]:
#         band_mask = (axis >= l) & (axis <= r) & (profile >= lev)
#         ax.fill_between(axis, DB_RANGE[0], profile,
#                         where=band_mask,
#                         color=fc, alpha=0.80, zorder=2,
#                         label=f'−{abs(DB_THRESH)} dB band')
#         # Dashed vertical lines at the two FWHM edges
#         ax.vlines([l, r], DB_RANGE[0], lev,
#                   color=pc, lw=0.85, ls=':', alpha=0.65, zorder=3)

#     # ── Profile curve ─────────────────────────────────────────────────────────
#     ax.plot(axis, profile, color=C_PROF, lw=1.5, zorder=5)

#     # ── –6 dB horizontal lines ────────────────────────────────────────────────
#     for lev, pc in [(lev1, PEAK_COLORS[0]), (lev2, PEAK_COLORS[1])]:
#         ax.hlines(lev, x_lo, x_hi, color=pc, lw=0.85, ls='--',
#                   alpha=0.60, zorder=3)

#     # ── FWHM double-headed arrows (staggered vertically so they don't overlap) ─
#     dy_span = 2 - DB_RANGE[0]
#     for (l, r, fwhm, pc), frac in zip(
#         [(l1, r1, fwhm1, PEAK_COLORS[0]),
#          (l2, r2, fwhm2, PEAK_COLORS[1])],
#         [0.07, 0.15],
#     ):
#         y_arr = DB_RANGE[0] + dy_span * frac
#         ax.annotate('', xy=(r, y_arr), xytext=(l, y_arr),
#                     arrowprops=dict(arrowstyle='<->', color=pc, lw=1.3,
#                                    shrinkA=0, shrinkB=0),
#                     zorder=6)
#         ax.text((l + r) / 2, y_arr - dy_span * 0.025,
#                 f'{fwhm * 1000:.0f} µm',
#                 ha='center', va='top', fontsize=8.5, color=pc,
#                 weight='bold', zorder=7)

#     # ── Valley marker ─────────────────────────────────────────────────────────
#     ax.plot(val_pos, val_amp, 'v', color=C_VALLEY, ms=6.5, zorder=6,
#             clip_on=False)
#     ax.text(val_pos, val_amp - dy_span * 0.03,
#             f'{val_amp:.1f} dB  —  {res_str}',
#             ha='center', va='top', fontsize=8, color=res_col,
#             weight='bold', zorder=7,
#             bbox=dict(boxstyle='round,pad=0.25', fc='white',
#                       alpha=0.80, ec=res_col, lw=0.6))

#     # ── Axis formatting ───────────────────────────────────────────────────────
#     ax.set_xlim(x_lo, x_hi)
#     ax.set_ylim(DB_RANGE[0], 2)
#     ax.set_xlabel(xlabel, labelpad=2)
#     ax.set_ylabel('Amplitude (dB)', labelpad=2)
#     ax.set_title(title, pad=4)
#     ax.xaxis.set_minor_locator(AutoMinorLocator(2))
#     ax.yaxis.set_minor_locator(MultipleLocator(5))
#     ax.tick_params(which='minor', length=1.8, width=0.5)

#     # dB labels outside the right spine
#     for lev, pc in [(lev1, PEAK_COLORS[0]), (lev2, PEAK_COLORS[1])]:
#         ax.text(x_hi + 0.01 * (x_hi - x_lo), lev,
#                 f'{DB_THRESH} dB',
#                 ha='left', va='center', fontsize=7.5, color=pc,
#                 clip_on=False)

#     # Console summary
#     tag = 'Axial' if is_axial else 'Lateral'
#     sep = abs(axis[p2_idx] - axis[p1_idx])
#     print(f'  {tag} FWHM peak-1 = {fwhm1 * 1000:.0f} µm')
#     print(f'  {tag} FWHM peak-2 = {fwhm2 * 1000:.0f} µm')
#     print(f'  {tag} separation  = {sep  * 1000:.0f} µm')
#     print(f'  {tag} valley      = {val_amp:.1f} dB  ->  {res_str}')


# print('\n-- Drawing lateral panel --')
# draw_profile_panel(
#     ax_lat, x_mm, lat_prof, lat_p1, lat_p2,
#     xlabel='Lateral (mm)',
#     title=f'Lateral profile  at  z = {tz_mm:.1f} mm',
# )

# print('\n-- Drawing axial panel --')
# draw_profile_panel(
#     ax_ax, z_mm, ax_prof, ax_p1, ax_p2,
#     xlabel='Depth (mm)',
#     title=f'Axial profile  at  x = {tx_mm:.1f} mm',
#     is_axial=True,
# )

# # Panel labels (a), (b), (c)
# for ax, lbl in [(ax_bm, '(a)'), (ax_lat, '(b)'), (ax_ax, '(c)')]:
#     ax.text(-0.09, 1.02, lbl, transform=ax.transAxes,
#             fontsize=11, weight='bold', color='#1A1A1A', va='bottom')

# # ─────────────────────────────────────────────────────────────────────────────
# # 6.  Save
# # ─────────────────────────────────────────────────────────────────────────────
# os.makedirs(OUT_DIR, exist_ok=True)
# stem = os.path.join(OUT_DIR, 'resolution_figure_v3')
# fig.savefig(stem + '.pdf', bbox_inches='tight')
# fig.savefig(stem + '.png', dpi=300, bbox_inches='tight')
# plt.show()
# print(f'\nSaved -> {stem}.pdf  /  {stem}.png')



# ------------------------------------------------------------------------------------------------


"""
resolution_figure.py  (v3)
--------------------------
Changes from v2
---------------
  - Profile extraction uses a ±PROFILE_MARGIN_MM depth/lateral window and
    takes the MAX over that window at each position — handles tilted phantoms
  - Profiles are normalised to their own peak (0 dB = local max), not the
    global image peak (which is biased by surface artefacts)
  - Lateral panel (b) shows ALL peaks in the full slice so each wire target
    is visible; the two closest peaks are highlighted; slight zoom-out so
    peaks stand out against the noise floor
  - Layout: right column wider, spacing tightened

MATLAB export block (add to end of your reconstruction script)
--------------------------------------------------------------
  target_z_mm =  30.0;
  target_x_mm =   2.0;
  [~, iz] = min(abs(z_axis*1000 - target_z_mm));
  [~, ix] = min(abs(x_axis*1000 - target_x_mm));
  export = struct();
  export.B_Mode_dB       = B_Mode_dB;
  export.x_axis_mm       = x_axis * 1000;
  export.z_axis_mm       = z_axis * 1000;
  export.target_iz       = iz;
  export.target_ix       = ix;
  export.target_z_mm     = target_z_mm;
  export.target_x_mm     = target_x_mm;
  save('resolution_export.mat', '-struct', 'export');

Usage
-----
  python resolution_figure.py
  python resolution_figure.py path/to/export.mat

Figure size presets
-------------------
  Half A4   : (13.5, 6.5)   <- default
  Two-thirds: (13.5, 8.2)
  Full col  : ( 8.5, 5.0)

Requirements:  pip install scipy numpy matplotlib
"""

import sys, os
import numpy as np
import scipy.io as sio
import scipy.signal as sig
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from   matplotlib.ticker import MultipleLocator, AutoMinorLocator

# ─────────────────────────────────────────────────────────────────────────────
# 0.  Configuration
# ─────────────────────────────────────────────────────────────────────────────
MAT_FILE  = sys.argv[1] if len(sys.argv) > 1 else 'resolution_export.mat'
OUT_DIR   = 'out'

FIGSIZE   = (13.5, 6.5)   # (width, height) in inches  <- edit to resize

DB_RANGE  = (-50, 0)      # dynamic range for B-mode display
CMAP      = 'inferno'     # B-mode colourmap

DB_THRESH = -6            # dB level for FWHM  (-6 dB = half-power)

# ── Profile extraction ────────────────────────────────────────────────────────
# Window half-width (mm) around the chosen slice line.
# At each lateral/axial position, the MAXIMUM over this window is taken,
# so slightly off-centre targets are captured correctly.
# Increase if your phantom tilt is large (try 1.0 → 2.0).
PROFILE_MARGIN_MM = 0.5   # ±mm around target_z / target_x

# ── Peak detection ────────────────────────────────────────────────────────────
# PEAK_MIN_HEIGHT  : absolute dB floor (relative to the per-slice peak = 0 dB).
#   After per-slice normalisation the real targets sit near 0 dB and noise
#   is far below.  Start with -20 and raise if noise peaks are still picked.
PEAK_MIN_HEIGHT     = -20  # dB relative to per-slice peak

# PEAK_MIN_PROMINENCE : peak must rise this many dB above its neighbours.
# PEAK_MIN_DISTANCE   : minimum gap between peaks in samples.
PEAK_MIN_PROMINENCE = 3
PEAK_MIN_DISTANCE   = 3

# Extra padding on EACH SIDE beyond the outermost detected peak, for the
# lateral panel (b).  Set to 0 to trim tightly.
LATERAL_PAD_MM = 1.5

# Zoom padding for the axial panel (c) around the two closest peaks.
AXIAL_ZOOM_PAD_MM = 3.0

PHANTOM_LABEL = 'CIRS 040GSE — resolution group'

# ─────────────────────────────────────────────────────────────────────────────
# 1.  Style
# ─────────────────────────────────────────────────────────────────────────────
mpl.rcParams.update({
    'font.family':      'serif',
    'font.serif':       ['DejaVu Serif', 'Times New Roman', 'Times'],
    'font.size':         9.5,
    'axes.labelsize':    9.5,
    'axes.titlesize':    10,
    'xtick.labelsize':   8.5,
    'ytick.labelsize':   8.5,
    'axes.linewidth':    0.7,
    'xtick.major.width': 0.7,
    'ytick.major.width': 0.7,
    'xtick.major.size':  3.5,
    'ytick.major.size':  3.5,
    'xtick.direction':  'in',
    'ytick.direction':  'in',
    'mathtext.fontset': 'cm',
    'pdf.fonttype':      42,
    'ps.fonttype':       42,
})

C_CROSS      = '#E05A2B'
C_VALLEY     = '#444444'
C_PROF       = '#1A1A1A'
PEAK_COLORS  = ['#1158A0', '#C8580A']
FILL_COLORS  = ['#C5DEF0', '#FAE0C5']

# ─────────────────────────────────────────────────────────────────────────────
# 2.  Load & extract profiles
# ─────────────────────────────────────────────────────────────────────────────
print(f'Loading: {MAT_FILE}')
d    = sio.loadmat(MAT_FILE, squeeze_me=True, struct_as_record=False)

B_dB  = np.asarray(d['B_Mode_dB'],   dtype=float)
x_mm  = np.asarray(d['x_axis_mm'],   dtype=float).ravel()
z_mm  = np.asarray(d['z_axis_mm'],   dtype=float).ravel()
iz    = int(d['target_iz']) - 1       # MATLAB 1-based → Python 0-based
ix    = int(d['target_ix']) - 1
tz_mm = float(d['target_z_mm'])
tx_mm = float(d['target_x_mm'])

print(f'  B-mode: {B_dB.shape}   slice: x={tx_mm:.1f} mm, z={tz_mm:.1f} mm')

# ── Profile extraction with margin ───────────────────────────────────────────
# For the LATERAL profile: take max over a depth window of ±PROFILE_MARGIN_MM
iz_lo = max(0, np.argmin(np.abs(z_mm - (tz_mm - PROFILE_MARGIN_MM))))
iz_hi = min(B_dB.shape[0] - 1,
            np.argmin(np.abs(z_mm - (tz_mm + PROFILE_MARGIN_MM))))
lat_prof_raw = np.max(B_dB[iz_lo : iz_hi + 1, :], axis=0)

# For the AXIAL profile: take max over a lateral window of ±PROFILE_MARGIN_MM
ix_lo = max(0, np.argmin(np.abs(x_mm - (tx_mm - PROFILE_MARGIN_MM))))
ix_hi = min(B_dB.shape[1] - 1,
            np.argmin(np.abs(x_mm - (tx_mm + PROFILE_MARGIN_MM))))
ax_prof_raw = np.max(B_dB[:, ix_lo : ix_hi + 1], axis=1)

# ── Per-slice normalisation (0 dB = this slice's own peak) ───────────────────
lat_prof = lat_prof_raw - np.max(lat_prof_raw)
ax_prof  = ax_prof_raw  - np.max(ax_prof_raw)

print(f'  Lateral window: iz {iz_lo}–{iz_hi}  '
      f'(z = {z_mm[iz_lo]:.1f}–{z_mm[iz_hi]:.1f} mm)')
print(f'  Axial   window: ix {ix_lo}–{ix_hi}  '
      f'(x = {x_mm[ix_lo]:.1f}–{x_mm[ix_hi]:.1f} mm)')

# ─────────────────────────────────────────────────────────────────────────────
# 3.  Peak-finding helpers
# ─────────────────────────────────────────────────────────────────────────────

def find_all_peaks(axis, profile):
    """Return indices of all detected peaks above the absolute height floor."""
    peaks, _ = sig.find_peaks(
        profile,
        height=PEAK_MIN_HEIGHT,
        prominence=PEAK_MIN_PROMINENCE,
        distance=PEAK_MIN_DISTANCE,
    )
    if len(peaks) == 0:
        print(f'  WARNING: no peaks found above PEAK_MIN_HEIGHT={PEAK_MIN_HEIGHT} dB. '
              f'Global max of profile = {profile.max():.1f} dB. '
              f'Lower PEAK_MIN_HEIGHT.')
    return peaks


def two_closest(peaks, axis):
    """From a list of peak indices, return the pair with smallest axis separation."""
    if len(peaks) < 2:
        return None, None
    best_gap, best_pair = np.inf, (peaks[0], peaks[1])
    for i in range(len(peaks)):
        for j in range(i + 1, len(peaks)):
            gap = abs(axis[peaks[i]] - axis[peaks[j]])
            if gap < best_gap:
                best_gap, best_pair = gap, (peaks[i], peaks[j])
    return best_pair


def fwhm_of_peak(axis, profile, peak_idx):
    """–6 dB full-width by linear interpolation at threshold crossings."""
    level = profile[peak_idx] + DB_THRESH

    left_seg   = profile[:peak_idx + 1]
    below_left = np.where(left_seg <= level)[0]
    if len(below_left) == 0:
        left = axis[0]
    else:
        j = below_left[-1]
        if j + 1 <= peak_idx:
            dy = profile[j+1] - profile[j]
            t  = (level - profile[j]) / dy if dy != 0 else 0.5
            left = axis[j] + t * (axis[j+1] - axis[j])
        else:
            left = axis[j]

    right_seg   = profile[peak_idx:]
    below_right = np.where(right_seg <= level)[0]
    if len(below_right) == 0:
        right = axis[-1]
    else:
        j = below_right[0] + peak_idx
        if j - 1 >= 0:
            dy = profile[j] - profile[j-1]
            t  = (level - profile[j-1]) / dy if dy != 0 else 0.5
            right = axis[j-1] + t * (axis[j] - axis[j-1])
        else:
            right = axis[j]

    return right - left, left, right, level


def valley_between(axis, profile, pa, pb):
    lo, hi = sorted([pa, pb])
    vi = lo + int(np.argmin(profile[lo : hi + 1]))
    return axis[vi], profile[vi]


# ─────────────────────────────────────────────────────────────────────────────
# 4.  Run analysis
# ─────────────────────────────────────────────────────────────────────────────
print('\n-- Lateral profile --')
lat_peaks        = find_all_peaks(x_mm, lat_prof)
lat_p1, lat_p2   = two_closest(lat_peaks, x_mm)

print('\n-- Axial profile --')
ax_peaks         = find_all_peaks(z_mm, ax_prof)
ax_p1, ax_p2     = two_closest(ax_peaks, z_mm)

# ─────────────────────────────────────────────────────────────────────────────
# 5.  Figure layout
# ─────────────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=FIGSIZE)
gs  = gridspec.GridSpec(
    2, 2,
    width_ratios=[0.85, 1.30],   # right (profile) column is wider
    height_ratios=[1, 1],
    hspace=0.48, wspace=0.32,    # tighter horizontal gap
    left=0.07, right=0.97, top=0.93, bottom=0.10,
)
ax_bm  = fig.add_subplot(gs[:, 0])
ax_lat = fig.add_subplot(gs[0, 1])
ax_ax  = fig.add_subplot(gs[1, 1])


# ── 5a. B-mode ────────────────────────────────────────────────────────────────
im = ax_bm.imshow(
    B_dB,
    aspect='auto',
    extent=[x_mm[0], x_mm[-1], z_mm[-1], z_mm[0]],
    vmin=DB_RANGE[0], vmax=DB_RANGE[1],
    cmap=CMAP, interpolation='bilinear',
)
cbar = fig.colorbar(im, ax=ax_bm, fraction=0.046, pad=0.03, aspect=38)
cbar.set_label('Normalised amplitude (dB)', fontsize=8.5)
cbar.ax.tick_params(labelsize=8)
cbar.set_ticks(range(DB_RANGE[0], DB_RANGE[1] + 1, 10))

# Shaded extraction window on B-mode (shows the ±margin band visually)
ax_bm.axhspan(tz_mm - PROFILE_MARGIN_MM, tz_mm + PROFILE_MARGIN_MM,
              color=C_CROSS, alpha=0.18, zorder=2)
ax_bm.axvspan(tx_mm - PROFILE_MARGIN_MM, tx_mm + PROFILE_MARGIN_MM,
              color=PEAK_COLORS[1], alpha=0.12, zorder=2)

# Centre lines
ax_bm.axhline(tz_mm, color=C_CROSS,        lw=0.9, ls='--', alpha=0.9,
              label=f'Lateral slice  z = {tz_mm:.0f} mm  (±{PROFILE_MARGIN_MM} mm)')
ax_bm.axvline(tx_mm, color=PEAK_COLORS[1], lw=0.9, ls=':',  alpha=0.9,
              label=f'Axial slice  x = {tx_mm:.1f} mm  (±{PROFILE_MARGIN_MM} mm)')
ax_bm.plot(tx_mm, tz_mm, 'o',
           color='none', mec=C_CROSS, mew=1.3, ms=10, zorder=5)

ax_bm.set_xlabel('Lateral (mm)', labelpad=3)
ax_bm.set_ylabel('Axial / Depth (mm)', labelpad=3)
ax_bm.set_title('B-mode reconstruction — DAS + CF', pad=5)
ax_bm.legend(fontsize=7, framealpha=0.45, loc='lower right',
             facecolor='black', labelcolor='white', edgecolor='none')
if PHANTOM_LABEL:
    ax_bm.text(0.02, 0.98, PHANTOM_LABEL, transform=ax_bm.transAxes,
               fontsize=7.5, color='white', va='top', ha='left',
               bbox=dict(boxstyle='round,pad=0.25', fc='#1A1A1A', alpha=0.55, ec='none'))


# ─────────────────────────────────────────────────────────────────────────────
# 6.  Profile panel drawing function
# ─────────────────────────────────────────────────────────────────────────────

def draw_profile(ax, axis, profile, all_peaks, p1_idx, p2_idx,
                 xlabel, title, zoom_to_two=False, pad_mm=2.0):
    """
    Draw a profile panel.

    Parameters
    ----------
    zoom_to_two : bool
        True  → zoom around the two closest peaks only  (axial panel)
        False → show all detected peaks in full          (lateral panel)
    pad_mm : float
        Extra padding on each side of the visible region.
    """
    dy_span = 2 - DB_RANGE[0]

    # ── X limits ──────────────────────────────────────────────────────────────
    if zoom_to_two and p1_idx is not None:
        if axis[p1_idx] > axis[p2_idx]:
            p1_idx, p2_idx = p2_idx, p1_idx
        _, l1, r1, _ = fwhm_of_peak(axis, profile, p1_idx)
        _, l2, r2, _ = fwhm_of_peak(axis, profile, p2_idx)
        x_lo = min(l1, l2) - pad_mm
        x_hi = max(r1, r2) + pad_mm
    else:
        # Show all detected peaks (plus a bit of margin)
        if len(all_peaks) > 0:
            x_lo = axis[all_peaks[0]]  - pad_mm
            x_hi = axis[all_peaks[-1]] + pad_mm
        else:
            x_lo, x_hi = axis[0], axis[-1]
        # Clamp to data range
        x_lo = max(x_lo, axis[0])
        x_hi = min(x_hi, axis[-1])

    # ── Highlight the two closest peaks with coloured bands ───────────────────
    if p1_idx is not None:
        if axis[p1_idx] > axis[p2_idx]:
            p1_idx, p2_idx = p2_idx, p1_idx
        fwhm1, l1, r1, lev1 = fwhm_of_peak(axis, profile, p1_idx)
        fwhm2, l2, r2, lev2 = fwhm_of_peak(axis, profile, p2_idx)
        val_pos, val_amp     = valley_between(axis, profile, p1_idx, p2_idx)
        is_resolved = (val_amp < lev1) and (val_amp < lev2)
        res_str = 'Resolved' if is_resolved else 'Unresolved'
        res_col = '#1C7A3A' if is_resolved else '#AA2222'

        for (l, r, lev, pc, fc) in [
            (l1, r1, lev1, PEAK_COLORS[0], FILL_COLORS[0]),
            (l2, r2, lev2, PEAK_COLORS[1], FILL_COLORS[1]),
        ]:
            band_mask = (axis >= l) & (axis <= r) & (profile >= lev)
            ax.fill_between(axis, DB_RANGE[0], profile,
                            where=band_mask, color=fc, alpha=0.80, zorder=2)
            ax.vlines([l, r], DB_RANGE[0], lev,
                      color=pc, lw=0.85, ls=':', alpha=0.65, zorder=3)

    # ── Profile curve ─────────────────────────────────────────────────────────
    ax.plot(axis, profile, color=C_PROF, lw=1.5, zorder=5)

    # ── Mark ALL detected peaks with small triangles ──────────────────────────
    for pk in all_peaks:
        ax.plot(axis[pk], profile[pk], 'v',
                color='#555555', ms=4.5, zorder=6, clip_on=True)

    if p1_idx is not None:
        # –6 dB threshold lines (only for the two highlighted peaks)
        for lev, pc in [(lev1, PEAK_COLORS[0]), (lev2, PEAK_COLORS[1])]:
            ax.hlines(lev, x_lo, x_hi, color=pc, lw=0.85, ls='--',
                      alpha=0.60, zorder=3)

        # FWHM arrows (staggered so labels don't overlap)
        for (l, r, fwhm, pc), frac in zip(
            [(l1, r1, fwhm1, PEAK_COLORS[0]),
             (l2, r2, fwhm2, PEAK_COLORS[1])],
            [0.07, 0.15],
        ):
            y_arr = DB_RANGE[0] + dy_span * frac
            ax.annotate('', xy=(r, y_arr), xytext=(l, y_arr),
                        arrowprops=dict(arrowstyle='<->', color=pc, lw=1.3,
                                        shrinkA=0, shrinkB=0),
                        zorder=6)
            ax.text((l + r) / 2, y_arr - dy_span * 0.025,
                    f'{fwhm * 1000:.0f} µm',
                    ha='center', va='top', fontsize=8.5,
                    color=pc, weight='bold', zorder=7)

        # Valley marker
        ax.plot(val_pos, val_amp, 'v', color=C_VALLEY, ms=7, zorder=6)
        ax.text(val_pos, val_amp - dy_span * 0.03,
                f'{val_amp:.1f} dB  —  {res_str}',
                ha='center', va='top', fontsize=8, color=res_col,
                weight='bold', zorder=7,
                bbox=dict(boxstyle='round,pad=0.25', fc='white',
                          alpha=0.80, ec=res_col, lw=0.6))

        # –6 dB labels outside right spine
        for lev, pc in [(lev1, PEAK_COLORS[0]), (lev2, PEAK_COLORS[1])]:
            ax.text(x_hi + 0.01 * (x_hi - x_lo), lev,
                    f'{DB_THRESH} dB',
                    ha='left', va='center', fontsize=7.5,
                    color=pc, clip_on=False)

        # Console output
        tag = 'Axial' if zoom_to_two else 'Lateral'
        sep = abs(axis[p2_idx] - axis[p1_idx])
        print(f'  {tag} FWHM peak-1 = {fwhm1 * 1000:.0f} µm')
        print(f'  {tag} FWHM peak-2 = {fwhm2 * 1000:.0f} µm')
        print(f'  {tag} separation  = {sep  * 1000:.0f} µm')
        print(f'  {tag} valley      = {val_amp:.1f} dB  ->  {res_str}')

    # ── Axis formatting ───────────────────────────────────────────────────────
    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(DB_RANGE[0], 2)
    ax.set_xlabel(xlabel, labelpad=2)
    ax.set_ylabel('Amplitude (dB, per-slice norm.)', labelpad=2)
    ax.set_title(title, pad=4)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(MultipleLocator(5))
    ax.tick_params(which='minor', length=1.8, width=0.5)

    # Horizontal reference at 0 dB (per-slice peak level)
    ax.axhline(0, color='#888888', lw=0.5, ls='-', alpha=0.4, zorder=1)


# ─────────────────────────────────────────────────────────────────────────────
# 7.  Draw panels
# ─────────────────────────────────────────────────────────────────────────────
print('\n-- Drawing lateral panel (all peaks shown) --')
draw_profile(
    ax_lat, x_mm, lat_prof, lat_peaks, lat_p1, lat_p2,
    xlabel='Lateral (mm)',
    title=f'Lateral profile  at  z = {tz_mm:.1f} mm  (±{PROFILE_MARGIN_MM} mm max projection)',
    zoom_to_two=False,
    pad_mm=LATERAL_PAD_MM,
)

print('\n-- Drawing axial panel (zoom to two closest peaks) --')
draw_profile(
    ax_ax, z_mm, ax_prof, ax_peaks, ax_p1, ax_p2,
    xlabel='Depth (mm)',
    title=f'Axial profile  at  x = {tx_mm:.1f} mm  (±{PROFILE_MARGIN_MM} mm max projection)',
    zoom_to_two=True,
    pad_mm=AXIAL_ZOOM_PAD_MM,
)

# Panel labels
for ax, lbl in [(ax_bm, '(a)'), (ax_lat, '(b)'), (ax_ax, '(c)')]:
    ax.text(-0.08, 1.02, lbl, transform=ax.transAxes,
            fontsize=11, weight='bold', color='#1A1A1A', va='bottom')

# ─────────────────────────────────────────────────────────────────────────────
# 8.  Save
# ─────────────────────────────────────────────────────────────────────────────
os.makedirs(OUT_DIR, exist_ok=True)
stem = os.path.join(OUT_DIR, 'resolution_figure_v4')
fig.savefig(stem + '.pdf', bbox_inches='tight')
fig.savefig(stem + '.png', dpi=300, bbox_inches='tight')
plt.show()
print(f'\nSaved -> {stem}.pdf  /  {stem}.png')