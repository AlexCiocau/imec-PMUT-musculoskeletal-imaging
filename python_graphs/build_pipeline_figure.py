# """
# Publication-quality figure: hydrophone-voltage to calibrated-pressure
# signal-processing pipeline (thesis Fig. 4.5).

# Outputs two files into ./out/:
#   - pressure_pipeline.pdf / .png            (clean, no caption — for LaTeX)
#   - pressure_pipeline_with_caption.pdf/.png (self-contained, for slides/drafts)

# == HOW TO PROVIDE THE DATA ============================================

# You need 6 numpy arrays (every "*_raw" / "*_rec" pair must share length):

#     t_raw, v_raw         time-domain hydrophone voltage  (us, V)
#     t_env_raw, v_env_raw envelope of v_raw               (us, V)
#     f_v, V_f             |FFT(v_raw)|                    (MHz, V)
#     f_p, P_f             |FFT(p_rec)| = |V_f / S(f)|     (MHz, Pa)
#     f_s, S_f             hydrophone sensitivity curve    (MHz, dB)
#     t_rec, p_rec         reconstructed pressure waveform (us, Pa)
#     t_env_rec, p_env_rec envelope of p_rec               (us, Pa)

# Pick whichever of the three loaders below matches your situation.
# """
# import os
# import numpy as np
# import matplotlib as mpl
# import matplotlib.pyplot as plt
# from matplotlib.patches import FancyArrowPatch


# # =====================================================================
# # DATA LOADING -- pick ONE of (A) / (B) / (C)
# # =====================================================================

# def load_data():
#     # ---- (A) RECOMMENDED: from a single .mat file you exported in MATLAB
#     # In MATLAB do (just once):
#     #   save('pipeline_data.mat', ...
#     #        't_raw','v_raw','t_env_raw','v_env_raw', ...
#     #        'f_v','V_f','f_p','P_f','f_s','S_f', ...
#     #        't_rec','p_rec','t_env_rec','p_env_rec', '-v7')
#     # Then uncomment:
#     #
#     # from scipy.io import loadmat
#     # m = loadmat('pipeline_data.mat', squeeze_me=True)
#     # return {k: np.asarray(m[k], dtype=float) for k in
#     #         ['t_raw','v_raw','t_env_raw','v_env_raw',
#     #          'f_v','V_f','f_p','P_f','f_s','S_f',
#     #          't_rec','p_rec','t_env_rec','p_env_rec']}

#     # ---- (B) From CSV files (portable, no scipy needed for .mat). Each
#     # CSV has two columns: x, y  (no header). File names self-explanatory.
#     #
#     # def csv(path): return np.loadtxt(path, delimiter=',').T
#     # t_raw,     v_raw     = csv('raw_signal.csv')
#     # t_env_raw, v_env_raw = csv('raw_signal_envelope.csv')
#     # f_v, V_f             = csv('voltage_fft.csv')
#     # f_p, P_f             = csv('pressure_fft.csv')
#     # f_s, S_f             = csv('sensitivity.csv')
#     # t_rec,     p_rec     = csv('reconstructed_pressure.csv')
#     # t_env_rec, p_env_rec = csv('reconstructed_pressure_envelope.csv')
#     # return dict(t_raw=t_raw, v_raw=v_raw, ...)   # fill in

#     # ---- (C) From the pickle produced by extract_fig_data.py
#     # (only if you can't re-run MATLAB and have to parse the .fig files).
#     import pickle
#     with open('extracted_data.pkl', 'rb') as f:
#         E = pickle.load(f)
#     # Each entry is a list of (offset, length, array) runs ordered by
#     # byte offset in the .fig blob; for time-domain plots the order is
#     # (env_x, env_y, signal_x, signal_y); for FFT plots it is (x, y).
#     return dict(
#         t_raw     = E['RawSignal'][2][2],
#         v_raw     = E['RawSignal'][3][2],
#         t_env_raw = E['RawSignal'][0][2],
#         v_env_raw = E['RawSignal'][1][2],
#         f_v       = E['VoltageFFT'][0][2],
#         V_f       = E['VoltageFFT'][1][2],
#         f_p       = E['PressureFFT'][0][2],
#         P_f       = E['PressureFFT'][1][2],
#         f_s       = E['MarioSensitivityCurve'][0][2],
#         S_f       = E['MarioSensitivityCurve'][1][2],
#         t_rec     = E['Reconstructed_TimeDomain_Pressure'][2][2],
#         p_rec     = E['Reconstructed_TimeDomain_Pressure'][3][2],
#         t_env_rec = E['Reconstructed_TimeDomain_Pressure'][0][2],
#         p_env_rec = E['Reconstructed_TimeDomain_Pressure'][1][2],
#     )


# # =====================================================================
# # STYLE
# # =====================================================================

# # Original Params
# # mpl.rcParams.update({
# #     'font.family': 'serif',
# #     'font.serif': ['DejaVu Serif', 'Times New Roman', 'Times'],
# #     'font.size': 10, 'axes.labelsize': 10, 'axes.titlesize': 10,
# #     'xtick.labelsize': 9, 'ytick.labelsize': 9,
# #     'axes.linewidth': 0.8,
# #     'xtick.major.width': 0.8, 'ytick.major.width': 0.8,
# #     'xtick.major.size': 3,    'ytick.major.size': 3,
# #     'xtick.direction': 'in',  'ytick.direction': 'in',
# #     'mathtext.fontset': 'cm',
# #     'axes.spines.top': False, 'axes.spines.right': False,
# #     'pdf.fonttype': 42, 'ps.fonttype': 42,   # editable text in PDF
# # })

# # Maarten's Suggestion
# mpl.rcParams.update({
#     'font.family': 'serif',
#     'font.serif': ['DejaVu Serif', 'Times New Roman', 'Times'],
#     'font.size': 14, 'axes.labelsize': 14, 'axes.titlesize': 14,
#     'xtick.labelsize': 12, 'ytick.labelsize': 12,
#     'axes.linewidth': 0.8,
#     'xtick.major.width': 0.8, 'ytick.major.width': 0.8,
#     'xtick.major.size': 3,    'ytick.major.size': 3,
#     'xtick.direction': 'in',  'ytick.direction': 'in',
#     'mathtext.fontset': 'cm',
#     'axes.spines.top': False, 'axes.spines.right': False,
#     'pdf.fonttype': 42, 'ps.fonttype': 42,   # editable text in PDF
# })

# # Colour palette
# C_VOLT   = '#1f4e79'   # voltage (dark blue)
# C_PRESS  = '#a8323b'   # pressure (deep red)
# C_SENS   = '#2f6f4f'   # sensitivity curve (green)
# C_ACCENT = '#0e5a73'   # operation arrows (teal)
# C_ENV    = '#cdd5dc'   # voltage envelope fill (pale grey-blue)
# C_ENV_P  = '#e6cdd1'   # pressure envelope fill (pale pink)


# # =====================================================================
# # FIGURE BUILDER
# # =====================================================================

# def build_figure(D, include_caption=False):
#     """D: dict of arrays returned by load_data()."""

#     # Plot panel (b)/(c) only show frequencies up to F_MAX
#     F_MAX = 30.0
#     mv = D['f_v'] <= F_MAX
#     mp = D['f_p'] <= F_MAX
#     ms = D['f_s'] <= F_MAX

#     # ---- canvas + manual axes layout -------------------------------
#     fig_h = 10.5 if include_caption else 9.5
#     fig = plt.figure(figsize=(8.0, fig_h))

#     LEFT, RIGHT = 0.045, 0.97
#     TOP = 0.97
#     BOTTOM = 0.105 if include_caption else 0.045

#     W_LABEL = 0.13           # left column: step name
#     W_ARROW = 0.22           # middle column: arrow + sensitivity inset
#     x_label_l, x_label_r = LEFT, LEFT + W_LABEL
#     x_arrow_l, x_arrow_r = x_label_r, x_label_r + W_ARROW
#     x_plot_l,  x_plot_r  = x_arrow_r, RIGHT

#     # vertical spacing: divide step gets a larger gap to fit the inset
#     H_TOTAL = TOP - BOTTOM
#     GAP_FFT, GAP_DIVIDE, GAP_IFFT = 0.058, 0.108, 0.058
#     H_PLOT = (H_TOTAL - GAP_FFT - GAP_DIVIDE - GAP_IFFT) / 4
#     y0 = TOP - H_PLOT
#     y1 = y0 - GAP_FFT    - H_PLOT
#     y2 = y1 - GAP_DIVIDE - H_PLOT
#     y3 = y2 - GAP_IFFT   - H_PLOT
#     W_PLOT = x_plot_r - x_plot_l
#     ax1 = fig.add_axes([x_plot_l, y0, W_PLOT, H_PLOT])
#     ax2 = fig.add_axes([x_plot_l, y1, W_PLOT, H_PLOT])
#     ax3 = fig.add_axes([x_plot_l, y2, W_PLOT, H_PLOT])
#     ax4 = fig.add_axes([x_plot_l, y3, W_PLOT, H_PLOT])

#     # ---- (a) gated toneburst (time domain) -------------------------
#     ax1.fill_between(D['t_env_raw'], D['v_env_raw'], -D['v_env_raw'],
#                      color=C_ENV, alpha=0.6, lw=0)
#     ax1.plot(D['t_raw'], D['v_raw'], color=C_VOLT, lw=0.7)
#     ax1.set_xlim(D['t_raw'].min(), D['t_raw'].max())
#     ax1.set_ylim(-1.85, 1.85)
#     ax1.set_xlabel(r'Time ($\mu$s)', labelpad=2)
#     ax1.set_ylabel('Voltage (V)',    labelpad=4)
#     ax1.text(0.985, 0.93, '(a)', transform=ax1.transAxes,
#              ha='right', va='top', fontsize=10, style='italic', color='#444')

#     # gate-window markers (boundaries inferred from the support of the
#     # reconstructed pressure -- replace if your gate is different)
#     t_gate_lo, t_gate_hi = D['t_rec'].min(), D['t_rec'].max()
#     for tg in (t_gate_lo, t_gate_hi):
#         ax1.axvline(tg, color='#666', ls='--', lw=0.8, alpha=0.85, zorder=2)
#     ax1.text(0.5*(t_gate_lo + t_gate_hi), 1.72, 'gate',
#              ha='center', va='center', fontsize=9, color='#666',
#              style='italic',
#              bbox=dict(facecolor='white', edgecolor='none', pad=1.2))

#     # ---- (b) voltage spectrum --------------------------------------
#     ax2.fill_between(D['f_v'][mv], 0, D['V_f'][mv],
#                      color=C_VOLT, alpha=0.18, lw=0)
#     ax2.plot(D['f_v'][mv], D['V_f'][mv], color=C_VOLT, lw=1.3)
#     ax2.set_xlim(0, F_MAX)
#     ax2.set_xlabel('Frequency (MHz)', labelpad=2)
#     ax2.set_ylabel(r'$|\,\widetilde{V}(f)\,|$  (V)', labelpad=4)
#     ax2.text(0.985, 0.93, '(b)', transform=ax2.transAxes,
#              ha='right', va='top', fontsize=10, style='italic', color='#444')

#     # ---- (c) pressure spectrum -------------------------------------
#     ax3.fill_between(D['f_p'][mp], 0, D['P_f'][mp],
#                      color=C_PRESS, alpha=0.18, lw=0)
#     ax3.plot(D['f_p'][mp], D['P_f'][mp], color=C_PRESS, lw=1.3)
#     ax3.set_xlim(0, F_MAX)
#     ax3.set_xlabel('Frequency (MHz)', labelpad=2)
#     ax3.set_ylabel(r'$|\,\widetilde{P}(f)\,|$  (Pa)', labelpad=4)
#     ax3.text(0.985, 0.93, '(c)', transform=ax3.transAxes,
#              ha='right', va='top', fontsize=10, style='italic', color='#444')

#     # ---- (d) reconstructed pressure (time domain) ------------------
#     ax4.fill_between(D['t_env_rec'], D['p_env_rec'], -D['p_env_rec'],
#                      color=C_ENV_P, alpha=0.6, lw=0)
#     ax4.plot(D['t_rec'], D['p_rec'], color=C_PRESS, lw=1.0)
#     ax4.set_xlim(D['t_rec'].min(), D['t_rec'].max())
#     ax4.set_xlabel(r'Time ($\mu$s)', labelpad=2)
#     ax4.set_ylabel('Pressure (Pa)',  labelpad=4)
#     ax4.text(0.985, 0.93, '(d)', transform=ax4.transAxes,
#              ha='right', va='top', fontsize=10, style='italic', color='#444')

#     # P+ / P- (peak positive and peak negative pressure)
#     P_plus, P_minus = float(D['p_rec'].max()), float(D['p_rec'].min())
#     ax4.set_ylim(P_minus - 6.0, P_plus + 6.0)
#     for y in (P_plus, P_minus):
#         ax4.axhline(y, color='#333', ls=':', lw=0.8, alpha=0.7, zorder=1)
#     x_lab = D['t_rec'].min() + 0.012*(D['t_rec'].max() - D['t_rec'].min())
#     ax4.text(x_lab, P_plus,  r'$\,P_{+}$', ha='left', va='bottom',
#              fontsize=10.5, color='#111')
#     ax4.text(x_lab, P_minus, r'$\,P_{-}$', ha='left', va='top',
#              fontsize=10.5, color='#111')

#     # ---- left-column step labels -----------------------------------
#     step_labels = ['Gated\ntoneburst', 'Voltage\nspectrum',
#                    'Pressure\nspectrum', 'Reconstructed\npressure']
#     x_label_centre = 0.5*(x_label_l + x_label_r)
#     for ax, label in zip([ax1, ax2, ax3, ax4], step_labels):
#         bb = ax.get_position()
#         fig.text(x_label_centre, bb.y0 + bb.height/2, label,
#                  ha='center', va='center', fontsize=11.5, weight='bold',
#                  color='#222')

#     # ---- operation arrows + sensitivity inset ----------------------
#     def draw_op(ax_top, ax_bot, label, inset=False, label_offset_x=0.045):
#         bb_t, bb_b = ax_top.get_position(), ax_bot.get_position()
#         y_top  = bb_t.y0 - 0.005
#         y_bot  = bb_b.y1 + 0.005
#         y_mid  = 0.5*(y_top + y_bot)
#         x_arrow = x_arrow_l + 0.22*(x_arrow_r - x_arrow_l)

#         fig.patches.append(FancyArrowPatch(
#             (x_arrow, y_top), (x_arrow, y_bot),
#             transform=fig.transFigure, arrowstyle='-|>',
#             mutation_scale=20, lw=2.4, color=C_ACCENT))

#         fig.text(x_arrow - label_offset_x, y_mid, label,
#                  ha='center', va='center', fontsize=12.5,
#                  weight='bold', color=C_ACCENT)

#         if inset:
#             # Original
#             # ins_x = x_arrow + 0.040
#             # ins_w = (x_arrow_r - ins_x) - 0.005
#             # ins_h = (y_top - y_bot) * 0.66
#             # ins_y = y_mid - ins_h/2 - 0.003

#             # Maarten's
#             ins_x = 0.18                    # left edge of inset (smaller = further left)
#             ins_w = 0.22                    # width
#             ins_h = gap_h * 0.95            # height (fraction of the divide-row gap)
#             ins_y = y_mid - ins_h/2         # vertical centre on the arrow midpoint

#             ax_ins = fig.add_axes([ins_x, ins_y, ins_w, ins_h])
#             ax_ins.plot(D['f_s'][ms], D['S_f'][ms], color=C_SENS, lw=1.2)
#             ax_ins.fill_between(D['f_s'][ms], D['S_f'].min()-2,
#                                 D['S_f'][ms], color=C_SENS, alpha=0.10, lw=0)
#             ax_ins.set_facecolor('#f8fbf6')
#             ax_ins.tick_params(labelsize=7, length=1.8, pad=1)
#             ax_ins.set_ylabel(r'$S$ (dB)', fontsize=7.5, labelpad=1)
#             ax_ins.set_xlim(D['f_s'][ms].min(), D['f_s'][ms].max())
#             ax_ins.set_ylim(D['S_f'][ms].min()-1, D['S_f'][ms].max()+2)
#             ax_ins.set_title('Sensitivity curve', fontsize=8, pad=2,
#                              color=C_SENS, weight='bold')
#             ax_ins.text(0.97, 0.06, r'$f$ (MHz)',
#                         transform=ax_ins.transAxes,
#                         ha='right', va='bottom', fontsize=7,
#                         color='#666', style='italic')
#             for sp in ax_ins.spines.values():
#                 sp.set_color('#999'); sp.set_linewidth(0.6)

#     draw_op(ax1, ax2, 'FFT')
#     draw_op(ax2, ax3, r'$\div\, S(f)$', inset=True, label_offset_x=0.038)
#     draw_op(ax3, ax4, 'IFFT')

#     # ---- optional embedded caption ---------------------------------
#     if include_caption:
#         fig.add_artist(plt.Line2D(
#             [LEFT, RIGHT], [BOTTOM - 0.025]*2,
#             transform=fig.transFigure, color='#cccccc', lw=0.5))
#         caption = (
#             r'$\bf{Figure\ 4.5:}$  '
#             'Signal-processing pipeline for reconstructing absolute '
#             'time-domain pressure waveforms.  The gated\n'
#             'hydrophone voltage (a) is transformed to the frequency '
#             'domain via FFT (b), normalized by the frequency-\n'
#             'dependent hydrophone sensitivity $S(f)$ to obtain the '
#             'pressure spectrum (c), and inverse-transformed (IFFT)\n'
#             'to recover the calibrated acoustic pressure waveform (d).'
#         )
#         fig.text(LEFT, 0.005, caption, ha='left', va='bottom',
#                  fontsize=8.7, color='#222')
#     return fig


# # =====================================================================
# # MAIN
# # =====================================================================

# if __name__ == '__main__':
#     D = load_data()
#     os.makedirs('out', exist_ok=True)

#     fig = build_figure(D, include_caption=False)
#     fig.savefig('out/pressure_pipeline.pdf', bbox_inches='tight')
#     fig.savefig('out/pressure_pipeline.png', dpi=300, bbox_inches='tight')
#     plt.close(fig)

#     fig = build_figure(D, include_caption=True)
#     fig.savefig('out/pressure_pipeline_with_caption.pdf', bbox_inches='tight')
#     fig.savefig('out/pressure_pipeline_with_caption.png',
#                 dpi=300, bbox_inches='tight')
#     plt.close(fig)

#     print('Wrote figures to ./out/')

















"""
Publication-quality figure: hydrophone-voltage to calibrated-pressure
signal-processing pipeline (thesis Fig. 4.5).

Outputs two files into ./out/:
  - pressure_pipeline.pdf / .png            (clean, no caption -- for LaTeX)
  - pressure_pipeline_with_caption.pdf/.png (self-contained, for slides/drafts)

== HOW TO PROVIDE THE DATA ============================================

You need 6 numpy arrays (every "*_raw" / "*_rec" pair must share length):

    t_raw, v_raw         time-domain hydrophone voltage  (us, V)
    t_env_raw, v_env_raw envelope of v_raw               (us, V)
    f_v, V_f             |FFT(v_raw)|                    (MHz, V)
    f_p, P_f             |FFT(p_rec)| = |V_f / S(f)|     (MHz, Pa)
    f_s, S_f             hydrophone sensitivity curve    (MHz, dB)
    t_rec, p_rec         reconstructed pressure waveform (us, Pa)
    t_env_rec, p_env_rec envelope of p_rec               (us, Pa)

Pick whichever of the three loaders below matches your situation.
"""
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch


# =====================================================================
# DATA LOADING -- pick ONE of (A) / (B) / (C)
# =====================================================================

def load_data():
    # ---- (A) RECOMMENDED: from a single .mat file you exported in MATLAB
    # In MATLAB do (just once):
    #   save('pipeline_data.mat', ...
    #        't_raw','v_raw','t_env_raw','v_env_raw', ...
    #        'f_v','V_f','f_p','P_f','f_s','S_f', ...
    #        't_rec','p_rec','t_env_rec','p_env_rec', '-v7')
    # Then uncomment:
    #
    # from scipy.io import loadmat
    # m = loadmat('pipeline_data.mat', squeeze_me=True)
    # return {k: np.asarray(m[k], dtype=float) for k in
    #         ['t_raw','v_raw','t_env_raw','v_env_raw',
    #          'f_v','V_f','f_p','P_f','f_s','S_f',
    #          't_rec','p_rec','t_env_rec','p_env_rec']}

    # ---- (B) From CSV files (portable, no scipy needed for .mat). Each
    # CSV has two columns: x, y  (no header). File names self-explanatory.
    #
    # def csv(path): return np.loadtxt(path, delimiter=',').T
    # t_raw,     v_raw     = csv('raw_signal.csv')
    # t_env_raw, v_env_raw = csv('raw_signal_envelope.csv')
    # f_v, V_f             = csv('voltage_fft.csv')
    # f_p, P_f             = csv('pressure_fft.csv')
    # f_s, S_f             = csv('sensitivity.csv')
    # t_rec,     p_rec     = csv('reconstructed_pressure.csv')
    # t_env_rec, p_env_rec = csv('reconstructed_pressure_envelope.csv')
    # return dict(t_raw=t_raw, v_raw=v_raw, ...)   # fill in

    # ---- (C) From the pickle produced by extract_fig_data.py
    # (only if you can't re-run MATLAB and have to parse the .fig files).
    import pickle
    with open('extracted_data.pkl', 'rb') as f:
        E = pickle.load(f)
    # Each entry is a list of (offset, length, array) runs ordered by
    # byte offset in the .fig blob; for time-domain plots the order is
    # (env_x, env_y, signal_x, signal_y); for FFT plots it is (x, y).
    return dict(
        t_raw     = E['RawSignal'][2][2],
        v_raw     = E['RawSignal'][3][2],
        t_env_raw = E['RawSignal'][0][2],
        v_env_raw = E['RawSignal'][1][2],
        f_v       = E['VoltageFFT'][0][2],
        V_f       = E['VoltageFFT'][1][2],
        f_p       = E['PressureFFT'][0][2],
        P_f       = E['PressureFFT'][1][2],
        f_s       = E['MarioSensitivityCurve'][0][2],
        S_f       = E['MarioSensitivityCurve'][1][2],
        t_rec     = E['Reconstructed_TimeDomain_Pressure'][2][2],
        p_rec     = E['Reconstructed_TimeDomain_Pressure'][3][2],
        t_env_rec = E['Reconstructed_TimeDomain_Pressure'][0][2],
        p_env_rec = E['Reconstructed_TimeDomain_Pressure'][1][2],
    )


# =====================================================================
# STYLE
# =====================================================================

mpl.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['DejaVu Serif', 'Times New Roman', 'Times'],
    'font.size': 14, 'axes.labelsize': 14, 'axes.titlesize': 14,
    'xtick.labelsize': 11, 'ytick.labelsize': 11,
    'axes.linewidth': 0.8,
    'xtick.major.width': 0.8, 'ytick.major.width': 0.8,
    'xtick.major.size': 3,    'ytick.major.size': 3,
    'xtick.direction': 'in',  'ytick.direction': 'in',
    'mathtext.fontset': 'cm',
    'axes.spines.top': False, 'axes.spines.right': False,
    'pdf.fonttype': 42, 'ps.fonttype': 42,   # editable text in PDF
})

# Colour palette
C_VOLT   = '#1f4e79'   # voltage (dark blue)
C_PRESS  = '#a8323b'   # pressure (deep red)
C_SENS   = '#2f6f4f'   # sensitivity curve (green)
C_ACCENT = '#0e5a73'   # operation arrows (teal)
C_ENV    = '#cdd5dc'   # voltage envelope fill (pale grey-blue)
C_ENV_P  = '#e6cdd1'   # pressure envelope fill (pale pink)


# =====================================================================
# FIGURE BUILDER
# =====================================================================

def build_figure(D, include_caption=False):
    """D: dict of arrays returned by load_data()."""

    # Plot panel (b)/(c) only show frequencies up to F_MAX
    F_MAX = 30.0
    mv = D['f_v'] <= F_MAX
    mp = D['f_p'] <= F_MAX
    ms = D['f_s'] <= F_MAX

    # ---- canvas + manual axes layout -------------------------------
    fig_h = 10.5 if include_caption else 9.5
    fig = plt.figure(figsize=(8.0, fig_h))

    LEFT, RIGHT = 0.045, 0.97
    TOP = 0.97
    BOTTOM = 0.105 if include_caption else 0.045

    W_LABEL = 0.13           # left column: step name
    W_ARROW = 0.22           # middle column: arrow + sensitivity inset
    x_label_l, x_label_r = LEFT, LEFT + W_LABEL
    x_arrow_l, x_arrow_r = x_label_r, x_label_r + W_ARROW
    x_plot_l,  x_plot_r  = x_arrow_r, RIGHT

    # vertical spacing: divide step gets a much larger gap so the
    # sensitivity inset can be drawn at a comfortable readable size
    H_TOTAL = TOP - BOTTOM
    GAP_FFT, GAP_DIVIDE, GAP_IFFT = 0.058, 0.18, 0.058
    H_PLOT = (H_TOTAL - GAP_FFT - GAP_DIVIDE - GAP_IFFT) / 4
    y0 = TOP - H_PLOT
    y1 = y0 - GAP_FFT    - H_PLOT
    y2 = y1 - GAP_DIVIDE - H_PLOT
    y3 = y2 - GAP_IFFT   - H_PLOT
    W_PLOT = x_plot_r - x_plot_l
    ax1 = fig.add_axes([x_plot_l, y0, W_PLOT, H_PLOT])
    ax2 = fig.add_axes([x_plot_l, y1, W_PLOT, H_PLOT])
    ax3 = fig.add_axes([x_plot_l, y2, W_PLOT, H_PLOT])
    ax4 = fig.add_axes([x_plot_l, y3, W_PLOT, H_PLOT])

    # ---- (a) gated toneburst (time domain) -------------------------
    ax1.fill_between(D['t_env_raw'], D['v_env_raw'], -D['v_env_raw'],
                     color=C_ENV, alpha=0.6, lw=0)
    ax1.plot(D['t_raw'], D['v_raw'], color=C_VOLT, lw=0.7)
    ax1.set_xlim(D['t_raw'].min(), D['t_raw'].max())
    ax1.set_ylim(-1.85, 1.85)
    ax1.set_xlabel(r'Time ($\mu$s)', labelpad=2)
    ax1.set_ylabel('Voltage (V)',    labelpad=4)
    ax1.text(0.985, 0.93, '(a)', transform=ax1.transAxes,
             ha='right', va='top', fontsize=10, style='italic', color='#444')

    # gate-window markers (boundaries inferred from the support of the
    # reconstructed pressure -- replace if your gate is different)
    t_gate_lo, t_gate_hi = D['t_rec'].min(), D['t_rec'].max()
    for tg in (t_gate_lo, t_gate_hi):
        ax1.axvline(tg, color='#666', ls='--', lw=0.8, alpha=0.85, zorder=2)
    ax1.text(0.5*(t_gate_lo + t_gate_hi), 1.72, 'gate',
             ha='center', va='center', fontsize=9, color='#666',
             style='italic',
             bbox=dict(facecolor='white', edgecolor='none', pad=1.2))

    # ---- (b) voltage spectrum --------------------------------------
    ax2.fill_between(D['f_v'][mv], 0, D['V_f'][mv],
                     color=C_VOLT, alpha=0.18, lw=0)
    ax2.plot(D['f_v'][mv], D['V_f'][mv], color=C_VOLT, lw=1.3)
    ax2.set_xlim(0, F_MAX)
    ax2.set_xlabel('Frequency (MHz)', labelpad=2)
    ax2.set_ylabel(r'$|\,\widetilde{V}(f)\,|$  (V)', labelpad=4)
    ax2.text(0.985, 0.93, '(b)', transform=ax2.transAxes,
             ha='right', va='top', fontsize=10, style='italic', color='#444')

    # ---- (c) pressure spectrum -------------------------------------
    ax3.fill_between(D['f_p'][mp], 0, D['P_f'][mp],
                     color=C_PRESS, alpha=0.18, lw=0)
    ax3.plot(D['f_p'][mp], D['P_f'][mp], color=C_PRESS, lw=1.3)
    ax3.set_xlim(0, F_MAX)
    ax3.set_xlabel('Frequency (MHz)', labelpad=2)
    ax3.set_ylabel(r'$|\,\widetilde{P}(f)\,|$  (Pa)', labelpad=4)
    ax3.text(0.985, 0.93, '(c)', transform=ax3.transAxes,
             ha='right', va='top', fontsize=10, style='italic', color='#444')

    # ---- (d) reconstructed pressure (time domain) ------------------
    ax4.fill_between(D['t_env_rec'], D['p_env_rec'], -D['p_env_rec'],
                     color=C_ENV_P, alpha=0.6, lw=0)
    ax4.plot(D['t_rec'], D['p_rec'], color=C_PRESS, lw=1.0)
    ax4.set_xlim(D['t_rec'].min(), D['t_rec'].max())
    ax4.set_xlabel(r'Time ($\mu$s)', labelpad=2)
    ax4.set_ylabel('Pressure (Pa)',  labelpad=4)
    ax4.text(0.985, 0.93, '(d)', transform=ax4.transAxes,
             ha='right', va='top', fontsize=10, style='italic', color='#444')

    # P+ / P- (peak positive and peak negative pressure)
    P_plus, P_minus = float(D['p_rec'].max()), float(D['p_rec'].min())
    ax4.set_ylim(P_minus - 6.0, P_plus + 6.0)
    for y in (P_plus, P_minus):
        ax4.axhline(y, color='#333', ls=':', lw=0.8, alpha=0.7, zorder=1)
    x_lab = D['t_rec'].min() + 0.012*(D['t_rec'].max() - D['t_rec'].min())
    ax4.text(x_lab, P_plus,  r'$\,P_{+}$', ha='left', va='bottom',
             fontsize=10.5, color='#111')
    ax4.text(x_lab, P_minus, r'$\,P_{-}$', ha='left', va='top',
             fontsize=10.5, color='#111')

    # ---- left-column step labels -----------------------------------
    step_labels = ['Gated\ntoneburst', 'Voltage\nspectrum',
                   'Pressure\nspectrum', 'Reconstructed\npressure']
    # x_label_centre = 0.5*(x_label_l + x_label_r)
    x_label_centre = 0.19
    for ax, label in zip([ax1, ax2, ax3, ax4], step_labels):
        bb = ax.get_position()
        fig.text(x_label_centre, bb.y0 + bb.height/2, label,
                 ha='center', va='center', fontsize=16, weight='bold',
                 color='#222')

    # ---- operation arrows + sensitivity inset ----------------------
    # Layout for the divide-row (between panels b and c):
    #
    #   [ "÷ S(f)" label ]  [ arrow ]  ......  [ sensitivity inset .... ]
    #
    # The inset is now drawn at a hardcoded position so it can be made
    # large and pushed well to the left of where the arrow would naturally
    # sit. The arrow position for the divide step is moved to the right
    # so the inset has free real estate to its left.
    def draw_op(ax_top, ax_bot, label, inset=False, label_offset_x=0.045,
                arrow_pos_frac=0.22):
        bb_t, bb_b = ax_top.get_position(), ax_bot.get_position()
        y_top  = bb_t.y0 - 0.005
        y_bot  = bb_b.y1 + 0.005
        y_mid  = 0.5*(y_top + y_bot)
        # arrow_pos_frac: 0 = left edge of arrow column, 1 = right edge
        x_arrow = x_arrow_l + arrow_pos_frac*(x_arrow_r - x_arrow_l)

        fig.patches.append(FancyArrowPatch(
            (x_arrow, y_top), (x_arrow, y_bot),
            transform=fig.transFigure, arrowstyle='-|>',
            mutation_scale=20, lw=2.4, color=C_ACCENT))

        fig.text(x_arrow - label_offset_x, y_mid, label,
                 ha='center', va='center', fontsize=12.5,
                 weight='bold', color=C_ACCENT)

        if inset:
            # Hardcoded figure-fraction box -- tweak these four numbers
            # to move/resize the sensitivity inset.
            ins_x = 0.03                    # left edge (smaller -> further left)
            ins_y = y_mid - 0.085           # vertical centre on the arrow midpoint
            ins_w = 0.24                    # width
            ins_h = 0.12                    # height
            ax_ins = fig.add_axes([ins_x, ins_y, ins_w, ins_h])
            ax_ins.plot(D['f_s'][ms], D['S_f'][ms], color=C_SENS, lw=1.4)
            ax_ins.fill_between(D['f_s'][ms], D['S_f'].min()-2,
                                D['S_f'][ms], color=C_SENS, alpha=0.10, lw=0)
            ax_ins.set_facecolor('#f8fbf6')
            ax_ins.tick_params(labelsize=9, length=2.5, pad=2)
            ax_ins.set_xlabel(r'Frequency (MHz)', fontsize=10, labelpad=2)
            ax_ins.set_ylabel(r'$S$ (dB)', fontsize=10, labelpad=3)
            ax_ins.set_xlim(D['f_s'][ms].min(), D['f_s'][ms].max())
            ax_ins.set_ylim(D['S_f'][ms].min()-1, D['S_f'][ms].max()+2)
            ax_ins.set_title('Sensitivity curve', fontsize=11, pad=4,
                             color=C_SENS, weight='bold')
            for sp in ax_ins.spines.values():
                sp.set_color('#999'); sp.set_linewidth(0.6)

    draw_op(ax1, ax2, 'FFT')
    # divide step: arrow pushed to the right side of its column, label
    # offset bumped so it sits just left of the arrow (and well clear
    # of the now-large inset on the left)
    draw_op(ax2, ax3, r'$\div\, S(f)$', inset=True,
            label_offset_x=0.045, arrow_pos_frac=0.85)
    draw_op(ax3, ax4, 'IFFT')

    # ---- optional embedded caption ---------------------------------
    if include_caption:
        fig.add_artist(plt.Line2D(
            [LEFT, RIGHT], [BOTTOM - 0.025]*2,
            transform=fig.transFigure, color='#cccccc', lw=0.5))
        caption = (
            r'$\bf{Figure\ 4.5:}$  '
            'Signal-processing pipeline for reconstructing absolute '
            'time-domain pressure waveforms.  The gated\n'
            'hydrophone voltage (a) is transformed to the frequency '
            'domain via FFT (b), normalized by the frequency-\n'
            'dependent hydrophone sensitivity $S(f)$ to obtain the '
            'pressure spectrum (c), and inverse-transformed (IFFT)\n'
            'to recover the calibrated acoustic pressure waveform (d).'
        )
        fig.text(LEFT, 0.005, caption, ha='left', va='bottom',
                 fontsize=8.7, color='#222')
    return fig


# =====================================================================
# MAIN
# =====================================================================

if __name__ == '__main__':
    D = load_data()
    os.makedirs('out', exist_ok=True)

    fig = build_figure(D, include_caption=False)
    fig.savefig('out/pressure_pipeline.pdf', bbox_inches='tight')
    fig.savefig('out/pressure_pipeline.png', dpi=300, bbox_inches='tight')
    plt.close(fig)

    fig = build_figure(D, include_caption=True)
    fig.savefig('out/pressure_pipeline_with_caption.pdf', bbox_inches='tight')
    fig.savefig('out/pressure_pipeline_with_caption.png',
                dpi=300, bbox_inches='tight')
    plt.close(fig)

    print('Wrote figures to ./out/')
