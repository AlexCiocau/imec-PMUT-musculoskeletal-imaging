# PMUT arrays for musculoskeletal ultrasound imaging

**Alexandru Ciocau & Maarten Luypaert**  
Master's Thesis — KU Leuven, in collaboration with imec  
Academic Year 2025–2026

---

## About

This repository contains the code developed as part of our master's thesis at KU Leuven, conducted at imec.



## Repository Structure

```
├── pressure_characteristics/         # compute acoustic pressure from hydrophone data
├── python_graphs/                    # generate publication figures
├── reconstruction_scripts/           # reconstruct B-mode images from RF data
├── miscellaneous/                    # standalone analysis scripts
└── volumetric_simulation/            # simulate and compare PMUT array geometries
    └── auxilliary_sim_files/         # helper scripts for beam visualisation
```

## Requirements

### Python (python_graphs/)
- Python 3.8+
- `numpy`
- `scipy`
- `matplotlib`

Install with:
```bash
pip install numpy scipy matplotlib
```

### MATLAB (all other folders)
- MATLAB R2021a or later
- **Signal Processing Toolbox** — required by `PSF_figure_simple.m`, `Compute_Pressure_Holistic_v4.m`, `L7_4_Recon.m`, `PMUT_coherence_v4.m`
- **Image Processing Toolbox** — required by `PSF_figure_simple.m`, `plot_3D_beam.m`, `plot_beam_profile.m`, `plot_pressure_map.m`, `L7_4_Recon.m`, `PMUT_coherence_v4.m`

> `ChipMUT_3Dsim.m` is fully self-contained and requires no toolboxes.


## Usage

### python_graphs/
All Python scripts are standalone and run directly from their folder. They produce PDF and PNG figures in a local `out/` directory (created automatically).

| Script | Description | Input |
|---|---|---|
| `build_das_real.py` | DAS reconstruction figure using real nylon-string RF data | `NylonString_Rcv_35.mat` |
| `build_das_reconstruction.py` | Diagrammatic DAS explanation figure | None (synthetic) |
| `build_iq_steps.py` | Step-by-step IQ demodulation figure | None (synthetic) |
| `build_matched_filter_v3.py` | Matched filter explanation figure | `rf_trace.mat` |
| `build_pipeline_figure.py` | Pressure signal-processing pipeline figure | `.mat`, `.csv`, or pickled data (configure in file) |
| `build_waveform.py` | Waveform figure from MATLAB `.fig` file | `waveform.fig` |
| `figure_directivity_portrait.py` | Array directivity polar patterns | None (synthetic) |
| `resolution_figure_v2.py` | Resolution analysis with –6 dB profiles | `resolution_export.mat` (default) or custom path |

```bash
# Run any script from within the python_graphs/ directory:
cd python_graphs
python build_das_reconstruction.py
python resolution_figure_v2.py                       # uses resolution_export.mat
python resolution_figure_v2.py path/to/export.mat    # custom .mat file
```

Scripts that require measured data (`.mat` or `.fig` files) expect those files in the working directory, or the path can be edited at the top of the script.

---

### pressure_characteristics/
MATLAB scripts for computing and visualizing acoustic pressure from water-tank hydrophone measurements.

- **`Compute_Pressure_Holistic_v4.m`** — converts a raw voltage waveform (`.txt`) to calibrated pressure using a hydrophone sensitivity curve (`.csv`); outputs spectrum metrics to the console and four figures.
- **`plot_beam_profile.m`** / **`plot_pressure_map.m`** — load a 2D scan file (`Scan_2D_UMSmap.txt`) and produce heatmap visualizations in linear and dB scale.
- **`plot_3D_beam.m`** — reconstructs a 3D beam volume from XZ and YZ 2D scans and renders isosurfaces at –3 dB and –6 dB.

Input file paths are hardcoded at the top of each script and must be updated to point to the local data files before running.

---

### reconstruction_scripts/
MATLAB beamforming pipelines that take raw multichannel RF data and produce B-mode images.

- **`PMUT_coherence_v4.m`** — full PMUT pipeline: channel reordering → matched filtering → IQ demodulation → DAS + coherence-factor beamforming. Loads a `.mat` file with raw RF data and displays B-mode images.
- **`L7_4_Recon.m`** — same pipeline adapted for the clinical L7-4 linear array probe. Also exports a `resolution_export.mat` file consumed by `python_graphs/resolution_figure_v2.py`.

Update the input file path at the top of each script to point to the appropriate acquisition `.mat` file.

---

### volumetric_simulation/
All simulation scripts are self-contained — no measured data files required.

- **`ChipMUT_3Dsim.m`** — Huygens–Fresnel 3D simulation of the current PMUT chip. Renders a 3D scatter plot of the intensity volume and 2D axial slices (XZ / YZ) with on-axis depth profiles. Array geometry and acoustic parameters are configured via variables at the top of the script.

- **`pmut_geometry_sweep.m`** — parametric sweep across alternative array geometries. Computes focal metrics (peak depth, FWHM, focal gain, depth of field) for each configuration and prints a comparison table. Optionally loads `alpha_cal.mat` for absolute pressure prediction (set `USE_MEASURED_CALIBRATION = true`). Requires the helper functions `pmut_metrics.m` and `pressure_budget.m` to be on the MATLAB path.

- **`pmut_geometry_sweep_optimized.m`** — extends the geometry sweep by steering the electronic focus across a range of depths (10–40 mm) for each configuration. Produces capability curves (intensity at intended focus vs. depth), snapshot comparisons at a fixed steered depth, and a summary table of best/worst performance per geometry.

#### auxilliary_sim_files/

- **`beam_hourglass_visualizer.m`** — 2D elevational beam analysis for multiple cell-layout configurations (`equal`, `parabolic`, `grouped`, `chirp`, `hybrid`, `multifocal_hybrid`). Plots the hourglass beam map, on-axis intensity, elevational FWHM vs. depth, and an intensity–resolution tradeoff scatter. Useful for comparing Fresnel-zone layout strategies.

- **`pmut_beam_3d_holistic.m`** — 3D volumetric field renderer for `equal`, `parabolic`, and `grouped` layouts. Produces orthogonal slice planes through the focal point and nested isosurfaces at –3, –6, and –12 dB to visualise the full 3D beam shape.

---

### miscellaneous/
- **`PSF_figure_simple.m`** — loads three water-tank phantom acquisitions and plots side-by-side B-mode heatmaps for PSF comparison. Update the three hardcoded file paths before running.


## Contact

For any questions regarding this work, feel free to reach out:

- Alexandru Ciocau — [alexandru-mihai.ciocau@student.kuleuven.be]
- Maarten Luypaert — [maarten.luypaert@student.kuleuven.be]

---

> This code is made publicly available for the purpose of academic evaluation.  
> Please contact the authors before reusing any part of this work.