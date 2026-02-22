# HEASoft-Code
Modular pipeline created to run NICER data through HEASoft.

HEASoft + NICER tools + nibackgen3C50 must be installed and on your PATH for the program to run.

## Config parameters

- **star_name**: Name of the target you wish to analyse.
- **orbital_period**: Orbital period of object (in days).
- **reference_epoch**: Phase-zero epoch (**T0**) in MJD (e.g. Tasc / inferior conjunction from literature).
- **working_dir**: Directory for the program to operate in. Use `"working"` for the directory containing the `.py` file, or provide an absolute path.
- **run_nicer_l3**: `true` to run `nicerl3` instead of `nicerl2 + XSELECT`.
- **ObsID_to_exclude**: ObsIDs to skip. Example: `["0060010101", "0060010102"]`. Use `"current"` to exclude ObsID already downloaded.
- **ObsID_specifically_chosen**: ObsIDs to use (downloads only these). Example: `["0060010101", "0060010102"]`.
- **number_of_datasets**: How many datasets to download. If fewer exist, the program uses what’s available.
- **lower_energy_limit**: Lower energy cut in PI (1 PI = 0.01 keV).
- **upper_energy_limit**: Upper energy cut in PI (1 PI = 0.01 keV).
- **bins**: Number of bins used for binned plots.
- **epoch_range**: Epoch width (days) used for stacking.
- **auto_resolve**: `true` to backtrack and run missing prerequisites automatically.
- **reprocess**: `true` to regenerate files even if they already exist.
- **flare_filtering**: `true` to apply simple 3σ clipping on count rate to flag/remove outliers (toy model).

## run_steps (set each to `true` or `false`)

- **download**: Downloads datasets based on target selection rules.
- **process_l2**: Runs `nicerl2` processing.
- **barycorr**: Runs barycentric correction on L2 products.
- **barycurve**: Extracts light curve from barycentered event file.
- **remove_background**: Estimates and subtracts NICER background using `nibackgen3C50`.
- **statistics**: Computes basic statistics used for uncertainties in plots.
- **time_plot**: Time-domain plot for each ObsID.
- **unbinned_plot**: Phase-folded (unbinned) plot for each ObsID.
- **binned_plot**: Phase-folded (binned) plot for each ObsID.
- **stacked_plot**: Stacks binned plots into epoch groups defined by `epoch_range`.

---

This program was developed by Caden Phillips, Third Year MPhys Astrophysics Student, University of Liverpool.
