# Two-Angle Dead Layer Calibration

Small ROOT/C++ project for estimating a silicon detector dead layer from two
measurements of the same alpha source:

- normal incidence, for example `0 deg`;
- angled incidence, for example `45 deg`.

The code is a cleaned source-only extraction from `AculUti`. It keeps the
classes and macros needed to rebuild the libraries, calculate the dead layer,
write reports, draw spectra with found peak positions, and export calibration
coefficients in the same `a b` / `b a` text format used by the one-spectrum
workflow.

## Requirements

- ROOT with `root-config` and `rootcint` available in `PATH`;
- `g++`;
- `gfortran`;
- Linux shell environment.

## Build

```bash
make all
```

This builds:

- `libTELoss.so`;
- `libAculData.so`.

To clean generated files:

```bash
make clean
```

## Quick test

```bash
LD_LIBRARY_PATH=$PWD:$LD_LIBRARY_PATH root -l -b -q 'macro/testTwoAngleDeadLayer.C()'
```

Expected result: the synthetic test recovers the generated dead layer and prints
`testTwoAngleDeadLayer: ok ...`.

## Main workflow

The full workflow assumes that peak positions for both measurements are already
stored in `channelsHist_low`, `channelsHist_middle`, and `channelsHist_high`
files. These files can be produced by the old one-spectrum peak search.

Example:

```bash
LD_LIBRARY_PATH=$PWD:$LD_LIBRARY_PATH root -l -b -q \
'macro/twoAngleFullCalibrationReport.C(
  "/path/to/after_L1500_junction_Ra226_30cm_00deg_0001.root",
  "/path/to/after_L1500_junction_Ra226_30cm_45deg_0001.root",
  "ampDSSD_L_X",
  45.,
  "AnalysisxTree",
  "NeEvent",
  "/path/to/si_calibreation/results",
  "results/two_angle_full_ampDSSD_L_X",
  0,
  15
)'
```

The macro writes:

- `two_angle_dead_layer_report.txt`;
- `two_angle_dead_layer_report.csv`;
- `two_angle_dead_layer_graphs.root`;
- `report_<normal_file_base>_<detector>.txt`;
- `calibCoeff_<normal_file_base>_<detector>.txt`;
- `calibCoeff_<normal_file_base>_<detector>_rev.txt`;
- raw spectra images with found peak positions;
- calibrated spectra images with found peak positions.

## Other entry points

Use prepared text peak files:

```bash
LD_LIBRARY_PATH=$PWD:$LD_LIBRARY_PATH root -l -b -q \
'macro/twoAngleDeadLayer.C("normal.peaks", "angled.peaks", "ampDSSD_L_X", 45.)'
```

Use existing `channelsHist_*` files without drawing raw/calibrated spectra:

```bash
LD_LIBRARY_PATH=$PWD:$LD_LIBRARY_PATH root -l -b -q \
'macro/twoAngleDeadLayerFromChannelHists.C(
  "/path/to/results/<normal_file_base>/ampDSSD_L_X",
  "<normal_file_base>",
  "/path/to/results/<angled_file_base>/ampDSSD_L_X",
  "<angled_file_base>",
  "ampDSSD_L_X",
  45.,
  0,
  15
)'
```

## Algorithm

For each strip the method compares the same alpha peaks measured at two angles.
For a trial dead layer `d`:

- the normal spectrum uses the path length `d`;
- the angled spectrum uses `d / cos(theta)`;
- `TELoss` converts source alpha energies to detected energies after silicon
  energy loss;
- one temporary linear response is fitted to both angle datasets;
- the dead layer is chosen by minimizing the combined residual RMS;
- after that, final calibration coefficients are fitted for the normal spectrum.

For the current Ra-226 test pair the two-angle method gave a more stable dead
layer estimate across strips than the one-spectrum method. The practical
recommendation is to measure at `0 deg` and at `45 deg`, use the two-angle data
to determine the dead layer, and then, if needed, use that fixed dead layer in
the one-spectrum calibration workflow.

## Source Layout

- `AculData/` - calibration classes and two-angle analysis methods;
- `TELoss/` - silicon energy-loss model;
- `macro/twoAngleFullCalibrationReport.C` - full two-file analysis;
- `macro/twoAngleDeadLayerFromChannelHists.C` - dead layer from existing
  channel histograms;
- `macro/twoAngleDeadLayer.C` - dead layer from prepared text peak files;
- `macro/testTwoAngleDeadLayer.C` - synthetic regression test;
- `examples/` - input format notes.
