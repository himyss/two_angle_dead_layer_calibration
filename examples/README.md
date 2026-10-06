# Input formats

## Prepared peak files

`macro/twoAngleDeadLayer.C` can read simple text peak files. Each non-comment line
contains one strip and the peak positions in ADC channels:

```text
# strip low middle high
0 812.3 1024.8 1420.5
1 809.7 1021.4 1417.9
```

The normal-angle and angled files must contain the same strips and the same number
of peaks.

## Channel histogram files

`macro/twoAngleDeadLayerFromChannelHists.C` and
`macro/twoAngleFullCalibrationReport.C` read the `channelsHist_low`,
`channelsHist_middle`, and `channelsHist_high` ROOT files produced by the old
one-spectrum workflow.

Expected layout:

```text
results/
  <normal_file_base>/<detector>/
    channelsHist_low_<normal_file_base>_<detector>.root
    channelsHist_middle_<normal_file_base>_<detector>.root
    channelsHist_high_<normal_file_base>_<detector>.root
  <angled_file_base>/<detector>/
    channelsHist_low_<angled_file_base>_<detector>.root
    channelsHist_middle_<angled_file_base>_<detector>.root
    channelsHist_high_<angled_file_base>_<detector>.root
```
