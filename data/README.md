# Data provenance

Three files, all derived from **MIT-BIH Arrhythmia Database record 208**.

| file | what it is |
|---|---|
| `ecg_mitbih_record208.csv` | the ECG signal, 108 000 samples at 360 Hz (5 minutes) |
| `annotations_record208_excerpt.csv` | 509 reference beat positions and their labels |
| `quality_markers_record208_excerpt.csv` | 14 annotator signal-quality markers |

## Excerpt

The signal is samples **422820–530820** of record 208, channel MLII. That's the 5-minute
excerpt distributed with SciPy as `scipy.datasets.electrocardiogram`.

That offset is verified, not assumed. `tools/fetch_annotations.py` searches the full 30-minute
record for the excerpt by minimising absolute sample difference, and finds it at 422820 with
a residual of exactly zero. The annotation sample indices are then re-based to the start of
the excerpt, which is what makes them directly comparable to detector output.

## Annotation symbols

`annotations_record208_excerpt.csv` contains only true **beat** annotations, filtered to the
AAMI beat labels `N L R e j A a J S V E F P / f Q`:

| symbol | meaning | count |
|---|---|---:|
| `N` | normal sinus beat | 358 |
| `V` | premature ventricular contraction | 93 |
| `F` | fusion of ventricular and normal beat | 56 |
| `Q` | unclassifiable | 2 |

`quality_markers_record208_excerpt.csv` contains the non-beat markers used by the error
analysis:

| symbol | meaning |
|---|---|
| `~` | change in signal quality |
| `\|` | isolated QRS-like artifact |

Rhythm-change markers (`+`) and comments are discarded, since they are not detectable events
and counting them either way would distort the score.

## Citation

> Moody, G.B. and Mark, R.G. "The impact of the MIT-BIH Arrhythmia Database."
> *IEEE Engineering in Medicine and Biology Magazine*, 20(3), 45–50, 2001.

> Goldberger, A., Amaral, L., Glass, L., Hausdorff, J., Ivanov, P.C., Mark, R., Mietus,
> J.E., Moody, G.B., Peng, C.K., and Stanley, H.E. "PhysioBank, PhysioToolkit, and
> PhysioNet: Components of a New Research Resource for Complex Physiologic Signals."
> *Circulation* 101(23), e215–e220, 2000.

## Licence

The MIT-BIH Arrhythmia Database is published on PhysioNet under the
**Open Data Commons Attribution License v1.0 (ODC-BY 1.0)**, which permits redistribution
with attribution. These extracts are redistributed under that licence so the project runs
without a download step.

> **Confirm the licence on the PhysioNet page before making this repository public.**
> PhysioNet terms are per-dataset and can change between versions. If redistribution is not
> permitted, delete these files and have the README point at `tools/fetch_annotations.py`
> instead.

Source: <https://physionet.org/content/mitdb/>
