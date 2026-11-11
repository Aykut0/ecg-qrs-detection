# ECG QRS Detection (Pan-Tompkins)

A training project: the Pan-Tompkins QRS detector written from scratch and scored against the beat labels cardiologists made on the same recording.

![Detected R-peaks and beat-to-beat heart rate](figures/ecg_rpeaks.png)

Pan-Tompkins (1985) is the classic algorithm for finding heartbeats in an ECG. I wanted to understand why each stage is there, and to check my implementation against ground truth instead of trusting a plot that looks right. The test signal is a 5-minute excerpt of MIT-BIH record 208 at 360 Hz.

The filtering follows the paper: a 5–15 Hz band-pass to keep QRS energy and reject baseline wander and T waves, a derivative, squaring, and a 150 ms moving-window integration that turns each QRS into a single hump.

![Pan-Tompkins stages on the real ECG](figures/ecg_pipeline.png)

The decision stage is where most from-scratch versions cut corners, so the script has two. `detect_static` uses one fixed threshold for the whole record. `detect_adaptive` follows the paper: running signal and noise peak estimates, a search-back when an RR interval runs long, and T-wave rejection based on slope. Scored against 509 annotated beats with the standard 150 ms tolerance:

| | detected | correct | false positives | missed | sensitivity | positive predictivity |
|---|---:|---:|---:|---:|---:|---:|
| static threshold | 522 | 499 | 23 | 10 | 98.04% | 95.59% |
| adaptive Pan-Tompkins | 499 | 496 | 3 | 13 | 97.45% | 99.40% |

The adaptive version cuts false positives from 23 to 3 and gives up a little sensitivity in noisy stretches.

![Validation against the reference annotations](figures/ecg_validation.png)

MIT-BIH annotators also mark signal quality (`~` where quality changes, `|` for artifacts that look like a QRS). 14 of the adaptive detector's 16 errors fall inside stretches they flagged: the single false positive at 19.9 s lands on a `|` artifact, and the eight consecutive misses between 209 and 214 s sit in a marked segment with about 5 mV of motion artifact. Two isolated premature ventricular contractions are real misses. Record 208 is a hard case, with 29% abnormal beats (93 PVCs and 56 fusion beats out of 509); the published Pan-Tompkins figures (99.3% sensitivity, 99.5% PPV) are averages over all 48 records.

## Run

```bash
pip install -r requirements.txt
python ecg_qrs_pan_tompkins.py
```

No network access is needed. The script reads the files in `data/`, writes the figures and `metrics.txt`, and should reproduce the table above exactly (499 detections, 3 false positives, PPV 99.40%).

`tools/fetch_annotations.py` shows how the reference files were made. It needs network access and `pip install wfdb`.

## Data

Samples 422820–530820 of MIT-BIH record 208 (channel MLII), the same excerpt SciPy ships as `scipy.datasets.electrocardiogram`. The offset is found by matching samples against the full record, with zero residual, which is what lines the annotations up with the detector output. Citation and license details are in `data/README.md`.

## Reference

Pan, J. and Tompkins, W.J. "A Real-Time QRS Detection Algorithm." IEEE Transactions on Biomedical Engineering, BME-32(3), 230-236, 1985.

## License

Code is MIT (see `LICENSE`). The MIT-BIH data has its own citation requirements.

---

Yaşar Aykut Altuntop · Koç University · [LinkedIn](https://linkedin.com/in/aykutaltuntop) · [GitHub](https://github.com/Aykut0)
