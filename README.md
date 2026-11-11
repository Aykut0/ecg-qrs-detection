# QRS Detection on a Real ECG (Pan-Tompkins)

> **A training project.** A from-scratch heartbeat detector, checked against the beat labels cardiologists made on the same recording.

![Detected R-peaks and beat-to-beat heart rate](figures/ecg_rpeaks.png)

This is a from-scratch implementation of the Pan-Tompkins QRS detector, the classic 1985
algorithm for finding heartbeats in an ECG. It's a learning exercise, not original
research. The algorithm is well known and plenty of people have implemented it before. What
I wanted out of this was to actually understand why it's built the way it is, and then hold
my own implementation to the same standard a real detector would be held to: check it
against ground truth instead of assuming it works because the plot looks right.

I ran it on a 5-minute excerpt of MIT-BIH record 208, at 360 Hz, and scored the output
against the beat-by-beat annotations that cardiologists made on that same recording. That
scoring is really the point of the project. A detector nobody has validated against a real
reference is just a demo.

## Two versions, because the difference is the interesting part

The filtering side of the algorithm is the same in both versions: band-pass at 5-15 Hz to
keep QRS energy while rejecting baseline wander and T waves, a derivative to emphasize
steep slopes, squaring so everything is positive and the steep parts get amplified more,
then a moving-window integration (150 ms, about one QRS width) that turns the spiky
derivative into one hump per beat. Every implementation gets this part right because it's
mostly signal-processing plumbing.

Where implementations usually cut corners is the decision stage, deciding which hump in
that integrated signal is actually a heartbeat. Most from-scratch versions just pick one
threshold for the whole recording and call it done. That's what I built first (it's in the
code as `detect_static`), and it's not what Pan and Tompkins actually specified. Their real
method keeps running estimates of a signal peak and a noise peak that update continuously,
uses the recent RR interval to decide whether a beat arrived too late and needs a
search-back to recover something that got missed, and rejects T waves by comparing how
steep a new peak is against the one before it. I implemented that too (`detect_adaptive`)
and kept both versions in the script so the difference could be measured instead of just
asserted.

It matters. Scored against 509 cardiologist-annotated beats with the standard 150 ms
matching tolerance:

| | detected | correct | false positives | missed | sensitivity | positive predictivity |
|---|---:|---:|---:|---:|---:|---:|
| static threshold | 522 | 499 | 23 | 10 | 98.04% | 95.59% |
| adaptive Pan-Tompkins | 499 | 496 | 3 | 13 | 97.45% | 99.40% |

Doing the decision stage properly cuts the false positives from 23 down to 3. Sensitivity
drops a little, because the adaptive version is more cautious in noisy stretches, and that
trade-off is right there in the table instead of hidden behind one summary number.

![Validation against the reference annotations](figures/ecg_validation.png)

## Where the remaining errors come from

I didn't stop at the count. MIT-BIH annotators mark more than just beats, they also flag
signal quality: a `~` where the quality changes, and a `|` for something that looks like a
QRS but isn't one. I checked every error the adaptive detector made against those markers.

Of its 16 total errors, 14 sit inside a stretch the annotators themselves flagged. The one
false positive fires at 19.9 seconds, right on a sample marked `|`, an artifact that looks
enough like a real beat that the algorithm can be forgiven for it. The eight beats missed in
a row between 209 and 214 seconds happen inside a segment bracketed by quality markers,
where the baseline is swinging through about 5 mV of motion artifact. Only two errors, both
isolated premature ventricular contractions, have no such excuse and are genuine misses.

I picked record 208 for a reason: 29% of its beats are abnormal (93 premature ventricular
contractions and 56 fusion beats among 509 total), so it's a meaningfully harder test than
an average MIT-BIH record. Published Pan-Tompkins numbers, sensitivity 99.3%, PPV 99.5%, are
averaged over all 48 records in the database, so this result should read against that
context rather than against those numbers directly.

## The data

The signal is samples 422820 to 530820 of MIT-BIH record 208, the same 5-minute excerpt
that ships inside SciPy as `scipy.datasets.electrocardiogram`. I didn't just assume that
offset; `tools/fetch_annotations.py` locates the excerpt inside the full 30-minute record by
matching samples directly, and finds it at 422820 with zero residual difference. That's what
makes the reference annotations line up with the detector's output sample-for-sample.

Full citation and licence details for the MIT-BIH database are in `data/README.md`. Since
the annotation files are small text extracts bundled here for convenience, please check the
current PhysioNet licence terms before treating this as a public reference, as terms can
change between dataset versions.

## Running it

```bash
pip install -r requirements.txt
python ecg_qrs_pan_tompkins.py
```

That's it, no network access needed. It reads the three files already in `data/`, writes
the figures into `figures/`, and writes `metrics.txt`. If you run it yourself the numbers
should come out exactly as in the table above (499 detections, 3 false positives, PPV
99.40%). If they don't match, something about the environment is different and that's worth
tracking down.

The ground-truth files in `data/` were generated once from PhysioNet and are already
included, so you don't need to regenerate them. If you want to see exactly how they were
made, or reproduce them yourself, `tools/fetch_annotations.py` does that: it needs network
access and the `wfdb` package, which isn't a dependency of the main script.

```bash
pip install wfdb
python tools/fetch_annotations.py
```

![The Pan-Tompkins stages on the real ECG](figures/ecg_pipeline.png)

## Why I built this

QRS detection sits in front of almost every ECG application: arrhythmia classification,
heart-rate variability, cardiac gating. I work on wearable ultrasound blood-pressure
monitoring at the Bio-Integrated Microdevices Laboratory at Koç University, where beat
detection triggers the acquisition window for each pulse. Getting that right on a noisy
recording, and actually knowing which beats it misses and why, is the gap between something
that works on a bench and something that works on a person. This project was my way of
closing that gap on the algorithm side before trusting it on hardware.

## Licence

Code is MIT, see `LICENSE`. The MIT-BIH data carries its own citation requirements, see
`data/README.md`.

## Reference

Pan, J. and Tompkins, W.J. "A Real-Time QRS Detection Algorithm." IEEE Transactions on
Biomedical Engineering, BME-32(3), 230-236, 1985.

---

**Author:** Yaşar Aykut Altuntop · MSc student, Biomedical Sciences and Engineering, Koç University  
[LinkedIn](https://linkedin.com/in/aykutaltuntop) · [GitHub](https://github.com/Aykut0)
