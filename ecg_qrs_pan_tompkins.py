"""
QRS detection on a real ECG — Pan & Tompkins (1985), built and measured
=======================================================================
Implements the Pan-Tompkins QRS detector from scratch and scores it against the
reference beat annotations made by cardiologists, using the standard evaluation
protocol (ANSI/AAMI EC57: a detection counts as correct if it falls within 150 ms
of an annotated beat).

Data
----
A 5-minute excerpt of **MIT-BIH record 208**, sampled at 360 Hz — samples
422820..530820 of the full record, the excerpt distributed with SciPy as
`scipy.datasets.electrocardiogram`. Record 208 is one of the harder records in the
database: of the 509 beats in this excerpt, 93 are premature ventricular
contractions and 56 are fusion beats, so 29 % of the beats are abnormal. It also
carries real baseline wander and stretches flagged by the annotators as noisy.

Two thresholding strategies
---------------------------
The filter chain (band-pass, derivative, squaring, moving-window integration) is
common to both. What differs is how a peak in the integrated signal is judged:

  `static`    a single threshold computed once from the whole record. Simple, and
              what most from-scratch implementations actually do.

  `adaptive`  what the paper actually specifies, and the reason the algorithm works:
              running signal and noise peak estimates (SPKI / NPKI) that update on
              every peak, a threshold derived from them, RR-interval limits with
              search-back to recover beats missed during a dropout, and T-wave
              discrimination by slope comparison.

Both are implemented and both are scored, because the difference between them is
the interesting result — see the table this script prints.

Run
---
    python ecg_qrs_pan_tompkins.py

Everything runs offline from the two files in `data/`. Writes three figures into
`figures/`, plus `metrics.txt`.

Author: Yasar Aykut Altuntop
"""
import os
import numpy as np
from scipy.signal import butter, filtfilt, find_peaks
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FS = 360.0
MATCH_TOLERANCE_S = 0.150          # ANSI/AAMI EC57

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
DATA = os.path.join(HERE, "data")


# ----------------------------------------------------------------- data

def load_ecg():
    return np.loadtxt(os.path.join(DATA, "ecg_mitbih_record208.csv"))


def _read_marker_csv(filename):
    samples, symbols = [], []
    with open(os.path.join(DATA, filename)) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("sample,"):
                continue
            s, y = line.split(",")
            samples.append(int(s))
            symbols.append(y)
    return np.array(samples), np.array(symbols)


def load_annotations():
    """Reference beat positions, in samples relative to the start of the excerpt."""
    return _read_marker_csv("annotations_record208_excerpt.csv")


def load_quality_markers():
    """Annotator signal-quality markers: '~' quality change, '|' QRS-like artifact."""
    return _read_marker_csv("quality_markers_record208_excerpt.csv")


# ------------------------------------------------------- the filter chain

def preprocess(x, fs):
    """Band-pass -> derivative -> square -> moving-window integration."""
    b, a = butter(3, [5.0 / (fs / 2), 15.0 / (fs / 2)], btype="band")
    bandpassed = filtfilt(b, a, x)
    derivative = np.gradient(bandpassed) * fs
    squared = derivative ** 2
    width = max(1, int(0.150 * fs))                     # 150 ms, ~one QRS width
    integrated = np.convolve(squared, np.ones(width) / width, mode="same")
    return bandpassed, derivative, squared, integrated


def refine_to_r_peak(bandpassed, candidates, fs):
    """Move each detection onto the largest local deflection, +/- 50 ms."""
    half = int(0.05 * fs)
    out = []
    for c in candidates:
        a0, a1 = max(0, c - half), min(len(bandpassed), c + half)
        out.append(a0 + int(np.argmax(np.abs(bandpassed[a0:a1]))))
    return np.unique(out)


# ------------------------------------------------------------ detectors

def detect_static(bandpassed, integrated, fs):
    """One threshold for the whole record."""
    threshold = 0.15 * np.mean(integrated) + 0.20 * np.median(integrated[integrated > 0])
    candidates, _ = find_peaks(integrated, distance=int(0.25 * fs), height=threshold)
    return refine_to_r_peak(bandpassed, candidates, fs), threshold


def detect_adaptive(bandpassed, integrated, fs):
    """Pan-Tompkins as specified: adaptive thresholds, RR limits, search-back, T-wave rejection."""
    peaks, _ = find_peaks(integrated, distance=int(0.200 * fs))     # 200 ms refractory period
    if len(peaks) == 0:
        return np.array([], dtype=int), np.array([])

    # initialise from the first two seconds
    warmup = integrated[:int(2 * fs)]
    spki = float(np.max(warmup))            # running estimate of a signal peak
    npki = float(np.mean(warmup))           # running estimate of a noise peak
    thr1 = npki + 0.25 * (spki - npki)
    thr2 = 0.5 * thr1

    slope = np.abs(np.gradient(bandpassed))
    half = int(0.02 * fs)

    def local_slope(i):
        return slope[max(0, i - half):i + half].max()

    qrs, rr_recent, rr_normal, rr_average = [], [], [], None
    threshold_trace = np.empty(len(peaks))

    for n, p in enumerate(peaks):
        peak_value = integrated[p]
        threshold_trace[n] = thr1
        accepted = peak_value > thr1

        # a peak within 360 ms of the last QRS with less than half its slope is a T wave
        if accepted and qrs and (p - qrs[-1]) < 0.360 * fs:
            if local_slope(p) < 0.5 * local_slope(qrs[-1]):
                accepted = False
                npki = 0.125 * peak_value + 0.875 * npki

        if not accepted:
            if peak_value <= thr1:
                npki = 0.125 * peak_value + 0.875 * npki
        else:
            # search-back: if this beat arrives too late, look for one we rejected
            if qrs and rr_average and (p - qrs[-1]) > 1.66 * rr_average:
                lo, hi = qrs[-1] + int(0.200 * fs), p
                missed = [q for q in peaks if lo < q < hi and integrated[q] > thr2]
                if missed:
                    recovered = max(missed, key=lambda q: integrated[q])
                    qrs.append(recovered)
                    spki = 0.25 * integrated[recovered] + 0.75 * spki    # lower weight on search-back
                    if len(qrs) > 1:
                        rr_recent = (rr_recent + [qrs[-1] - qrs[-2]])[-8:]

            qrs.append(p)
            spki = 0.125 * peak_value + 0.875 * spki
            if len(qrs) > 1:
                rr = qrs[-1] - qrs[-2]
                rr_recent = (rr_recent + [rr])[-8:]
                if rr_average is None or 0.92 * rr_average < rr < 1.16 * rr_average:
                    rr_normal = (rr_normal + [rr])[-8:]
                rr_average = float(np.mean(rr_normal)) if rr_normal else float(np.mean(rr_recent))

        thr1 = npki + 0.25 * (spki - npki)
        thr2 = 0.5 * thr1

    return refine_to_r_peak(bandpassed, np.unique(qrs), fs), threshold_trace


# ----------------------------------------------------------- evaluation

def score(detections, reference, fs, tolerance_s=MATCH_TOLERANCE_S):
    """Match detections to reference beats; return counts and per-detection labels."""
    tolerance = tolerance_s * fs
    matched = np.zeros(len(reference), dtype=bool)
    labels = []
    for d in detections:
        k = int(np.argmin(np.abs(reference - d)))
        if abs(reference[k] - d) <= tolerance and not matched[k]:
            matched[k] = True
            labels.append("TP")
        else:
            labels.append("FP")
    tp = int(matched.sum())
    fp = len(detections) - tp
    fn = int((~matched).sum())
    return dict(
        tp=tp, fp=fp, fn=fn,
        sensitivity=100.0 * tp / (tp + fn),
        ppv=100.0 * tp / (tp + fp),
        labels=np.array(labels),
        missed=reference[~matched],
    )


ERROR_CONTEXT_S = 3.0


def error_context(result, quality, fs, window_s=ERROR_CONTEXT_S):
    """How many of this detector's errors fall in a stretch the annotators flagged?

    The annotators mark changes in signal quality ('~') and isolated QRS-like
    artifacts ('|'). An error near one of those is a consequence of the recording,
    not of the algorithm, so counting them separates the two.
    """
    false_pos = result["detections"][result["labels"] == "FP"]
    errors = np.concatenate([false_pos, result["missed"]])
    if len(errors) == 0 or len(quality) == 0:
        return dict(total=len(errors), flagged=0, unexplained=np.array([], dtype=int))
    distance = np.array([np.min(np.abs(quality - e)) for e in errors]) / fs
    near = distance <= window_s
    return dict(total=len(errors), flagged=int(near.sum()), unexplained=errors[~near])


# -------------------------------------------------------------- figures

def figure_pipeline(t, ecg, bandpassed, squared, integrated, threshold, rpeaks, window):
    w = (t >= window[0]) & (t <= window[1])
    rw = rpeaks[(t[rpeaks] >= window[0]) & (t[rpeaks] <= window[1])]
    fig, ax = plt.subplots(4, 1, figsize=(10, 7), sharex=True)
    ax[0].plot(t[w], ecg[w], color="#374151", lw=1)
    ax[0].set_ylabel("raw ECG\n[mV]")
    ax[0].set_title("Pan-Tompkins stages on real ECG (MIT-BIH record 208)")
    ax[1].plot(t[w], bandpassed[w], color="#2563eb", lw=1)
    ax[1].set_ylabel("band-pass\n5-15 Hz")
    ax[2].plot(t[w], squared[w], color="#7c3aed", lw=1)
    ax[2].set_ylabel("derivative²")
    ax[3].plot(t[w], integrated[w], color="#059669", lw=1)
    ax[3].axhline(threshold, color="#dc2626", ls="--", lw=1, label="threshold")
    ax[3].plot(t[rw], integrated[rw], "o", color="#dc2626", ms=5)
    ax[3].set_ylabel("integrated\n(150 ms)")
    ax[3].set_xlabel("time [s]")
    ax[3].legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "ecg_pipeline.png"), dpi=140)
    plt.close(fig)


def figure_rpeaks(t, ecg, rpeaks, window):
    w = (t >= window[0]) & (t <= window[1])
    rw = rpeaks[(t[rpeaks] >= window[0]) & (t[rpeaks] <= window[1])]
    rr = np.diff(t[rpeaks])
    keep = (rr > 0.3) & (rr < 2.0)
    hr = 60.0 / rr[keep]
    hr_t = t[rpeaks][1:][keep]

    fig, ax = plt.subplots(2, 1, figsize=(10, 5.5))
    ax[0].plot(t[w], ecg[w], color="#374151", lw=1)
    ax[0].plot(t[rw], ecg[rw], "v", color="#dc2626", ms=8, label="detected R-peak")
    ax[0].set_title(f"Detected R-peaks on the real ECG ({window[1]-window[0]:.0f} s window)")
    ax[0].set_ylabel("ECG [mV]")
    ax[0].set_xlabel("time [s]")
    ax[0].legend(loc="upper right")
    ax[1].plot(hr_t, hr, ".", color="#0891b2", ms=3)
    ax[1].set_title(f"Beat-to-beat heart rate over 5 minutes (mean {np.mean(hr):.0f} bpm)")
    ax[1].set_xlabel("time [s]")
    ax[1].set_ylabel("HR [bpm]")
    ax[1].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "ecg_rpeaks.png"), dpi=140)
    plt.close(fig)


def _error_panel(ax, t, ecg, adaptive, reference, quality, span, title):
    lo, hi = span
    w = (t >= lo) & (t <= hi)

    for q in quality:
        if lo <= t[q] <= hi:
            ax.axvline(t[q], color="#b45309", ls=":", lw=1.2, zorder=0)

    ax.plot(t[w], ecg[w], color="#374151", lw=0.9, zorder=1)

    ref_w = reference[(t[reference] >= lo) & (t[reference] <= hi)]
    ax.plot(t[ref_w], ecg[ref_w], "o", mfc="none", mec="#059669", ms=11, mew=1.4, zorder=2)

    det, lab = adaptive["detections"], adaptive["labels"]
    tp = det[(lab == "TP") & (t[det] >= lo) & (t[det] <= hi)]
    fp = det[(lab == "FP") & (t[det] >= lo) & (t[det] <= hi)]
    ax.plot(t[tp], ecg[tp], "v", color="#2563eb", ms=7, zorder=3)
    ax.plot(t[fp], ecg[fp], "X", color="#dc2626", ms=10, zorder=4)

    miss = adaptive["missed"]
    miss = miss[(t[miss] >= lo) & (t[miss] <= hi)]
    ax.plot(t[miss], ecg[miss], "s", mfc="none", mec="#dc2626", ms=14, mew=1.8, zorder=4)

    ax.set_title(title, fontsize=10)
    ax.set_xlabel("time [s]")
    ax.set_xlim(lo, hi)


def figure_validation(t, ecg, results, reference, quality):
    """Two error cases in context, plus the head-to-head bars."""
    adaptive = results["adaptive"]
    fig = plt.figure(figsize=(11, 7.8))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1], hspace=0.62, wspace=0.16)

    ax_fp = fig.add_subplot(gs[0, 0])
    _error_panel(ax_fp, t, ecg, adaptive, reference, quality, (17.5, 22.5),
                 "False positive at 19.9 s — on an artifact the annotators marked '|'")
    ax_fp.set_ylabel("ECG [mV]")

    ax_fn = fig.add_subplot(gs[0, 1])
    _error_panel(ax_fn, t, ecg, adaptive, reference, quality, (207.5, 215.5),
                 "Eight missed beats, 209-214 s — inside a stretch marked noisy")

    handles = [
        plt.Line2D([], [], ls="none", marker="o", mfc="none", mec="#059669", ms=10, mew=1.4,
                   label="annotated beat"),
        plt.Line2D([], [], ls="none", marker="v", color="#2563eb", ms=7, label="correct (TP)"),
        plt.Line2D([], [], ls="none", marker="X", color="#dc2626", ms=9, label="false positive"),
        plt.Line2D([], [], ls="none", marker="s", mfc="none", mec="#dc2626", ms=12, mew=1.8,
                   label="missed beat"),
        plt.Line2D([], [], ls=":", color="#b45309", lw=1.2, label="annotator quality marker"),
    ]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.495),
               ncol=5, fontsize=8.5, frameon=False)

    ax = fig.add_subplot(gs[1, :])
    names = ["static threshold", "adaptive Pan-Tompkins"]
    se = [results["static"]["sensitivity"], results["adaptive"]["sensitivity"]]
    pv = [results["static"]["ppv"], results["adaptive"]["ppv"]]
    fpc = [results["static"]["fp"], results["adaptive"]["fp"]]
    y = np.arange(len(names))
    ax.barh(y - 0.19, se, height=0.36, color="#0891b2")
    ax.barh(y + 0.19, pv, height=0.36, color="#dc2626")
    for i in range(len(names)):
        ax.text(90.15, y[i] - 0.19, "sensitivity", va="center", ha="left",
                color="white", fontsize=8.5, fontweight="bold")
        ax.text(90.15, y[i] + 0.19, "positive predictivity", va="center", ha="left",
                color="white", fontsize=8.5, fontweight="bold")
        ax.text(se[i] + 0.12, y[i] - 0.19, f"{se[i]:.2f} %", va="center", fontsize=9)
        ax.text(pv[i] + 0.12, y[i] + 0.19, f"{pv[i]:.2f} %   ({fpc[i]} false positives)",
                va="center", fontsize=9)
    ax.set_yticks(y, names)
    ax.set_xlim(90, 101.5)
    ax.set_xlabel("percent")
    ax.set_title("Scored against 509 cardiologist-annotated beats, 150 ms matching tolerance "
                 "(ANSI/AAMI EC57)", fontsize=10)
    ax.grid(axis="x", alpha=0.3)
    ax.set_axisbelow(True)
    fig.savefig(os.path.join(FIG, "ecg_validation.png"), dpi=140, bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------------------ main

def main():
    os.makedirs(FIG, exist_ok=True)
    ecg = load_ecg()
    t = np.arange(len(ecg)) / FS
    reference, symbols = load_annotations()
    quality, _ = load_quality_markers()

    bandpassed, _, squared, integrated = preprocess(ecg, FS)

    static_det, static_thr = detect_static(bandpassed, integrated, FS)
    adaptive_det, _ = detect_adaptive(bandpassed, integrated, FS)

    results = {}
    for name, det in [("static", static_det), ("adaptive", adaptive_det)]:
        r = score(det, reference, FS)
        r["detections"] = det
        r["context"] = error_context(r, quality, FS)
        results[name] = r

    rr = np.diff(t[adaptive_det])
    hr = 60.0 / rr[(rr > 0.3) & (rr < 2.0)]

    abnormal = int(np.sum(~np.isin(symbols, list("NLRej"))))
    lines = [
        "QRS detection on MIT-BIH record 208, 5-minute excerpt (real data)",
        f"record            : {len(ecg)/FS:.0f} s @ {FS:.0f} Hz, samples 422820-530820 of record 208",
        f"reference beats   : {len(reference)}  ({abnormal} abnormal - PVC, fusion or unclassifiable)",
        f"matching tolerance: {MATCH_TOLERANCE_S*1000:.0f} ms (ANSI/AAMI EC57)",
        "",
        f"{'detector':24s}{'det':>6}{'TP':>6}{'FP':>5}{'FN':>5}{'Se %':>9}{'PPV %':>9}",
    ]
    for name in ("static", "adaptive"):
        r = results[name]
        label = "static threshold" if name == "static" else "adaptive Pan-Tompkins"
        lines.append(f"{label:24s}{len(r['detections']):6d}{r['tp']:6d}{r['fp']:5d}"
                     f"{r['fn']:5d}{r['sensitivity']:9.2f}{r['ppv']:9.2f}")
    lines += [
        "",
        f"heart rate (adaptive) : mean {np.mean(hr):.1f} bpm, "
        f"{np.percentile(hr,5):.0f}-{np.percentile(hr,95):.0f} bpm (5-95th pct)",
        "",
        "Published Pan-Tompkins performance is Se 99.3 %, PPV 99.5 % averaged over all 48",
        "MIT-BIH records. Record 208 is materially harder than that average.",
        "",
        f"Where the errors fall (within {ERROR_CONTEXT_S:.0f} s of an annotator '~' or '|' marker)",
    ]
    for name in ("static", "adaptive"):
        c = results[name]["context"]
        label = "static threshold" if name == "static" else "adaptive Pan-Tompkins"
        pct = 100.0 * c["flagged"] / c["total"] if c["total"] else 0.0
        lines.append(f"  {label:24s}{c['flagged']:3d} of {c['total']:3d} errors "
                     f"({pct:.0f} %) sit in a stretch the annotators flagged")
    unexplained = results["adaptive"]["context"]["unexplained"]
    lines.append(f"  the {len(unexplained)} adaptive errors outside a flagged stretch are at "
                 f"t = {', '.join(f'{u/FS:.1f} s' for u in unexplained)}")

    with open(os.path.join(HERE, "metrics.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))

    figure_pipeline(t, ecg, bandpassed, squared, integrated, static_thr, adaptive_det, (10, 15))
    figure_rpeaks(t, ecg, adaptive_det, (10, 20))
    figure_validation(t, ecg, results, reference, quality)


if __name__ == "__main__":
    main()
