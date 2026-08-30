"""
Regenerate the ground-truth annotation file from PhysioNet.

The repository already ships `data/annotations_record208_excerpt.csv`, so you do NOT
need to run this — it exists to document exactly how that file was produced, and to let
anyone reproduce it from the original source.

What it does:

  1. downloads MIT-BIH record 208 (full 30 min) and its reference annotations,
  2. locates the 5-minute excerpt used by this project inside the full record by direct
     sample matching against `data/ecg_mitbih_record208.csv`,
  3. keeps only true *beat* annotations (AAMI beat labels), discarding rhythm markers,
     noise markers and comments,
  4. writes the beat positions, re-indexed to the excerpt, to the CSV,
  5. writes the annotators' signal-quality markers ('~' quality change, '|' isolated
     QRS-like artifact) to a second CSV, used by the error analysis.

Requires network access and the `wfdb` package (not a dependency of the main script):

    pip install wfdb
    python tools/fetch_annotations.py
"""
import os
import numpy as np

# AAMI beat labels. Everything else in an MIT-BIH annotation file is a non-beat marker:
# '+' rhythm change, '~' signal-quality change, '|' isolated QRS-like artifact, '"' comment.
BEAT_SYMBOLS = set("NLRejAaJSVEFP/fQ")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")


def main():
    import wfdb

    excerpt = np.loadtxt(os.path.join(DATA, "ecg_mitbih_record208.csv"))
    record = wfdb.rdrecord("208", pn_dir="mitdb")
    full = record.p_signal[:, 0]           # channel 0 = MLII
    fs = record.fs

    offset = locate_excerpt(full, excerpt)
    print(f"excerpt located at sample {offset} ({offset / fs:.1f} s into record 208)")

    ann = wfdb.rdann("208", "atr", pn_dir="mitdb")
    inside = (ann.sample >= offset) & (ann.sample < offset + len(excerpt))
    samples = ann.sample[inside] - offset
    symbols = np.array(ann.symbol)[inside]

    is_beat = np.array([s in BEAT_SYMBOLS for s in symbols])
    beats, labels = samples[is_beat], symbols[is_beat]
    print(f"{inside.sum()} annotations in excerpt, of which {len(beats)} are beats")

    out = os.path.join(DATA, "annotations_record208_excerpt.csv")
    with open(out, "w") as fh:
        fh.write("# MIT-BIH record 208 reference beat annotations\n")
        fh.write(f"# excerpt = samples {offset}..{offset + len(excerpt)} of the full record, fs = {fs:g} Hz\n")
        fh.write("# sample index is relative to the start of the excerpt\n")
        fh.write("sample,symbol\n")
        for s, y in zip(beats, labels):
            fh.write(f"{s},{y}\n")
    print("wrote", out)

    # Signal-quality markers, used to explain where the detector's errors fall.
    # '~' marks a change in signal quality, '|' an isolated QRS-like artifact.
    is_quality = np.isin(symbols, ["~", "|"])
    q_samples, q_labels = samples[is_quality], symbols[is_quality]
    print(f"{len(q_samples)} signal-quality markers")

    out_q = os.path.join(DATA, "quality_markers_record208_excerpt.csv")
    with open(out_q, "w") as fh:
        fh.write("# MIT-BIH record 208 annotator signal-quality markers\n")
        fh.write("# '~' = signal-quality change, '|' = isolated QRS-like artifact\n")
        fh.write("# sample index is relative to the start of the excerpt\n")
        fh.write("sample,symbol\n")
        for s, y in zip(q_samples, q_labels):
            fh.write(f"{s},{y}\n")
    print("wrote", out_q)


def locate_excerpt(full, excerpt, probe_len=2000, block=20000):
    """Find where `excerpt` begins inside `full` by minimising absolute difference."""
    from numpy.lib.stride_tricks import sliding_window_view

    probe = excerpt[:probe_len]
    windows = sliding_window_view(full, probe_len)
    best, best_err = None, np.inf
    for start in range(0, windows.shape[0], block):
        chunk = windows[start:start + block]
        err = np.abs(chunk - probe).sum(axis=1)
        i = int(np.argmin(err))
        if err[i] < best_err:
            best_err, best = err[i], start + i
    if best_err > 1e-6:
        raise SystemExit(f"excerpt not found in the full record (residual {best_err:.4g})")
    return best


if __name__ == "__main__":
    main()
