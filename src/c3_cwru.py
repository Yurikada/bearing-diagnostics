"""C3 ループ4: CWRU known-answer test。

宣言済み手順(エンベロープを見る前に帯域を確定):
1. 各欠陥ファイルの PSD / 正常97.mat の PSD 比を 1kHz 幅・500Hz 刻みの窓で評価し、
   中央値比が最大の窓を帯域として機械的に採用する(500-5500Hz)。
2. その帯域で envelope_spectrum (全長 ~10s -> Δf ~0.1Hz)。
3. 賭けシート(2026-07-27 宣言)と照合:
   105 IR   -> 162.2 Hz, 側帯 ±fr=29.95 Hz
   118 Ball -> 141.2 Hz (2xBSF), 側帯 ±FTF=11.93 Hz
   130 OR   -> 107.4 Hz, 側帯なし(極めて弱い)
   照合窓は ±2%。線の強さは ±0.5Hz 積分 / 床(20-300Hzの中央値ベース)。
"""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.io

sys.path.insert(0, str(Path(__file__).parent))
from lib.envelope import envelope_spectrum
from lib.spectral import welch_psd

ROOT = Path(__file__).resolve().parents[1]
FS = 12000.0
FR = 29.95
CASES = [
    ("105.mat", "IR 0.007in", 162.2, FR, "#0072B2"),
    ("118.mat", "Ball 0.007in", 141.2, 11.93, "#009E73"),
    ("130.mat", "OR 0.007in @6:00", 107.4, None, "#D55E00"),
]


def de_signal(fname):
    mat = scipy.io.loadmat(ROOT / "data/cwru" / fname)
    key = [k for k in mat if k.endswith("_DE_time")][0]
    return mat[key].ravel()


def pick_band(x_fault, x_normal):
    """正常/欠陥の Welch PSD 比(中央値)が最大の 1kHz 窓を返す。"""
    f, pf = welch_psd(x_fault, FS, nseg=2048)
    _, pn = welch_psd(x_normal, FS, nseg=2048)
    best, best_r = None, -1.0
    for lo in np.arange(500.0, 4600.0, 500.0):
        m = (f >= lo) & (f < lo + 1000.0)
        r = float(np.median(pf[m] / pn[m]))
        if r > best_r:
            best, best_r = lo, r
    return best, best + 1000.0, best_r


def line_strength(freqs, amp, f0, tol=0.02):
    """±tol 窓内のピークを探し (peak_freq, line, line/floor) を返す。"""
    df = freqs[1] - freqs[0]
    m = (freqs >= f0 * (1 - tol)) & (freqs <= f0 * (1 + tol))
    k = np.argmax(amp[m])
    fpk = freqs[m][k]
    kk = int(round(fpk / df))
    w = max(int(round(0.5 / df)), 1)
    line = float(amp[kk - w : kk + w + 1].sum())
    base = (freqs >= 20) & (freqs <= 300)
    floor = float(np.median(amp[base])) * (2 * w + 1)
    return fpk, line, line / floor


def main():
    xn = de_signal("97.mat")
    fig, axes = plt.subplots(3, 1, figsize=(10, 8.5), sharex=True)

    for ax, (fname, label, f_pred, f_side, color) in zip(axes, CASES):
        x = de_signal(fname)
        lo, hi, ratio = pick_band(x, xn[: x.size])
        freqs, amp = envelope_spectrum(x, FS, lo, hi)
        fpk, line, lf = line_strength(freqs, amp, f_pred)
        dev = (fpk - f_pred) / f_pred * 100

        print(f"{fname} {label:18s} band {lo:.0f}-{hi:.0f}Hz (ratio {ratio:6.1f})"
              f"  main: pred {f_pred:.1f} -> meas {fpk:.2f} Hz ({dev:+.2f}%)  line/floor {lf:8.1f}")
        if f_side:
            for s in (-1, +1):
                fs_, _, lfs = line_strength(freqs, amp, fpk + s * f_side)
                print(f"    sideband {'-' if s < 0 else '+'}{f_side:.2f}Hz:"
                      f" meas {fs_:.2f} Hz  line/floor {lfs:6.1f}")
        else:
            for s in (-1, +1):
                _, _, lfs = line_strength(freqs, amp, fpk + s * FR)
                print(f"    (check) sideband {'-' if s < 0 else '+'}fr: line/floor {lfs:6.1f}")

        m = freqs <= 250
        ax.plot(freqs[m], amp[m], color=color, lw=0.9)
        ax.axvline(f_pred, color="#555555", ls="--", lw=0.8)
        if f_side:
            for s in (-1, +1):
                ax.axvline(f_pred + s * f_side, color="#aaaaaa", ls=":", lw=0.8)
        ax.set_title(f"{fname}  {label}   band {lo:.0f}-{hi:.0f} Hz"
                     f"   (dashed = predicted {f_pred} Hz"
                     + (f", dotted = +/-{f_side} Hz sidebands)" if f_side else ", no sidebands predicted)"),
                     fontsize=9, loc="left")
        ax.set_ylabel("envelope amp [g]")
        ax.grid(True, color="#eeeeee", lw=0.5)
        ax.spines[["top", "right"]].set_visible(False)
    axes[-1].set_xlabel("frequency [Hz]")
    fig.tight_layout()
    out = ROOT / "docs/figs/C3-cwru-envelope.png"
    fig.savefig(out, dpi=150)
    print("figure:", out.name)


if __name__ == "__main__":
    sys.exit(main() or 0)
