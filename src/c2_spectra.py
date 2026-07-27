"""C2 ループ4: B1 スペクトルの3時点比較 (Set2, welch_psd 自前実装で計算)。

固定: welch_psd(nseg=2048, Hann, 50%) = C2 合意。各時点は対象日を中心に
6スナップショットの PSD を平均(分散をさらに低減)。変えるもの: 時点のみ。

時点(C1 実測から宣言):
- day 1.0 : ベースライン
- day 3.8 : RMS 立ち上がり(3.70日)直後
- day 5.0 : 尖度立ち上がり(4.79日)後の最初のプラトー

賭け(描く前に宣言済み):
- Q1: day 3.8 で離散線が既に立っている(B)。潤滑仮説(広帯域のみ=A)と対立。
- Q2: day 5.0 の BPFO 高調波に fr=33.3Hz の側帯波は立たない(外輪固定・荷重一定)。
- fs 判定: BPFO 線の実測位置。20.48kHz 採用が正しければ 236.4Hz の -2%〜0%
  (すべり)近傍。もし真の fs が 20.00kHz なら見かけ +2.4% の 242.1Hz 側に出る。
"""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from lib.spectral import welch_psd

ROOT = Path(__file__).resolve().parents[1]
FS = 20480.0
NSEG = 2048
BPFO = 236.4
FR = 33.33
EPOCHS = [("day 1.0 (baseline)", 1.0, "#999999"),
          ("day 3.8 (RMS onset)", 3.8, "#0072B2"),
          ("day 5.0 (kurt onset)", 5.0, "#D55E00")]


def epoch_psd(files, days, target, n_avg=6):
    idx = int(np.argmin(np.abs(days - target)))
    lo = max(idx - n_avg // 2, 0)
    psds = []
    for f in files[lo : lo + n_avg]:
        x = np.loadtxt(f)[:, 0]
        freqs, p = welch_psd(x, fs=FS, nseg=NSEG)
        psds.append(p)
    return freqs, np.mean(psds, axis=0)


def main():
    files = sorted((ROOT / "data/ims/2nd_test").iterdir())
    from datetime import datetime
    t = [datetime.strptime(f.name, "%Y.%m.%d.%H.%M.%S") for f in files]
    days = np.array([(ti - t[0]).total_seconds() / 86400 for ti in t])

    curves = [(label, *epoch_psd(files, days, d)[1:2], c) for label, d, c in EPOCHS]
    freqs = epoch_psd(files, days, 1.0)[0]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7))
    for (label, p, c) in curves:
        ax1.semilogy(freqs, p, color=c, lw=0.8, label=label)
        ax2.semilogy(freqs, p, color=c, lw=1.0, label=label)
    ax1.set_xlim(0, FS / 2)
    ax1.set_xlabel("frequency [Hz]")
    ax1.set_title("B1 Welch PSD (nseg=2048, Hann, 6-snapshot average) - full band",
                  fontsize=10, loc="left")
    for k in range(1, 4):
        ax2.axvline(k * BPFO, color="#555555", lw=0.7, ls="--", zorder=1)
        ax2.annotate(f"{k}xBPFO", xy=(k * BPFO, 1), xycoords=("data", "axes fraction"),
                     xytext=(2, -10), textcoords="offset points", fontsize=7, color="#555555")
        for s in (-1, 1):
            ax2.axvline(k * BPFO + s * FR, color="#bbbbbb", lw=0.5, ls=":", zorder=1)
    ax2.set_xlim(0, 800)
    ax2.set_xlabel("frequency [Hz]")
    ax2.set_title("zoom 0-800 Hz  (dashed = declared BPFO harmonics, dotted = +/-fr sideband positions)",
                  fontsize=10, loc="left")
    for ax in (ax1, ax2):
        ax.set_ylabel("PSD [g$^2$/Hz]")
        ax.grid(True, color="#eeeeee", lw=0.5)
        ax.spines[["top", "right"]].set_visible(False)
    ax1.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    out = ROOT / "docs/figs/C2-b1-spectra.png"
    fig.savefig(out, dpi=150)
    print("figure:", out.name)

    # 定量: 広帯域レベルと BPFO 線
    df = freqs[1] - freqs[0]
    band = (freqs > 2000) & (freqs < 6000)
    print("\nbroadband level 2-6 kHz (median PSD):")
    ref = None
    for label, p, _ in curves:
        med = float(np.median(p[band]))
        ref = ref or med
        print(f"  {label:22s} {med:.3e}  (x{med / ref:.1f} vs baseline)")

    print("\nBPFO line search 225-245 Hz (peak freq, prominence vs local floor):")
    zone = (freqs >= 225) & (freqs <= 245)
    floor_zone = (freqs >= 180) & (freqs <= 300)
    for label, p, _ in curves:
        kpk = np.argmax(p[zone])
        fpk = freqs[zone][kpk]
        line = float(np.sum(p[zone][max(kpk - 1, 0) : kpk + 2]) * df)  # +/-1 bin 積分(リーケージ対策)
        floor = float(np.median(p[floor_zone]) * 3 * df)
        print(f"  {label:22s} peak at {fpk:6.1f} Hz  line/floor = {line / floor:7.2f}")

    print("\nsideband check at day 5.0: PSD(BPFO+/-fr) vs local floor (should be ~1 if no sideband):")
    p50 = curves[2][1]
    fl = float(np.median(p50[floor_zone]))
    for f0 in (BPFO - FR, BPFO + FR):
        k = int(round(f0 / df))
        val = float(p50[k - 1 : k + 2].sum() / 3)
        print(f"  {f0:6.1f} Hz: ratio = {val / fl:6.2f}")


if __name__ == "__main__":
    sys.exit(main() or 0)
