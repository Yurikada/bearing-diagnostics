"""C4: IMS Set2 B1 のエンベロープ BPFO 線トレンド (回顧帯域 vs 因果帯域)。

合意済みパラメータ (2026-07-27):
- 線: 236Hz±2% 窓でピーク探索、ピーク±1bin 積分。床: 50-500Hz 中央値ベース。
- 回顧帯域: 3000-6500Hz 固定 (lookahead あり、と明記して報告)。
- 因果帯域: 1kHz 窓×500Hz 刻み (500-9200) の帯域内尖度最大 → 直近5点の中央値。
  その時点までの情報しか使わない (lookahead なし)。
- 除外: rms < 0.01 g (機械停止)。
出力: data/c4_set2_envtrend.npz、docs/figs/C4-envelope-trend.png
"""

import sys
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from lib.envelope import analytic
from lib.indicators import kurtosis, rms

ROOT = Path(__file__).resolve().parents[1]
FS = 20480.0
CACHE = ROOT / "data/c4_set2_envtrend.npz"
FIG = ROOT / "docs/figs/C4-envelope-trend.png"
BPFO_ZONE = (232.0, 241.0)          # 236 Hz ± 2%
FLOOR_ZONE = (50, 500)
RETRO_BAND = (3000.0, 6500.0)
CAUSAL_LOS = np.arange(500.0, 9300.0, 500.0)   # 1kHz 窓の下端
RMS_FLOOR = 0.01


def env_amp_from_mask(X, keep):
    """rfft 済み X に帯域マスク keep を掛け、包絡線の片側振幅スペクトルを返す。"""
    Xb = np.where(keep, X, 0.0)
    xb = np.fft.irfft(Xb, n=20480)
    env = np.abs(analytic(xb))
    env -= env.mean()
    return np.abs(np.fft.rfft(env)) * 2.0 / env.size


def line_over_floor(amp):
    """BPFO 窓内ピークを ±1bin 積分し、レイリー床に対する比と位置を返す。Δf=1Hz。"""
    lo, hi = int(BPFO_ZONE[0]), int(BPFO_ZONE[1])
    k = lo + int(np.argmax(amp[lo : hi + 1]))
    line = float(amp[k - 1 : k + 2].sum())
    floor = float(np.median(amp[FLOOR_ZONE[0] : FLOOR_ZONE[1]])) * 3
    return line / floor, float(line), k


def compute():
    files = sorted((ROOT / "data/ims/2nd_test").iterdir())
    t0 = datetime.strptime(files[0].name, "%Y.%m.%d.%H.%M.%S")
    freqs = np.fft.rfftfreq(20480, d=1.0 / FS)
    recent_los = []
    rows = []
    for i, f in enumerate(files):
        x = np.loadtxt(f)[:, 0]
        day = (datetime.strptime(f.name, "%Y.%m.%d.%H.%M.%S") - t0).total_seconds() / 86400
        if rms(x) < RMS_FLOOR:
            rows.append((day, *([np.nan] * 7)))
            continue
        x = x - x.mean()
        X = np.fft.rfft(x)

        keep = (freqs >= RETRO_BAND[0]) & (freqs <= RETRO_BAND[1])
        r_lf, r_abs, r_pos = line_over_floor(env_amp_from_mask(X, keep))

        # 因果帯域: 現時点のスナップショットだけで帯域内尖度を評価
        kurts = []
        for lo in CAUSAL_LOS:
            m = (freqs >= lo) & (freqs < lo + 1000.0)
            kurts.append(kurtosis(np.fft.irfft(np.where(m, X, 0.0), n=20480)))
        recent_los.append(float(CAUSAL_LOS[int(np.argmax(kurts))]))
        lo_used = float(np.median(recent_los[-5:]))   # 直近5点の中央値 (過去のみ)
        m = (freqs >= lo_used) & (freqs < lo_used + 1000.0)
        c_lf, c_abs, c_pos = line_over_floor(env_amp_from_mask(X, m))

        rows.append((day, r_lf, r_abs, r_pos, c_lf, c_abs, c_pos, lo_used))
        if (i + 1) % 200 == 0:
            print(f"{i + 1}/{len(files)}", flush=True)
    arr = np.array(rows, dtype=float)
    np.savez(CACHE, days=arr[:, 0], retro_lf=arr[:, 1], retro_abs=arr[:, 2],
             retro_pos=arr[:, 3], causal_lf=arr[:, 4], causal_abs=arr[:, 5],
             causal_pos=arr[:, 6], causal_lo=arr[:, 7])
    print("cached:", CACHE.name)


def plot():
    z = np.load(CACHE)
    d = z["days"]
    fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)

    ax = axes[0]
    ax.semilogy(d, z["retro_lf"], color="#D55E00", lw=0.9, label="retrospective band 3-6.5 kHz (lookahead)")
    ax.semilogy(d, z["causal_lf"], color="#0072B2", lw=0.9, label="causal band (online kurtogram, median-5)")
    ax.axhline(1.0, color="#bbbbbb", lw=0.8, ls="--")
    ax.set_ylabel("BPFO line / floor")
    ax.legend(fontsize=8, frameon=False, loc="upper left")

    ax = axes[1]
    ax.plot(d, z["retro_pos"], ".", ms=2, color="#D55E00")
    ax.plot(d, z["causal_pos"], ".", ms=2, color="#0072B2", alpha=0.5)
    ax.axhline(236, color="#555555", lw=0.7, ls="--")
    ax.axhline(237, color="#009E73", lw=0.7, ls=":")
    ax.set_ylabel("peak freq in 232-241 Hz window [Hz]")
    ax.set_ylim(231, 242)

    ax = axes[2]
    ax.plot(d, z["causal_lo"], ".", ms=2, color="#0072B2")
    ax.set_ylabel("causal band lower edge [Hz]")
    ax.set_xlabel("days since start")

    for ax in axes:
        for x0, ls in ((3.70, "--"), (4.79, ":")):
            ax.axvline(x0, color="#888888", lw=0.8, ls=ls)
        ax.grid(True, color="#eeeeee", lw=0.5)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_title("IMS Set2 B1: envelope BPFO line trend  "
                      "(vertical dashed = C1 RMS onset 3.70d, dotted = kurtosis onset 4.79d)",
                      fontsize=10, loc="left")
    fig.tight_layout()
    fig.savefig(FIG, dpi=150)
    print("figure:", FIG.name)


if __name__ == "__main__":
    if not CACHE.exists() or "--recompute" in sys.argv:
        compute()
    plot()
