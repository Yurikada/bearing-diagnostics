"""C7 転移テスト: IMS Set1 (8ch, B3内輪+B4転動体ラベル) への同一パイプライン適用。

合意 (2026-07-27):
- 主対象 ch5-8 (B3/B4 の x/y)、健全対照 ch1-4。因果帯域を主、回顧を参考。
- ターゲット窓: B3 = BPFI 296.9±2%、B4 = 2xBSF 279.8±2%。側帯: ±33.33 / ±14.77 Hz。
- 帰属: 周波数分離型 max-channel (各ターゲット窓の線強度を ch 間比較)。
- ギャップ >30min をまたぐ連続判定は切る (C5 のロック仕様と同じ理由)。
- 床規則 rms < 0.01 g。事前の予想は scheme.md / C7 文書に記録済み。

stage 1 (このファイルの compute_indicators): 8ch の RMS・尖度トレンド。
"""

import sys
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from lib.indicators import kurtosis, rms

ROOT = Path(__file__).resolve().parents[1]
SET1 = ROOT / "data/ims/1st_test"
CACHE = ROOT / "data/c7_set1_indicators.npz"
FIG = ROOT / "docs/figs/C7-set1-trend.png"
RMS_FLOOR = 0.01


def compute_indicators():
    files = sorted(SET1.iterdir())
    t0 = datetime.strptime(files[0].name, "%Y.%m.%d.%H.%M.%S")
    rows = []
    for i, f in enumerate(files):
        x = np.loadtxt(f)
        day = (datetime.strptime(f.name, "%Y.%m.%d.%H.%M.%S") - t0).total_seconds() / 86400
        row = [day]
        for ch in range(8):
            v = x[:, ch]
            if rms(v) < RMS_FLOOR:
                row += [np.nan, np.nan]
            else:
                row += [rms(v), kurtosis(v)]
        rows.append(row)
        if (i + 1) % 300 == 0:
            print(f"{i + 1}/{len(files)}", flush=True)
    arr = np.array(rows)
    np.savez(CACHE, days=arr[:, 0],
             **{f"rms_ch{c + 1}": arr[:, 1 + 2 * c] for c in range(8)},
             **{f"kurt_ch{c + 1}": arr[:, 2 + 2 * c] for c in range(8)})
    print("cached:", CACHE.name)


BEARING_COLORS = {1: "#999999", 2: "#0072B2", 3: "#D55E00", 4: "#009E73"}


def plot():
    z = np.load(CACHE)
    d = z["days"]
    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
    for ax, key, ylabel in ((axes[0], "rms", "RMS [g]"), (axes[1], "kurt", "Kurtosis (Pearson)")):
        if key == "kurt":
            ax.axhline(3.0, color="#bbbbbb", lw=0.8, ls="--")
        for ch in range(1, 9):
            b = (ch + 1) // 2                       # ch1,2->B1 ... ch7,8->B4
            style = "-" if ch % 2 == 1 else ":"
            ax.plot(d, z[f"{key}_ch{ch}"], style, lw=1.0 if b in (3, 4) else 0.7,
                    color=BEARING_COLORS[b], label=f"ch{ch} (B{b}{'x' if ch % 2 else 'y'})")
        ax.set_ylabel(ylabel)
        ax.grid(True, color="#eeeeee", lw=0.5)
        ax.spines[["top", "right"]].set_visible(False)
    axes[1].set_yscale("log")
    axes[0].legend(fontsize=7, frameon=False, ncol=4)
    axes[1].set_xlabel("days since start (wall time; 24 gaps >30 min incl. 6.2-day gap)")
    axes[0].set_title("IMS Set1: time-domain indicators, 8 channels "
                      "(labels: B3 inner race + B4 roller at end of test)", fontsize=10, loc="left")
    fig.tight_layout()
    fig.savefig(FIG, dpi=150)
    print("figure:", FIG.name)


CACHE_ENV = ROOT / "data/c7_set1_envelope.npz"
W_B3 = (291, 303)      # BPFI 296.9 +/- 2%
W_B4 = (274, 285)      # 2xBSF 279.8 +/- 2%
RETRO_BAND = (3000.0, 6500.0)
CAUSAL_LOS = None      # set below


def _line(amp, w):
    """窓 w 内ピークの±1bin積分 / レイリー床(50-500Hz中央値ベース)。"""
    lo, hi = w
    k = lo + int(np.argmax(amp[lo : hi + 1]))
    line = float(amp[k - 1 : k + 2].sum())
    floor = float(np.median(amp[50:500])) * 3
    return line / floor, k


def compute_envelope():
    """全8ch: 因果帯域(各ch独立スキャン+中央値5)で B3窓/B4窓の線床比・位置。
    ch5-8 は回顧帯域(3-6.5kHz, lookahead参考)も併記。"""
    from c4_trend import env_amp_from_mask
    from lib.indicators import kurtosis
    causal_los = np.arange(500.0, 9300.0, 500.0)
    files = sorted(SET1.iterdir())
    t0 = datetime.strptime(files[0].name, "%Y.%m.%d.%H.%M.%S")
    freqs = np.fft.rfftfreq(20480, d=1.0 / 20480.0)
    keep_retro = (freqs >= RETRO_BAND[0]) & (freqs <= RETRO_BAND[1])
    recent = {ch: [] for ch in range(8)}
    rows = []
    ncol_per_ch = 5   # causal: lo, lfB3, posB3, lfB4, posB4
    for i, f in enumerate(files):
        dat = np.loadtxt(f)
        day = (datetime.strptime(f.name, "%Y.%m.%d.%H.%M.%S") - t0).total_seconds() / 86400
        row = [day]
        for ch in range(8):
            x = dat[:, ch]
            if rms(x) < RMS_FLOOR:
                row += [np.nan] * ncol_per_ch
                if ch >= 4:
                    row += [np.nan] * 4
                continue
            x = x - x.mean()
            X = np.fft.rfft(x)
            kurts = [kurtosis(np.fft.irfft(np.where((freqs >= lo) & (freqs < lo + 1000.0), X, 0.0),
                                           n=20480)) for lo in causal_los]
            recent[ch].append(float(causal_los[int(np.argmax(kurts))]))
            lo_used = float(np.median(recent[ch][-5:]))
            m = (freqs >= lo_used) & (freqs < lo_used + 1000.0)
            amp = env_amp_from_mask(X, m)
            lf3, p3 = _line(amp, W_B3)
            lf4, p4 = _line(amp, W_B4)
            row += [lo_used, lf3, p3, lf4, p4]
            if ch >= 4:
                ampr = env_amp_from_mask(X, keep_retro)
                r3, rp3 = _line(ampr, W_B3)
                r4, rp4 = _line(ampr, W_B4)
                row += [r3, rp3, r4, rp4]
        rows.append(row)
        if (i + 1) % 200 == 0:
            print(f"{i + 1}/{len(files)}", flush=True)
    arr = np.array(rows)
    out = {"days": arr[:, 0]}
    col = 1
    for ch in range(8):
        for name in ("clo", "clf3", "cpos3", "clf4", "cpos4"):
            out[f"{name}_ch{ch + 1}"] = arr[:, col]; col += 1
        if ch >= 4:
            for name in ("rlf3", "rpos3", "rlf4", "rpos4"):
                out[f"{name}_ch{ch + 1}"] = arr[:, col]; col += 1
    np.savez(CACHE_ENV, **out)
    print("cached:", CACHE_ENV.name)


if __name__ == "__main__":
    if not CACHE.exists() or "--recompute" in sys.argv:
        compute_indicators()
    if "--envelope" in sys.argv:
        compute_envelope()
    else:
        plot()
