"""C1: IMS Set2 全期間(984スナップショット×4ch)の時間領域指標トレンド。

固定するもの: 指標定義(lib.indicators、C1で照合済み)、fs=20.48kHz(C0合意)、
対象=Set2全ファイル・全4ch。変えるもの: 時間(スナップショット順)。
出力: data/c1_set2_indicators.npz(キャッシュ)、docs/figs/C1-set2-trend.png。
"""

import sys
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from lib.indicators import crest_factor, kurtosis, rms

ROOT = Path(__file__).resolve().parents[1]
IMS2 = ROOT / "data/ims/2nd_test"
CACHE = ROOT / "data/c1_set2_indicators.npz"
FIG = ROOT / "docs/figs/C1-set2-trend.png"

# Okabe-Ito(CVD安全)。注目系列=軸受1(外輪故障ラベル)を強調、健全3chは抑える
COLORS = {1: "#D55E00", 2: "#0072B2", 3: "#009E73", 4: "#999999"}
LW = {1: 1.6, 2: 0.9, 3: 0.9, 4: 0.9}


def compute():
    files = sorted(IMS2.iterdir())
    t, out = [], {k: {ch: [] for ch in range(4)} for k in ("rms", "kurt", "cf")}
    for i, f in enumerate(files):
        x = np.loadtxt(f)
        t.append(datetime.strptime(f.name, "%Y.%m.%d.%H.%M.%S"))
        for ch in range(4):
            out["rms"][ch].append(rms(x[:, ch]))
            out["kurt"][ch].append(kurtosis(x[:, ch]))
            out["cf"][ch].append(crest_factor(x[:, ch]))
        if (i + 1) % 200 == 0:
            print(f"{i + 1}/{len(files)}", flush=True)
    t0 = t[0]
    days = np.array([(ti - t0).total_seconds() / 86400 for ti in t])
    np.savez(
        CACHE,
        days=days,
        **{f"{k}_ch{ch + 1}": np.array(v[ch]) for k, v in out.items() for ch in range(4)},
    )
    print("cached:", CACHE.name)


def _dodged_positions(vals, ax, min_px=13):
    """端ラベルのy位置(データ座標)を、表示ピクセルで min_px 以上離す。"""
    order = np.argsort(vals)
    pix = np.array([ax.transData.transform((0.0, v))[1] for v in vals], dtype=float)
    for a, b in zip(order[:-1], order[1:]):
        if pix[b] - pix[a] < min_px:
            pix[b] = pix[a] + min_px
    return {i: ax.transData.inverted().transform((0.0, pix[i]))[1] for i in range(len(vals))}


def plot():
    z = np.load(CACHE)
    days = z["days"]
    panels = [
        ("rms", "RMS [g]", None),
        ("kurt", "Kurtosis (Pearson)", 3.0),
        ("cf", "Crest factor", None),
    ]
    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    for ax, (key, ylabel, refline) in zip(axes, panels):
        if refline is not None:
            ax.axhline(refline, color="#bbbbbb", lw=0.8, ls="--", zorder=1)
            ax.annotate("Gaussian = 3", xy=(0.01, refline), xycoords=("axes fraction", "data"),
                        va="bottom", fontsize=8, color="#888888")
        for ch in (2, 3, 4, 1):  # 軸受1を最後に描いて前面へ
            y = z[f"{key}_ch{ch}"]
            ax.plot(days, y, color=COLORS[ch], lw=LW[ch], zorder=3 if ch == 1 else 2,
                    label=f"B{ch}")
        ax.set_ylabel(ylabel)
        ax.grid(True, color="#eeeeee", lw=0.6)
        ax.spines[["top", "right"]].set_visible(False)
        fig.canvas.draw()  # transData を確定させてからラベル位置を計算
        finals = np.array([z[f"{key}_ch{ch}"][-1] for ch in (1, 2, 3, 4)])
        ypos = _dodged_positions(finals, ax)
        for i, ch in enumerate((1, 2, 3, 4)):
            ax.annotate(f"B{ch}", xy=(days[-1], ypos[i]), xytext=(5, 0),
                        textcoords="offset points", va="center", fontsize=8,
                        color=COLORS[ch], fontweight="bold" if ch == 1 else "normal",
                        annotation_clip=False)
    handles, labels = axes[0].get_legend_handles_labels()
    order = [labels.index(f"B{ch}") for ch in (1, 2, 3, 4)]
    axes[0].legend([handles[i] for i in order],
                   ["B1 (failed: outer race)", "B2", "B3", "B4"],
                   loc="upper left", fontsize=8, frameon=False, ncol=4)
    axes[0].set_title(
        "IMS Set2 run-to-failure: time-domain indicators, all 4 bearings\n"
        "(B1 = outer-race failure label; 984 snapshots, 10 min apart, 2004-02-12 to 02-19)",
        fontsize=10, loc="left",
    )
    axes[-1].set_xlabel("days since start")
    fig.align_ylabels(axes)
    fig.tight_layout()
    fig.savefig(FIG, dpi=150)
    print("figure:", FIG.name)


if __name__ == "__main__":
    if not CACHE.exists() or "--recompute" in sys.argv:
        compute()
    plot()
