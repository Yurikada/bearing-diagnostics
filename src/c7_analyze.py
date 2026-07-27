"""C7 stage 3: Set1 エンベロープ結果の判定材料。

規則 (合意済み):
- ベースライン = day 7.11-8.64 の100ファイル (本人設計の二段構え、同等性チェック済み)。
- しきい値 = ベースライン q99.9 (C5 の env 推奨) と +5σ を併記。
- 発報確定 = 3点連続。ただし >30min のギャップをまたぐ連続は切る。
- 帰属 = 周波数分離型 max-channel (窓ごとに全8chの線/床を比較)。
- 側帯検査: B3 = BPFI±33.33、B4 = 2xBSF±14.77 (指定エポックの実スペクトルで測る)。
"""

import sys
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from c4_trend import env_amp_from_mask
from c7_set1 import SET1, W_B3, W_B4, _line

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "docs/figs/C7-set1-envelope.png"
BPFI, FR, BSF2, FTF = 296.9, 33.33, 279.8, 14.77


def gap_aware_confirm(days, y, thr):
    """3点連続超過(ギャップ>30minで連続を切る)の確定時刻。"""
    over = np.where(np.isfinite(y), y > thr, False)
    run = 0
    for i in range(len(over)):
        if i > 0 and (days[i] - days[i - 1]) > 30 / 1440:
            run = 0
        run = run + 1 if over[i] else 0
        if run >= 3:
            return float(days[i])
    return None


def clusters(days, y, thr):
    over = np.where(np.isfinite(y), y > thr, False)
    run = 0
    out = []
    prev_fired = False
    for i in range(len(over)):
        if i > 0 and (days[i] - days[i - 1]) > 30 / 1440:
            run = 0
        run = run + 1 if over[i] else 0
        fired = run >= 3
        if fired and not prev_fired:
            out.append(float(days[i]))
        prev_fired = fired
    return np.array(out)


def main():
    z = np.load(ROOT / "data/c7_set1_envelope.npz")
    d = z["days"]
    base = (d >= 7.11) & (d <= 8.64)

    def thr_of(y, q=99.9):
        b = y[np.isfinite(y) & base]
        return float(np.percentile(b, q)), float(b.mean() + 5 * b.std())

    print("== onsets (baseline day 7.11-8.64, gap-aware 3-consecutive confirm) ==")
    for label, key, chs in (("B3 window (BPFI)", "clf3", (5, 6)), ("B4 window (2xBSF)", "clf4", (7, 8))):
        for ch in chs:
            y = z[f"{key}_ch{ch}"]
            tq, tk = thr_of(y)
            cq = gap_aware_confirm(d, y, tq)
            ck = gap_aware_confirm(d, y, tk)
            print(f"  {label} ch{ch}: q99.9 -> {cq if cq is None else round(cq, 2)}"
                  f"  | k=5 -> {ck if ck is None else round(ck, 2)}")

    print("\n== false alarms on healthy B1/B2 (ch1-4, own baselines, q99.9) ==")
    for key, label in (("clf3", "B3win"), ("clf4", "B4win")):
        for ch in (1, 2, 3, 4):
            y = z[f"{key}_ch{ch}"]
            tq, _ = thr_of(y)
            cd = clusters(d, y, tq)
            if len(cd):
                print(f"  {label} ch{ch}: {len(cd)} clusters (first at day {cd[0]:.1f},"
                      f" early<20d: {(cd < 20).sum()})")
            else:
                print(f"  {label} ch{ch}: 0 clusters")

    print("\n== frequency-separated attribution (post-onset epochs) ==")
    for key, label, own, epoch in (("clf3", "B3win", (5, 6), (31.5, 34.0)),
                                   ("clf4", "B4win", (7, 8), (28.0, 31.0))):
        m = (d >= epoch[0]) & (d <= epoch[1])
        stack = np.vstack([np.where(np.isfinite(z[f"{key}_ch{c}"]), z[f"{key}_ch{c}"], -np.inf)
                           for c in range(1, 9)])
        win = np.argmax(stack[:, m], axis=0) + 1
        frac_own = float(np.isin(win, own).mean())
        vals, counts = np.unique(win, return_counts=True)
        top = sorted(zip(counts.tolist(), vals.tolist()), reverse=True)[:3]
        print(f"  {label} day {epoch[0]}-{epoch[1]}: own-ch({own}) wins {frac_own * 100:.1f}%"
              f"  top: {[(f'ch{v}', c) for c, v in top]}")

    print("\n== position lock (target vs measured mode, post-onset) ==")
    for key, ch, target, epoch in (("cpos3", 5, 297, (31.5, 34.0)), ("cpos4", 7, 280, (28.0, 31.0))):
        p = z[f"{key}_ch{ch}"]
        m = (d >= epoch[0]) & (d <= epoch[1]) & np.isfinite(p)
        vals, counts = np.unique(p[m].astype(int), return_counts=True)
        top = sorted(zip(counts.tolist(), vals.tolist()), reverse=True)[:3]
        print(f"  ch{ch} target {target}: n={int(m.sum())}  top: {[(v, c) for c, v in top]}")

    print("\n== sideband check on raw envelope spectra (6-file average) ==")
    files = sorted(SET1.iterdir())
    t0 = datetime.strptime(files[0].name, "%Y.%m.%d.%H.%M.%S")
    fdays = np.array([(datetime.strptime(f.name, "%Y.%m.%d.%H.%M.%S") - t0).total_seconds() / 86400
                      for f in files])
    freqs = np.fft.rfftfreq(20480, d=1.0 / 20480.0)

    def epoch_amp(target_day, ch, band):
        i = int(np.argmin(np.abs(fdays - target_day)))
        acc = None
        for f in files[max(i - 3, 0) : i + 3]:
            x = np.loadtxt(f)[:, ch - 1]
            x = x - x.mean()
            keep = (freqs >= band[0]) & (freqs <= band[1])
            a = env_amp_from_mask(np.fft.rfft(x), keep)
            acc = a if acc is None else acc + a
        return acc / 6

    def report_lines(amp, main_f, side_f, name):
        lf_main, pos = _line(amp, (int(main_f * 0.98), int(main_f * 1.02)))
        parts = [f"{name}: main {pos}Hz lf={lf_main:.1f}"]
        for s in (-1, +1):
            f0 = pos + s * side_f
            lf_s, ps = _line(amp, (int(round(f0 - 3)), int(round(f0 + 3))))
            parts.append(f"side{'-' if s < 0 else '+'} {ps}Hz lf={lf_s:.1f}")
        print("   " + "  ".join(parts))

    for day0, ch, main_f, side_f, name in ((32.5, 5, BPFI, FR, "B3(ch5) day32.5 BPFI+/-fr"),
                                           (30.3, 7, BSF2, FTF, "B4(ch7) day30.3 2xBSF+/-FTF")):
        lo = float(np.nanmedian(z[f"clo_ch{ch}"][(fdays >= day0 - 0.2) & (fdays <= day0 + 0.2)]))
        amp = epoch_amp(day0, ch, (lo, lo + 1000.0))
        report_lines(amp, main_f, side_f, name + f" [band {lo:.0f}-{lo + 1000:.0f}]")

    fig, axes = plt.subplots(3, 1, figsize=(11, 9), sharex=True)
    for ax, key, own, tgt_label in ((axes[0], "clf3", (5, 6), "B3 window (BPFI 296.9)"),
                                    (axes[1], "clf4", (7, 8), "B4 window (2xBSF 279.8)")):
        for ch in range(1, 9):
            b = (ch + 1) // 2
            bold = ch in own
            ax.semilogy(d, z[f"{key}_ch{ch}"], lw=1.2 if bold else 0.5,
                        color={1: "#999999", 2: "#0072B2", 3: "#D55E00", 4: "#009E73"}[b],
                        alpha=1.0 if bold else 0.45,
                        label=f"ch{ch}" if bold else None)
        ax.set_ylabel(f"{tgt_label}\nline/floor")
        ax.legend(fontsize=8, frameon=False, loc="upper left")
        ax.grid(True, color="#eeeeee", lw=0.5)
        ax.spines[["top", "right"]].set_visible(False)
    ax = axes[2]
    ax.plot(d, z["cpos3_ch5"], ".", ms=2, color="#D55E00", label="ch5 peak in B3 window")
    ax.plot(d, z["cpos4_ch7"], ".", ms=2, color="#009E73", label="ch7 peak in B4 window")
    ax.axhline(297, color="#D55E00", lw=0.7, ls="--")
    ax.axhline(280, color="#009E73", lw=0.7, ls="--")
    ax.set_ylabel("window peak freq [Hz]")
    ax.set_xlabel("days since start")
    ax.legend(fontsize=8, frameon=False)
    ax.grid(True, color="#eeeeee", lw=0.5)
    ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_title("IMS Set1: causal-band envelope line trends per target window",
                      fontsize=10, loc="left")
    fig.tight_layout()
    fig.savefig(FIG, dpi=150)
    print("figure:", FIG.name)


if __name__ == "__main__":
    sys.exit(main() or 0)
