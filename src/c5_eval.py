"""C5: しきい値設計と評価。

合意済み (2026-07-27):
- 検出器: RMS / 尖度 / エンベロープ線床比(回顧・因果) / FLI(同一3bin窓にm点連続)。
- しきい値: kσ (k=2..10, 0.5刻み) と分位点 (baseline 99% / 99.9%)。発報=3点連続。
- 誤報: B2/B4 (壊れていない) に同じ検出器を適用し発報クラスタ数を数える。
  エンベロープは回顧帯域 3-6.5kHz・236Hz 窓のまま (帰属誤りの検査)。
- 故障基準時刻: rms_ch1 > 10x baseline 中央値の最初のスナップショット。
- lookahead: retro=あり / causal・RMS・尖度・FLI(retro帯域使用なので FLI もあり)。
"""

import sys
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from c4_trend import env_amp_from_mask, line_over_floor, RETRO_BAND, RMS_FLOOR
from lib.indicators import rms

ROOT = Path(__file__).resolve().parents[1]
FS = 20480.0
CACHE24 = ROOT / "data/c5_ch24_envtrend.npz"
FIG = ROOT / "docs/figs/C5-operating-curves.png"
KS = np.arange(2.0, 10.5, 0.5)
MS = range(3, 9)


def compute_ch24():
    files = sorted((ROOT / "data/ims/2nd_test").iterdir())
    t0 = datetime.strptime(files[0].name, "%Y.%m.%d.%H.%M.%S")
    freqs = np.fft.rfftfreq(20480, d=1.0 / FS)
    keep = (freqs >= RETRO_BAND[0]) & (freqs <= RETRO_BAND[1])
    rows = []
    for i, f in enumerate(files):
        dat = np.loadtxt(f)
        day = (datetime.strptime(f.name, "%Y.%m.%d.%H.%M.%S") - t0).total_seconds() / 86400
        row = [day]
        for ch in (1, 3):                     # ch2, ch4 (0-index)
            x = dat[:, ch]
            if rms(x) < RMS_FLOOR:
                row += [np.nan, np.nan]
                continue
            x = x - x.mean()
            lf, _, pos = line_over_floor(env_amp_from_mask(np.fft.rfft(x), keep))
            row += [lf, pos]
        rows.append(row)
        if (i + 1) % 200 == 0:
            print(f"{i + 1}/{len(files)}", flush=True)
    arr = np.array(rows)
    np.savez(CACHE24, days=arr[:, 0], lf_ch2=arr[:, 1], pos_ch2=arr[:, 2],
             lf_ch4=arr[:, 3], pos_ch4=arr[:, 4])
    print("cached:", CACHE24.name)


def detect_day(days, y, thr):
    """3点連続で thr 超過する最初の day (なければ None)。NaN は不検出扱い。"""
    over = np.where(np.isfinite(y), y > thr, False)
    for i in range(len(over) - 2):
        if over[i] and over[i + 1] and over[i + 2]:
            return float(days[i])
    return None


def alarm_cluster_days(days, y, thr):
    """3点連続超過の塊の開始 day の配列。"""
    over = np.where(np.isfinite(y), y > thr, False)
    det = over[:-2] & over[1:-1] & over[2:]
    starts = det & ~np.concatenate([[False], det[:-1]])
    return days[: len(starts)][starts]


def alarm_clusters(y, thr, days=None):
    """3点連続超過の塊の数。"""
    over = np.where(np.isfinite(y), y > thr, False)
    det = over[:-2] & over[1:-1] & over[2:]
    return int(np.sum(det & ~np.concatenate([[False], det[:-1]])))


def lock_runs(pos):
    """各時点で終わる「max-min<=2bin」の最長連続長。NaN でリセット。"""
    L = np.zeros(len(pos), dtype=int)
    for i in range(len(pos)):
        if not np.isfinite(pos[i]):
            continue
        lo = hi = pos[i]
        n = 1
        j = i - 1
        while j >= 0 and np.isfinite(pos[j]):
            lo2, hi2 = min(lo, pos[j]), max(hi, pos[j])
            if hi2 - lo2 > 2:
                break
            lo, hi, n = lo2, hi2, n + 1
            j -= 1
        L[i] = n
    return L


def fli_detect_day(days, pos, m):
    L = lock_runs(pos)
    idx = np.where(L >= m)[0]
    return float(days[max(idx[0] - m + 1, 0)]) if idx.size else None


def fli_cluster_days(days, pos, m):
    L = lock_runs(pos)
    det = L >= m
    starts = det & ~np.concatenate([[False], det[:-1]])
    return days[starts]


def main():
    z1 = np.load(ROOT / "data/c1_set2_indicators.npz")
    z4 = np.load(ROOT / "data/c4_set2_envtrend.npz")
    z24 = np.load(CACHE24)
    d = z1["days"]

    rms1 = z1["rms_ch1"]
    base_med = np.median(rms1[d < 3.0])
    # 当初宣言の 10x は実データ最大 0.725g が届かず不成立 -> 5x に改訂 (2026-07-27, 開示済み)
    fail_idx = np.where(np.where(np.isfinite(rms1), rms1, 0) > 5 * base_med)[0][0]
    t_fail = float(d[fail_idx])
    print(f"failure reference = day {t_fail:.3f} (rms > 5x baseline median {base_med:.4f};"
          " NOTE: revised from declared 10x which no snapshot reached)")

    detectors = {  # name -> (B1 series, [(健全ch名, series), ...], lookahead)
        "RMS": (rms1, [("B2", z1["rms_ch2"]), ("B4", z1["rms_ch4"])], "no"),
        "kurtosis": (z1["kurt_ch1"], [("B2", z1["kurt_ch2"]), ("B4", z1["kurt_ch4"])], "no"),
        "env retro": (z4["retro_lf"], [("B2", z24["lf_ch2"]), ("B4", z24["lf_ch4"])], "YES"),
        "env causal": (z4["causal_lf"], [("B2", z24["lf_ch2"]), ("B4", z24["lf_ch4"])], "no"),
    }

    print("\n== k-sigma sweep: lead time [d] on B1 / false-alarm clusters on B2+B4 ==")
    curves = {}
    for name, (y1, healthy, la) in detectors.items():
        base = np.isfinite(y1) & (d < 3.0)
        mu, sd = y1[base].mean(), y1[base].std()
        leads, fas = [], []
        for k in KS:
            td = detect_day(d, y1, mu + k * sd)
            leads.append(t_fail - td if td is not None else np.nan)
            fa = 0
            for _, yh in healthy:
                bh = np.isfinite(yh) & (d < 3.0)
                fa += alarm_clusters(yh, yh[bh].mean() + k * yh[bh].std())
            fas.append(fa)
        curves[name] = (np.array(leads), np.array(fas), la)
        k5 = np.where(KS == 5.0)[0][0]
        early = late = 0
        for _, yh in healthy:
            bh = np.isfinite(yh) & (d < 3.0)
            cd = alarm_cluster_days(d, yh, yh[bh].mean() + 5.0 * yh[bh].std())
            early += int(np.sum(cd < 4.5))
            late += int(np.sum(cd >= 4.5))
        print(f"  {name:10s} (lookahead {la:3s}): k=5 -> lead {leads[k5]:.2f} d,"
              f" FA {fas[k5]} (early<4.5d: {early} / late: {late})"
              f" | k=10 -> lead {leads[-1]:.2f} d, FA {fas[-1]}")

    print("\n== quantile thresholds (baseline 99% / 99.9%) ==")
    for name, (y1, healthy, _) in detectors.items():
        base = np.isfinite(y1) & (d < 3.0)
        for q in (99.0, 99.9):
            thr = np.percentile(y1[base], q)
            td = detect_day(d, y1, thr)
            fa = sum(alarm_clusters(yh, np.percentile(yh[np.isfinite(yh) & (d < 3.0)], q))
                     for _, yh in healthy)
            lead = t_fail - td if td is not None else np.nan
            print(f"  {name:10s} q{q}: lead {lead:.2f} d, FA {fa}")

    print("\n== FLI (same-3bin-window run of m) ==")
    fli = {"B1": z4["retro_pos"], "B2": z24["pos_ch2"], "B4": z24["pos_ch4"]}
    fli_curve = []
    for m in MS:
        td = fli_detect_day(d, fli["B1"], m)
        lead = t_fail - td if td is not None else np.nan
        cd = np.concatenate([fli_cluster_days(d, fli["B2"], m), fli_cluster_days(d, fli["B4"], m)])
        fa, fe, fl = len(cd), int(np.sum(cd < 4.5)), int(np.sum(cd >= 4.5))
        # B1 側: 最初のロックが「本物の劣化(3.69日〜)」か「ベースライン期の偶然」かを開示
        fli_curve.append((m, lead, fa))
        print(f"  m={m}: B1 first lock at day {td if td is not None else float('nan'):.2f}"
              f" (lead {lead:.2f} d), FA(B2+B4) {fa} (early: {fe} / late: {fl})")

    print("\n== k-sigma calibration on skewed baselines (empirical exceedance) ==")
    for name, (y1, _, _) in detectors.items():
        b = y1[np.isfinite(y1) & (d < 3.0)]
        mu, sd = b.mean(), b.std()
        for k in (2, 3):
            emp = float(np.mean(b > mu + k * sd))
            gauss = {2: 2.28e-2, 3: 1.35e-3}[k]
            print(f"  {name:10s} P(x>mu+{k}sigma): empirical {emp:.4f} vs Gaussian {gauss:.4f}")

    colors = {"RMS": "#999999", "kurtosis": "#009E73", "env retro": "#D55E00", "env causal": "#0072B2"}
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
    for name, (leads, fas, la) in curves.items():
        ax1.plot(KS, leads, "-o", ms=3, lw=1.2, color=colors[name],
                 label=name + (" (lookahead)" if la == "YES" else ""))
        ax2.plot(KS, fas, "-o", ms=3, lw=1.2, color=colors[name])
    m_arr, lead_arr, fa_arr = zip(*fli_curve)
    ax1.plot([KS[0] + 0.1], [lead_arr[0]], "*", ms=12, color="#CC79A7")
    ax1.annotate(f"FLI m=3 (FA {fa_arr[0]})", xy=(KS[0] + 0.2, lead_arr[0]), fontsize=8,
                 color="#CC79A7", va="bottom")
    ax1.set_xlabel("k (threshold = baseline mean + k sigma)")
    ax1.set_ylabel("lead time before failure [days]")
    ax1.legend(fontsize=8, frameon=False)
    ax2.set_xlabel("k")
    ax2.set_ylabel("false-alarm clusters on healthy B2+B4")
    for ax in (ax1, ax2):
        ax.grid(True, color="#eeeeee", lw=0.6)
        ax.spines[["top", "right"]].set_visible(False)
    ax1.set_title("lead time vs threshold (B1)", fontsize=10, loc="left")
    ax2.set_title("wrong-attribution alarms (B2/B4)", fontsize=10, loc="left")
    fig.tight_layout()
    fig.savefig(FIG, dpi=150)
    print("figure:", FIG.name)


if __name__ == "__main__":
    if not CACHE24.exists() or "--recompute" in sys.argv:
        compute_ch24()
    main()
