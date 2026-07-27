"""C1: 自前実装 indicators.py と scipy/numpy の突き合わせ。

照合対象:
- rms          <-> numpy の母集団標準偏差 np.std(x) (平均除去後の定義なので一致するはず)
- kurtosis     <-> scipy.stats.kurtosis (既定値のまま / 規約明示 fisher=False, bias=True の両方)
- crest_factor <-> numpy 手組み (max|x-mean| / std)

信号: 合成ガウス + CWRU 実スナップショット + IMS Set2 実スナップショット。
判定閾値 rtol = 1e-12 (1/N vs 1/(N-1) の取り違え ~1e-4 を検出できる厳しさ)。
"""

import sys
from pathlib import Path

import numpy as np
import scipy.io
import scipy.stats

sys.path.insert(0, str(Path(__file__).parent))
from lib.indicators import crest_factor, kurtosis, rms

ROOT = Path(__file__).resolve().parents[1]
RTOL = 1e-12


def rel(a, b):
    return abs(a - b) / max(abs(a), abs(b))


def load_signals():
    signals = {}
    rng = np.random.default_rng(42)
    signals["gauss(N=20480)"] = rng.normal(0.0, 1.0, 20480)

    mat = scipy.io.loadmat(ROOT / "data/cwru/97.mat")
    key = [k for k in mat if k.endswith("_DE_time")][0]
    signals[f"CWRU 97.mat {key}"] = mat[key].ravel()[:20480]

    ims_dir = ROOT / "data/ims/2nd_test"
    first = sorted(p.name for p in ims_dir.iterdir())[0]
    signals[f"IMS Set2 {first} ch1"] = np.loadtxt(ims_dir / first)[:, 0]
    return signals


def main():
    all_ok = True
    for name, x in load_signals().items():
        print(f"--- {name}  (N={x.size}) ---")
        rows = [
            ("rms            vs np.std(x)               ", rms(x), float(np.std(x))),
            ("kurtosis       vs scipy DEFAULT           ", kurtosis(x), float(scipy.stats.kurtosis(x))),
            ("kurtosis       vs scipy fisher=F, bias=T  ", kurtosis(x),
             float(scipy.stats.kurtosis(x, fisher=False, bias=True))),
            ("crest_factor   vs numpy manual            ", crest_factor(x),
             float(np.max(np.abs(x - x.mean())) / np.std(x))),
        ]
        for label, ours, ref in rows:
            r = rel(ours, ref)
            ok = r < RTOL
            if "DEFAULT" not in label:  # DEFAULT 行は規約差の展示であり合否対象外
                all_ok &= ok
            print(f"{label} ours={ours:+.12f} ref={ref:+.12f} rel={r:.2e} {'MATCH' if ok else 'MISMATCH'}")
    print("NOTE: scipy DEFAULT row is expected to disagree - diagnosis is part of C1.")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
