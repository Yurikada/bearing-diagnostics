"""C3: envelope.py の照合。

1) analytic の手トレース (N=4 cos) 突き合わせ。
2) analytic <-> scipy.signal.hilbert (乱数・実データ)。丸め精度一致を予測。
3) 合成 known-answer: 3kHz 搬送波を 107.4 Hz のパルス列で AM 変調した信号の
   エンベロープスペクトルに 107.4 Hz の線(と高調波)が立つこと。
   生スペクトルには 107.4 Hz の線が「立たない」ことも同時に確認する。
"""

import sys
from pathlib import Path

import numpy as np
import scipy.signal

sys.path.insert(0, str(Path(__file__).parent))
from lib.envelope import analytic, envelope_spectrum

ROOT = Path(__file__).resolve().parents[1]


def main():
    ok = True

    z = analytic([1.0, 0.0, -1.0, 0.0])
    print("trace: analytic([1,0,-1,0]) =", np.round(z, 12), " |z| =", np.round(np.abs(z), 12))

    rng = np.random.default_rng(3)
    sigs = {"randn(4096)": rng.normal(size=4096),
            "randn(4095) odd N": rng.normal(size=4095)}
    x = np.loadtxt(sorted((ROOT / "data/ims/2nd_test").iterdir())[0])[:, 0]
    sigs["IMS Set2 first ch1"] = x
    for name, s in sigs.items():
        ours = analytic(s)
        ref = scipy.signal.hilbert(s)
        r = float(np.abs(ours - ref).max() / np.abs(ref).max())
        print(f"analytic vs scipy.hilbert  rel_max = {r:.2e}  ({name})")
        ok &= r < 1e-12

    # 合成 known-answer
    fs, T = 12000.0, 10.0
    N = int(fs * T)
    t = np.arange(N) / fs
    fd, fc = 107.4, 3000.0
    impacts = (np.sin(np.pi * fd * t) ** 20)          # fd Hz の鈍いパルス列
    xsyn = impacts * np.sin(2 * np.pi * fc * t) + 0.05 * rng.normal(size=N)
    freqs, amp = envelope_spectrum(xsyn, fs, 2500.0, 3500.0)
    df = freqs[1] - freqs[0]
    k = int(round(fd / df))
    line = amp[k - 1 : k + 2].sum()
    floor = np.median(amp[int(50 / df) : int(500 / df)]) * 3
    print(f"synthetic AM: envelope line at {freqs[np.argmax(amp[int(20/df):int(500/df)])] + 20:.1f} Hz"
          f"  line/floor = {line / floor:.1f}")
    ok &= line / floor > 10
    Xraw = np.abs(np.fft.rfft(xsyn))
    raw_line = Xraw[k - 1 : k + 2].sum()
    raw_floor = np.median(Xraw[int(50 / df) : int(500 / df)]) * 3
    print(f"raw spectrum at {fd} Hz: line/floor = {raw_line / raw_floor:.2f} (should be ~1: no line)")

    print("ALL OK" if ok else "SOME CHECK FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
