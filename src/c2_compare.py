"""C2: spectral.py と numpy/scipy の突き合わせ + リーケージの定量。

予測(実装前宣言):
- dft_naive / fft_radix2 <-> numpy.fft: 丸めのみ(~1e-12相対)。
- hann <-> scipy.signal.get_window('hann', N) (fftbins=True): 厳密一致。
- welch_psd <-> scipy.signal.welch(window='hann', nperseg, noverlap=nperseg//2,
  detrend='constant', scaling='density'): 丸めのみ。
- 236.4 Hz 正弦波(格子から0.4ずれ)は矩形窓でしっかり漏れる。
"""

import sys
from pathlib import Path

import numpy as np
import scipy.signal

sys.path.insert(0, str(Path(__file__).parent))
from lib.spectral import dft_naive, fft_radix2, hann, welch_psd

ROOT = Path(__file__).resolve().parents[1]


def rel_l2(a, b):
    return float(np.linalg.norm(a - b) / np.linalg.norm(b))


def main():
    rng = np.random.default_rng(7)
    ok = True

    # 1) FFT 実装同士
    y = rng.normal(size=2048)
    Y = np.fft.fft(y)
    for name, val in (("fft_radix2", fft_radix2(y)), ("dft_naive", dft_naive(y))):
        r = float(np.abs(val - Y).max() / np.abs(Y).max())
        print(f"{name:11s} vs numpy.fft : rel_max = {r:.2e}")
        ok &= r < 1e-10

    # 2) 窓
    w_scipy = scipy.signal.get_window("hann", 2048, fftbins=True)
    r = rel_l2(hann(2048), w_scipy)
    print(f"hann        vs scipy      : rel_l2  = {r:.2e}")
    ok &= r < 1e-14

    # 3) Welch (合成 + IMS 実データ)
    sigs = {"gauss": rng.normal(size=20480)}
    ims = sorted((ROOT / "data/ims/2nd_test").iterdir())[0]
    sigs["IMS Set2 first ch1"] = np.loadtxt(ims)[:, 0]
    for name, x in sigs.items():
        f1, p1 = welch_psd(x, fs=20480.0, nseg=2048)
        f2, p2 = scipy.signal.welch(x, fs=20480.0, window="hann", nperseg=2048,
                                    noverlap=1024, detrend="constant",
                                    scaling="density")
        r = rel_l2(p1, p2)
        print(f"welch_psd   vs scipy      : rel_l2  = {r:.2e}  ({name})")
        ok &= r < 1e-12
        # Parseval 検算: PSD の積分 ~ 分散
        var_psd = float(np.sum(p1) * (f1[1] - f1[0]))
        var_sig = float(np.var(x))
        print(f"  Parseval: integral(PSD)df = {var_psd:.6f}  var(x) = {var_sig:.6f}"
              f"  ratio = {var_psd / var_sig:.4f}")

    # 4) リーケージ定量: 236.4 Hz (格子から 0.4 ずれ) vs 236.0 Hz (格子上)
    fs, N = 20480.0, 20480
    t = np.arange(N) / fs
    print("leakage of a unit sine, rectangular window, N=1s:")
    for f0 in (236.0, 236.4):
        X = np.abs(np.fft.rfft(np.sin(2 * np.pi * f0 * t))) ** 2
        total = X.sum()
        kpk = int(np.argmax(X))
        frac_pk = X[kpk] / total
        frac_3 = X[max(kpk - 1, 0) : kpk + 2].sum() / total
        far = X[kpk + 10] / X[kpk]
        print(f"  f0={f0:6.1f} Hz: peak bin {frac_pk * 100:5.1f}% of energy,"
              f" +/-1 bin {frac_3 * 100:5.1f}%,  10 bins away {10 * np.log10(far):6.1f} dB")

    print("ALL OK" if ok else "SOME CHECK FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
