"""周波数領域 (C2): DFT/FFT・窓関数・Welch PSD。

規約 (C2 実装前合意で宣言):
- dft_naive: 定義式そのまま。検証専用 (N <= 4096)。
- fft_radix2: N は 2 のべき乗のみ。
- hann: periodic 版 (scipy の sym=False 相当)。
- welch_psd: 片側 density [unit^2/Hz]。セグメント毎に平均除去、Hann 窓、
  窓補正 U = mean(w^2)、DC/Nyquist 以外を 2 倍、重なり 50%。
  セグメント FFT は自前の fft_radix2 を通す (nseg は 2 のべき乗を要求)。
"""

import numpy as np


def dft_naive(x):
    """定義式 X_k = sum_n x_n exp(-i 2 pi k n / N) の行列積。検証専用。"""
    x = np.asarray(x, dtype=complex)
    N = x.size
    if N > 4096:
        raise ValueError("dft_naive is O(N^2); use N <= 4096")
    n = np.arange(N)
    W = np.exp(-2j * np.pi * np.outer(n, n) / N)
    return W @ x


def fft_radix2(x):
    """radix-2 Danielson-Lanczos。偶数番 E と奇数番 O の N/2 点 DFT から合成する。"""
    x = np.asarray(x, dtype=complex)
    N = x.size
    if N & (N - 1):
        raise ValueError("fft_radix2 requires N to be a power of 2")
    if N == 1:
        return x.copy()
    E = fft_radix2(x[0::2])
    O = fft_radix2(x[1::2])
    tw = np.exp(-2j * np.pi * np.arange(N // 2) / N)
    return np.concatenate([E + tw * O, E - tw * O])


def hann(N):
    """Hann 窓 (periodic 版): w_n = (1 - cos(2 pi n / N)) / 2。"""
    return 0.5 * (1.0 - np.cos(2.0 * np.pi * np.arange(N) / N))


def welch_psd(x, fs, nseg, overlap=0.5):
    """Welch PSD (片側 density)。freqs, psd を返す。

    セグメント毎: 平均除去 -> Hann 窓 -> |FFT|^2 / (fs * nseg * U)。
    U = mean(w^2) は窓がエネルギーを削る分の補正。
    片側化: DC と Nyquist 以外を 2 倍 (nseg 偶数を前提)。
    """
    x = np.asarray(x, dtype=float)
    step = int(nseg * (1.0 - overlap))
    if step < 1:
        raise ValueError("overlap too large")
    w = hann(nseg)
    U = np.mean(w * w)
    n_half = nseg // 2
    acc = np.zeros(n_half + 1)
    count = 0
    for start in range(0, x.size - nseg + 1, step):
        seg = x[start : start + nseg]
        seg = (seg - seg.mean()) * w
        X = fft_radix2(seg)[: n_half + 1]
        p = (np.abs(X) ** 2) / (fs * nseg * U)
        p[1:-1] *= 2.0
        acc += p
        count += 1
    if count == 0:
        raise ValueError("signal shorter than one segment")
    freqs = np.arange(n_half + 1) * (fs / nseg)
    return freqs, acc / count
