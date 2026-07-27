"""エンベロープ解析 (C3): 解析信号・帯域通過・エンベロープスペクトル。

規約 (C3 実装前合意):
- analytic: FFT 法 (負周波数 0、正周波数 x2、DC/Nyquist x1)。scipy.signal.hilbert と照合。
- bandpass: FFT マスク (零位相 brick-wall)。帯域端の Gibbs リンギングは
  1 秒定常スナップショット前提で受容する (宣言済みトレードオフ)。
- envelope_spectrum: bandpass -> analytic -> |z| -> 平均除去 -> 片側振幅スペクトル
  (2/N スケール)。線の読み取りは ±2% 窓・±1bin 積分 (C2 で決めた規則)。
- 帯域選定は正常/欠陥の PSD 比だけで行い、エンベロープ結果を見てから回さない。
"""

import numpy as np


def analytic(x):
    """解析信号 z = x + i x_hat を FFT 法で作る。|z| が瞬時包絡線。"""
    x = np.asarray(x, dtype=float)
    N = x.size
    X = np.fft.fft(x)
    h = np.zeros(N)
    h[0] = 1.0
    if N % 2 == 0:
        h[N // 2] = 1.0
        h[1 : N // 2] = 2.0
    else:
        h[1 : (N + 1) // 2] = 2.0
    return np.fft.ifft(X * h)


def bandpass(x, fs, lo, hi):
    """FFT マスクによる零位相帯域通過。[lo, hi] Hz の外の bin を 0 にする。"""
    x = np.asarray(x, dtype=float)
    N = x.size
    X = np.fft.rfft(x)
    f = np.fft.rfftfreq(N, d=1.0 / fs)
    X[(f < lo) | (f > hi)] = 0.0
    return np.fft.irfft(X, n=N)


def envelope_spectrum(x, fs, lo, hi):
    """帯域 [lo, hi] の包絡線の片側振幅スペクトル (freqs, amp) を返す。"""
    xb = bandpass(x, fs, lo, hi)
    env = np.abs(analytic(xb))
    env = env - env.mean()
    N = env.size
    amp = np.abs(np.fft.rfft(env)) * 2.0 / N
    freqs = np.fft.rfftfreq(N, d=1.0 / fs)
    return freqs, amp
