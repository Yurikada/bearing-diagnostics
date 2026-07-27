"""時間領域の健全性指標 (C1)。

規約 (C1 実装前合意で宣言):
- 3指標とも平均除去後 (AC成分) に対して定義する。
- 尖度は Pearson 流 (ガウス = 3、-3 しない)、モーメントは 1/N (バイアス補正なし)。
"""

import numpy as np


def rms(x):
    """実効値 sqrt(m2)。m2 = (1/N) sum((x - mean)^2)。"""
    x = np.asarray(x, dtype=float)
    ac = x - x.mean()
    return float(np.sqrt(np.mean(ac * ac)))


def kurtosis(x):
    """Pearson 尖度 m4 / m2^2。m_k = (1/N) sum((x - mean)^k)。ガウスなら 3。"""
    x = np.asarray(x, dtype=float)
    ac = x - x.mean()
    m2 = np.mean(ac * ac)
    if m2 == 0.0:
        raise ValueError("zero-variance signal: kurtosis undefined")
    m4 = np.mean(ac**4)
    return float(m4 / (m2 * m2))


def crest_factor(x):
    """波高率 max|x - mean| / rms(x)。分子・分母とも平均除去後。"""
    x = np.asarray(x, dtype=float)
    r = rms(x)
    if r == 0.0:
        raise ValueError("zero-variance signal: crest factor undefined")
    ac = x - x.mean()
    return float(np.max(np.abs(ac)) / r)
