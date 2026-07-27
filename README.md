# bearing-diagnostics

転がり軸受の振動データを使った異常検知・診断のケーススタディ。

公開データ2種を役割分担で使う:

- **CWRU Bearing Data Center**(人工欠陥・正解既知)— 手法の known-answer test
- **NASA IMS Bearing Dataset**(run-to-failure・自然劣化)— 本番の検出問題

コア数学(FFT/PSD、Hilbert変換によるエンベロープ解析、健全性指標、閾値設計)は
`src/lib/` に自前実装し、numpy/scipy と突き合わせる。

## ステージ計画

| Stage | 内容 | 状態 |
|-------|------|------|
| C0 | 観測設計とデータセット検分(幾何→故障特徴周波数を先に計算) | ✅ |
| C1 | 時間領域指標(RMS・尖度・波高率)自前実装、IMS 全期間トレンド | ✅ |
| C2 | DFT/FFT・窓関数・Welch PSD 自前実装 ↔ numpy/scipy 照合 | ✅ |
| C3 | エンベロープ解析(Hilbert)自前実装、CWRU 既知欠陥で照合 | ✅ |
| C4 | 検証済み手法を IMS run-to-failure に適用、劣化開始検出 | ✅ |
| C5 | 閾値設計と評価(検出リードタイム vs 誤報率) | ✅ |
| C6 | ケーススタディレポート(docs/casestudy.html) | 次 |

## 記録

- **`docs/casestudy.html`** — ケーススタディ本体(自己完結・図埋め込み)
- `scheme.md` — 設計判断ログ(Built / Observed / Decision / Consequence)
- `docs/C*-*.md` — 各ステージの数式・図・コード対応の学習文書
- `data/` — 生データ(gitignore、出典と検分結果は docs 側に記録)
