# bearing-diagnostics

**入口:** [公開レポート](https://yurikada.github.io/bearing-diagnostics/) · [観測設計・データ出典](docs/C0-observation.md) · [独立runへの転移](docs/C7-transfer.md)

成果を読むだけならデータ取得は不要です。処理を再実行するときは、下記の公開データと各ステージのキャッシュを用意します。

転がり軸受の振動データを使った異常検知・診断のケーススタディ。

## 最終ケーススタディ

**[GitHub Pagesで最終ドキュメントを開く](https://yurikada.github.io/bearing-diagnostics/)**

リポジトリをcloneせず、自己完結HTMLの最終レポートをブラウザで確認できる。

公開データ2種を役割分担で使う:

- **CWRU Bearing Data Center**(人工欠陥・正解既知)— 手法の known-answer test
- **NASA IMS Bearing Dataset**(run-to-failure・自然劣化)— 本番の検出問題

コア数学(FFT/PSD、Hilbert変換によるエンベロープ解析、健全性指標、閾値設計)は
`src/lib/` に自前実装し、numpy/scipy と突き合わせる。

## 主要結果(但し書きとセットで引用のこと)

| 結果 | 値 | 条件 |
|---|---|---|
| 検出リードタイム | 3.03日(回顧帯域) / 2.99日(因果帯域) | IMS Set2、事後定義の評価基準まで。警報確定基準 |
| 誤報 | 0件(2層合成・全4ch比較) | Set2開発区間。合成の定義依存性は本文に開示 |
| 周波数指紋 | 236Hz固定 288/288点 | 公開ラベル(外輪)と整合。BPFOを事前アンカーに使用 |
| ライブラリ照合 | rel ≤ 4×10⁻¹³ | 自前実装8種 vs scipy/numpy、規約明示 |
| 転移テスト(独立run) | 部分的成功 | 原理・B4指紋・誤報清潔さは転移/B3主線検出・帰属は非転移([C7](docs/C7-transfer.md)) |

## ステージと文書

| Stage | 内容 | 文書 |
|-------|------|------|
| C0 | 観測設計とデータセット検分(幾何→故障特徴周波数を先に計算) | [C0-observation](docs/C0-observation.md) |
| C1 | 時間領域指標(RMS・尖度・波高率)自前実装、IMS 全期間トレンド | [C1-indicators](docs/C1-indicators.md) |
| C2 | DFT/FFT・窓関数・Welch PSD 自前実装 ↔ numpy/scipy 照合 | [C2-spectral](docs/C2-spectral.md) |
| C3 | エンベロープ解析(Hilbert)自前実装、CWRU 既知欠陥で照合 | [C3-envelope](docs/C3-envelope.md) |
| C4 | 検証済み手法を IMS run-to-failure に適用、劣化開始検出 | [C4-detection](docs/C4-detection.md) |
| C5 | 閾値設計と評価(検出リードタイム vs 誤報率) | [C5-evaluation](docs/C5-evaluation.md) |
| C6 | ケーススタディレポート | [casestudy.html](https://yurikada.github.io/bearing-diagnostics/) |
| C7 | 転移テスト(Set1 内輪+転動体ラベルへの同一パイプライン適用) | [C7-transfer](docs/C7-transfer.md) |

## 再現手順

```
python -m venv .venv
.venv\Scripts\pip install -r requirements.txt
```

1. **データ取得**(`data/` はgitignore。出典URL・サイズ・SHA256は [C0-observation](docs/C0-observation.md) に記録):
   - CWRU: 97 / 105 / 118 / 130.mat → `data/cwru/`
   - NASA IMS: `IMS.zip` を展開 → `data/ims/{1st_test, 2nd_test}/`
2. **実行順**(各ステージのスクリプトは前段のキャッシュ `data/*.npz` を利用):
   `c1_compare → c1_trend → c2_compare → c2_spectra → c3_compare → c3_cwru →
   c4_trend → c5_eval → c7_set1 → c7_set1 --envelope → c7_analyze → c6_build`
3. 図は `docs/figs/`、最終HTMLは `docs/casestudy.html` に再生成される。

## 記録

- **[公開版ケーススタディ](https://yurikada.github.io/bearing-diagnostics/)** — ケーススタディ本体(自己完結・図埋め込み)
- `docs/casestudy.html` — 公開版の生成物
- `scheme.md` — 設計判断ログ(Built / Observed / Decision / Consequence)
- `docs/C*-*.md` — 各ステージの数式・図・コード対応の学習文書
- `data/` — 生データ(gitignore、出典と検分結果は docs 側に記録)

## 評価範囲

- リードタイム3.03日は、IMS Set2で事後定義した評価基準時刻までの回顧評価である。
- 2層判定後の誤報0件は、同じSet2の開発区間で観測された結果である。独立run(Set1)への
  転移テストは docs/C7-transfer.md — 転移した要素と転移しなかった要素を仕分けて記録している。
- lookaheadなしの因果版ではリードタイム2.99日だった。実設備への一般化はC7以降の転移テスト対象とする。

## License

MIT
