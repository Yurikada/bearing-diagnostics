# scheme.md — bearing-diagnostics 判断ログ

形式: Built / Observed / Decision / Consequence(ステージ単位で追記)

## C0 準備(2026-07-26)

- Built: リポジトリ骨格(src/lib, docs/figs, data, tests)、venv。
- Decision: **CWRU(正解既知)で手法を固定してから IMS(run-to-failure)に適用する二段構え。**
  geo-change-detection B4 の「選択を測定可能な時点まで遅らせる」の応用。正解のない
  データに対して手法を選ぶと、選択が評価に混入して検証できなくなる。
- Decision: 生データは `data/` に置き gitignore。出典 URL・ファイル名・サイズ・
  チェックサムは `docs/C0-observation.md` に記録し、再現可能性を担保する。
- Decision: 故障特徴周波数(BPFO/BPFI/BSF/FTF)は**データを見る前に**軸受幾何と
  回転数から計算して宣言する。スペクトルを見てから「合う周波数」を探すと、
  B5 で学んだ「調整した自由度」が増えて一致の情報量が消えるため。
