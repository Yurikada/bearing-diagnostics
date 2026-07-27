"""C6: 自己完結ケーススタディ HTML の生成。

docs/figs/ の PNG を base64 で埋め込み、docs/casestudy.html を出力する。
構成は B6 で確立した4層: ①成果 ②パイプライン ③学習の過程(外れた予想も記録する)
④転移可能な知見(本人選別: ラベル≠シグネチャ / 検出器の因果性 / 2層防御 /
感度=物理×ベースライン)。
"""

import base64
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FIGS = ROOT / "docs/figs"
OUT = ROOT / "docs/casestudy.html"


def img64(name):
    data = base64.b64encode((FIGS / name).read_bytes()).decode()
    return f"data:image/png;base64,{data}"


HTML = """<!DOCTYPE html>
<html lang="ja">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>転がり軸受の振動診断 — 公開データと自前実装で作る2層防御の異常検知</title>
<style>
  :root {{ --ink:#1a1a1a; --sub:#555; --line:#e0e0e0; --accent:#b34700; --bg:#fdfdfc; }}
  * {{ box-sizing:border-box; }}
  body {{ font-family:"Hiragino Kaku Gothic ProN","Yu Gothic",Meiryo,sans-serif;
         color:var(--ink); background:var(--bg); margin:0; line-height:1.85; }}
  main {{ max-width:900px; margin:0 auto; padding:2rem 1.2rem 4rem; }}
  h1 {{ font-size:1.55rem; line-height:1.5; border-bottom:3px solid var(--accent); padding-bottom:.6rem; }}
  h2 {{ font-size:1.25rem; margin-top:3rem; border-left:5px solid var(--accent); padding-left:.6rem; }}
  h3 {{ font-size:1.05rem; margin-top:2rem; }}
  p, li {{ font-size:.95rem; }}
  .sub {{ color:var(--sub); font-size:.9rem; }}
  .cards {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(250px,1fr)); gap:.8rem; margin:1.2rem 0; }}
  .card {{ border:1px solid var(--line); border-radius:8px; padding: .9rem 1rem; background:#fff; }}
  .card b {{ font-size:1.35rem; color:var(--accent); }}
  .card .s {{ font-size:.82rem; color:var(--sub); }}
  table {{ border-collapse:collapse; width:100%; margin:1rem 0; font-size:.85rem; background:#fff; }}
  th, td {{ border:1px solid var(--line); padding:.45rem .6rem; text-align:left; vertical-align:top; }}
  th {{ background:#f4f1ec; }}
  figure {{ margin:1.4rem 0; }}
  figure img {{ width:100%; border:1px solid var(--line); border-radius:6px; }}
  figcaption {{ font-size:.82rem; color:var(--sub); margin-top:.4rem; }}
  .note {{ background:#fbf6ef; border-left:4px solid var(--accent); padding:.7rem 1rem; font-size:.88rem; margin:1rem 0; }}
  .lost {{ color:#a00; font-weight:bold; }}
  .won {{ color:#06710a; font-weight:bold; }}
  code {{ background:#f2f0ec; padding:.1rem .35rem; border-radius:4px; font-size:.85em; }}
  footer {{ margin-top:4rem; padding-top:1rem; border-top:1px solid var(--line); font-size:.82rem; color:var(--sub); }}
</style>
</head>
<body>
<main>

<h1>転がり軸受の振動診断 —<br>公開データと自前実装で作る「2層防御」の異常検知</h1>
<p class="sub">CWRU Bearing Data(人工欠陥・正解既知)で手法を固定し、NASA IMS Bearing Dataset
(run-to-failure)に適用して検出器を設計・評価したケーススタディ。信号処理のコア数学
(FFT・Welch PSD・Hilbert 変換・時間領域指標)はすべて <code>src/lib/</code> に自前実装し、
scipy/numpy と突き合わせて検証した。2026-07 制作。</p>

<h2>① 成果</h2>
<div class="cards">
  <div class="card"><b>3.03日</b><div class="s">事後定義の評価基準時刻(RMS 5×基準、
  6.74日)までのリードタイム(警報確定 3.708日)。lookahead なしの因果版で 2.99日</div></div>
  <div class="card"><b>誤報 0件(開発区間)</b><div class="s">開発に使用した Set2 の健全軸受
  2本×6.8日で観測された誤報0件(2層の判定適用後)。独立 run では未検証</div></div>
  <div class="card"><b>236 Hz 固定 288/288</b><div class="s">4.5≤day&lt;6.5 の有効288点すべてで
  包絡線ピークが事前計算の BPFO に一致し、公開ラベルの外輪損傷と整合(BPFO を事前アンカーに
  使うため独立診断ではない)</div></div>
  <div class="card"><b>照合 rel ≤ 4×10⁻¹³</b><div class="s">自前実装(DFT/FFT/Hann/Welch/
  Hilbert/RMS/尖度/波高率)を規約明示で scipy/numpy と照合(手法別: FFT系 ~10⁻¹⁶、
  素朴DFT ~10⁻¹³)</div></div>
</div>
<p>副産物として、検証プロセス自体が
<b>ベンチマークのラベルとシグネチャの不一致1件</b>・
<b>公式Readmeとデータ実体の食い違い3件</b>を検出した。</p>
<div class="note"><b>但し書き(この数字を引用する場合は必ず添える)</b>:
リードタイム 3.03日は事後選定帯域(lookahead)による上限値であり、基準時刻(6.74日)は
実故障時刻ではなく事後定義した評価基準(定義変更込みで開示済み)。10分毎1秒の間欠観測のため
真の劣化開始より保守的。run-to-failure 1本(n=1)の結果であり一般化には追試を要する。
警報は「振幅ゲート×周波数ロック」の複合判定で、帰属には全チャンネル比較を伴う。</div>

<h2>② パイプラインと技術的選択</h2>
<table>
<tr><th>Stage</th><th>内容</th><th>核となる技術的選択</th></tr>
<tr><td>C0</td><td>観測設計とデータ検分</td><td>故障特徴周波数(BPFO 236.4Hz 等)を<b>データを見る前に</b>幾何から手計算し予想値とする。検算則 BPFO+BPFI=Z·f<sub>r</sub> ほか</td></tr>
<tr><td>C1</td><td>時間領域指標(RMS・尖度・波高率)</td><td>定義規約(平均除去・Pearson・1/N)を実装前に明示 → scipy と rel=0 一致</td></tr>
<tr><td>C2</td><td>DFT/FFT・リーケージ・Welch PSD</td><td>セグメント長 2048 は Hann 主葉幅×BPFO/BPFI 間隔から導出。fs の矛盾(20k vs 20.48k)を櫛の間隔実測で決着</td></tr>
<tr><td>C3</td><td>エンベロープ解析(Hilbert)</td><td>正解既知の CWRU で known-answer test。帯域選定は PSD 比最大の手続きに固定(自由度の事前封鎖)</td></tr>
<tr><td>C4</td><td>IMS run-to-failure へ適用</td><td>回顧帯域(上限)と因果帯域(因果的な帯域尖度スキャン+中央値5)の2本立てで lookahead を開示</td></tr>
<tr><td>C5</td><td>しきい値設計と評価</td><td>運用曲線(リードタイム vs 誤報)+誤報の時刻内訳+kσ較正。FLI-v2 と max-channel 規則の2層防御を検証</td></tr>
</table>

<figure><img src="{fig_c1}" alt="C1 trend">
<figcaption>C1: Set2 全期間の時間領域指標。教科書の「尖度が先・RMSが後」はこのデータでは
逆順(RMS 3.70日 vs 尖度 4.79日)。B3(緑)は全期間最も衝撃的だが壊れていない。</figcaption></figure>

<figure><img src="{fig_c2}" alt="C2 spectra">
<figcaption>C2: B1 の3時点 Welch PSD。day 5.0 に 3〜6.5kHz の共振帯へ 236Hz 間隔の櫛が成長。
櫛の間隔が fs=20.48kHz 採用を確定させた(すべり 0.17%)。</figcaption></figure>

<figure><img src="{fig_c3}" alt="C3 CWRU">
<figcaption>C3: CWRU known-answer test。内輪は主線 −0.31% + ±f<sub>r</sub> 側帯(荷重域変調の
幾何モデルと整合)。外輪は +0.20%・線/床105。玉はラベルに反しシグネチャ不在(→④主張1)。</figcaption></figure>

<figure><img src="{fig_c4}" alt="C4 trend">
<figcaption>C4: エンベロープ BPFO 線トレンド。3.69日を境に窓内ピークが 236Hz に固定される
(中段)。この変化が検出と診断の両方の根拠になる。</figcaption></figure>

<figure><img src="{fig_c5}" alt="C5 curves">
<figcaption>C5: 運用曲線。しきい値 k はリードタイムと誤報率の交換レートであり、
1点の検出時刻報告では性能を主張できない。</figcaption></figure>

<h2>③ 学習の過程 — 外れた予想も記録する</h2>
<p>各ステージで「測る前に予想を立ててから測る」を徹底した。<b>事前の予想・仮説は
12件で、的中4・外れ8</b>(うち2件はエージェント側の誤りを学習者が計算で正したもの)。
外れた予想は、原因の分析を通じて理解の修正につながった。</p>
<table>
<tr><th>事前の予想(測る前)</th><th>結果</th><th>外れから得られた知見</th></tr>
<tr><td>教科書ストーリー「尖度が先・RMSが後」は綺麗に現れない</td><td><span class="won">的中</span>(ただし想定と別機構)</td><td>RMS が1日先行。感度は物理でなくベースラインの静けさ(CV1.4%)が支配(→④主張4)</td></tr>
<tr><td>RMS 先行の原因は潤滑膜切れ→広帯域上昇(仮説)</td><td><span class="lost">棄却</span></td><td>実体は離散的な櫛の成長。仮説は棄却されたが「広帯域か離散線か」の判別設計を生んだ</td></tr>
<tr><td>236.4Hz は格子に乗りリーケージしない(エージェント出題)</td><td><span class="lost">誤り</span>(本人が指摘)</td><td>0.4bin ずれでピークbinに57%しか残らない。問題文も検証対象</td></tr>
<tr><td>240.0Hz の常在線は電源系 4×60Hz(エージェント仮説)</td><td><span class="lost">誤り</span></td><td>Δf=10Hz の bin 量子化の罠。Δf=1Hz で測ると 237Hz=B3 の BPFO 伝播だった</td></tr>
<tr><td>外輪欠陥に ±f<sub>r</sub> 側帯は立たない</td><td><span class="won">IMS で的中</span> / CWRU では弱く存在</td><td>側帯は二値でなく「荷重に占める回転同期成分の比率」に比例する連続量</td></tr>
<tr><td>エンベロープは RMS より数時間早い</td><td><span class="lost">外れ</span>(確定時刻基準で RMS と同等)</td><td>レイリー床の CV 21% が集約ゲインを相殺。④主張4の再演</td></tr>
<tr><td>B3 の 237Hz 線がベースラインを汚す(「237Hz の壁」)</td><td><span class="lost">不発</span></td><td>加算成分と変調成分は別物 — 帯域通過は加算線を消す。外れたことで両者の区別が明確になった</td></tr>
<tr><td>FLI(周波数ロック)は健全軸受で絶対に鳴らない</td><td><span class="lost">外れ(143件)</span></td><td>偶然ロックは一様床モデルで1試行16%(実装定義: 任意位置3bin幅×3点連続)と高頻度で、試行の反復により頻発(多重比較)。アンカーで検定空間を絞り解消(→④主張3)</td></tr>
<tr><td>CWRU 玉欠陥で 2×BSF 線が立つ</td><td><span class="lost">外れ</span></td><td>採用した帯域探索と指標では 2×BSF 成分を検出できず(→④主張1)。切り分け実験を自設計して判定</td></tr>
<tr><td>内輪欠陥は BPFI ± f<sub>r</sub> 側帯を伴う</td><td><span class="won">的中</span></td><td>荷重域変調の幾何モデルが実データで閉じた</td></tr>
<tr><td>fs=20.48kHz 採用(Readme の矛盾に対する判断)</td><td><span class="won">的中</span></td><td>櫛間隔 236Hz が離散仮説(20.48k: 231.7〜236.4 / 20k: 237.3〜242.1)を判別</td></tr>
<tr><td>故障基準 rms&gt;10×ベースライン</td><td><span class="lost">不成立</span>(最大0.725g&lt;0.77g)</td><td>5×に開示付き改訂。評価パラメータの事後変更は必ず開示する</td></tr>
</table>

<h3>降格した知見(主軸から外したが記録に値するもの)</h3>
<ul>
<li><b>乖離の3スケール分離</b>: ライブラリ照合の乖離は定義差 O(1)/補正差 O(1/N)/丸め
10⁻¹⁶ に分離する。照合閾値は「小さく」ではなく分離したい2スケールの空隙(10⁻¹²)に置く。</li>
<li><b>多重比較と検定空間の直接削減</b>: Bonferroni が閾値を α/m に切り詰めるのに対し、
ドメイン知識(事前計算の BPFO)で無関係な仮説を事前排除する方が検出力を犠牲にしない。</li>
<li><b>メタデータは検分してから使う</b>: IMS 公式 Readme とデータ実体の食い違い3件
(fs 表記、Set3 のフォルダ名、Set3 のファイル数と期間)。公開ベンチマークでも「何が・
どう測られたか」の検分が分析より先。</li>
<li><b>側帯は連続量</b>: 「外輪=側帯なし」という教科書規則は、荷重の回転同期成分比率に
比例する連続量として書き直すと3つのデータ(内輪36%/CWRU外輪14%/IMS外輪≈0%)に整合。</li>
<li><b>測定点の選定が検出可能性を決める</b>: 低周波は構造を遠くまで伝わり(B3 の BPFO 線が
全チャンネルを汚染)、高周波は局所に留まる(B1 の初期欠陥は B3 のセンサでは検出できなかった — 本 run の実測)。センサ配置はアルゴリズム以前の設計変数。</li>
</ul>

<h2>④ 転移可能な知見(主張4本)</h2>

<h3>主張1: ラベルの存在は、シグネチャの励起の十分条件ではない</h3>
<p>CWRU の玉欠陥ファイル(118.mat)は「玉に人工欠陥あり」とラベルされているが、
尖度×帯域のグリッド探索(500Hz〜Nyquist)で全帯域の尖度≈3・BSF 族の線は全帯域で床レベル
であり、採用した手法の範囲では欠陥シグネチャを検出できなかった(物理的不存在までは示さない)。判定表(どの結果ならどの解釈か)を先に書いてから測定し、
文献(Smith &amp; Randall 2015 — 当該記録を診断不成功に分類)と独立に整合した。<b>教師あり学習のベンチマークにおいて、
ラベルを無条件に Ground Truth と扱う前に、物理シグネチャの励起を検証する工程が要る。</b>
これは異常検知に限らず、ラベル付きデータ一般への教訓。</p>

<h3>主張2: 「いつ検出できたか」の主張には、検出器の因果性の開示が要る</h3>
<p>最適な帯域(3〜6.5kHz)は day5 の櫛を見て初めて分かる — それを使って「day 3.7 に検出
できた」と主張するのは、未来の情報から作った部品で過去の検出力を測る「先回り
(lookahead)」であり、オンライン運用の性能証明にならない。
対策は回顧(上限値)と因果(各時点までの情報のみで帯域選定)の<b>2本立て報告</b>で、
両者の差(3.03日 vs 2.99日)自体が「オンライン化のコスト」という情報になる。
<b>時系列上の検出・予測性能を語るすべての仕事で、検出器の部品がいつの情報から
作られたかの開示は、結果の数字と同格に重要である。</b></p>

<h3>主張3: 周波数の同一性と発生源の位置は、直交する2つの検査である</h3>
<p>周波数ロック検査(FLI-v2: 事前計算の BPFO への3点連続ロック×振幅ゲート)は「BPFO 近傍の
周期成分の存在」を支持するが、「どの軸受か」は識別しない — 伝播してきた線も周波数は
同じだから。位置の検査は空間比較(max-channel 規則: 自チャンネルが全チャンネル中最大の
ときのみ発報)が担う。本 run では発生源 B1 の線/床比が常に最大であり、帰属誤りを抑制できた(16件→0件、検出時刻は
不変)。ただしセンサ感度差・伝達経路・同時故障ではこの前提が崩れるため、チャンネル間校正が前提になる。<b>多センサ系の異常検知は「何が起きたか」と「どこで
起きたか」を別の検査に分けると、それぞれ単純な部品で頑健になる。</b></p>

<h3>主張4: 検出感度は、物理的応答とベースライン変動の積で決まる</h3>
<p>教科書の「尖度が先・RMS が後」は本データで逆転した(RMS 3.70日 vs 尖度 4.79日)。
RMS が先行した主因は物理ではなく、ベースライン変動係数 1.4% という変動の小ささである(+5σ が
わずか +7% の変化に相当する)。同じ理屈で、S/N を濃縮するはずのエンベロープ線/床比も
レイリー床の CV 21% が制約となり RMS と同等にとどまった。<b>指標の優劣は固有の性質ではなく、
その現場のベースライン統計との積で決まる。実プラントでは RMS ベースラインの安定は成立しにくい
(負荷・プロセス変動)ため、この積の構造ごと評価してから指標を選ぶ必要がある。</b></p>

<footer>
<p><b>データ</b>: CWRU Bearing Data Center(Case Western Reserve University)/
IMS Bearing Dataset(NASA Prognostics Data Repository, University of Cincinnati IMS Center)。
いずれも公開データ。出典・チェックサム・検分記録は <code>docs/C0-observation.md</code>。</p>
<p><b>実装</b>: Python + numpy/scipy/matplotlib。コア数学は <code>src/lib/</code>
(indicators / spectral / envelope)に自前実装し、確立ライブラリと照合。
各ステージの数式文書は <code>docs/C0〜C5</code>、設計判断ログは <code>scheme.md</code>。</p>
<p><b>制作過程</b>: 本ケーススタディは伴走学習ワークフロー(概念説明→自前実装→コード
読解→ライブラリ照合→図→数式文書→定着チェック→commit)で制作した。予想の設定と
定着チェックはすべて学習者本人が行い、その戦績(的中4・外れ8)を③に記載している。</p>
</footer>

</main>
</body>
</html>
"""


def main():
    html = HTML.format(
        fig_c1=img64("C1-set2-trend.png"),
        fig_c2=img64("C2-b1-spectra.png"),
        fig_c3=img64("C3-cwru-envelope.png"),
        fig_c4=img64("C4-envelope-trend.png"),
        fig_c5=img64("C5-operating-curves.png"),
    )
    OUT.write_text(html, encoding="utf-8")
    print(f"wrote {OUT.name}: {OUT.stat().st_size / 1e6:.2f} MB")


if __name__ == "__main__":
    sys.exit(main() or 0)
