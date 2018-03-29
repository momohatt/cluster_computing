# 体験活動プログラム - クラスター計算機の試作と並列計算

@ [奥田・橋本研究室](cc.u-tokyo.ac.jp/support/kosyu/materials/1-1.html)


## 活動記録

| 日 | やったこと |
| --- | --- |
| 3/19 | 環境構築 |
| 3/20 | 環境構築 |
| 3/22 | 環境構築, pi計算 |
| 3/23 | FrontISTRのインストール |
| 3/26 | CG法の実装 |
| 3/27 | Gauss-Seidel, SOR法の実装 |
| 3/28 | ICCG法実装 |
| 3/29 | ICCG法実装 |
| 3/30 |  |

## 実装したものとか
* モンテカルロ法を使ったpiの計算
* CG法 (Conjugate Gradient, 共役勾配法)
* SOR法 (Successive Over-Relaxation, 逐次加速緩和法)
* ICCG法 (CG法の前処理として不完全コレスキー分解を用いる)


## メモ
* NIS(Network Information Service)の説明
    * <https://docs.oracle.com/cd/E19455-01/806-2721/6jbtufgne/index.html>
* NFS(Network File System)

* 並列化の概念についての説明
    * [東京大学情報基盤センター スーパーコンピューティング部門 講習会ビデオ](https://www.cc.u-tokyo.ac.jp/support/kosyu/materials.html)

* 中島先生の数値解析のスライド
    * <http://nkl.cc.u-tokyo.ac.jp/13n/SolverIterative.pdf>
    * <http://nkl.cc.u-tokyo.ac.jp/13n/SolverPrecond.pdf>
