# 体験活動プログラム - クラスター計算機の試作と並列計算

@ [奥田・橋本研究室](cc.u-tokyo.ac.jp/support/kosyu/materials/1-1.html)


## 活動記録

| 日 | やったこと |
| --- | --- |
| 3/19 | 環境構築 |
| 3/20 | 環境構築 |
| 3/22 | 環境構築, pi計算(OpenMPIの練習) |
| 3/23 | FrontISTRのインストール |
| 3/26 | CG法の実装 |
| 3/27 | Gauss-Seidel, SOR法の実装 |
| 3/28 | ICCG法実装 |
| 3/29 | ICCG法, SSORなど |
| 3/30 |  |

## 実装した線形ソルバーの種類
### 反復法
#### 定常法
* SOR法 (Successive Over-Relaxation, 逐次加速緩和法)
* Gauss-Seidel法

#### 非定常法
* CG法 (Conjugate Gradient, 共役勾配法)
    * 点ヤコビ前処理
    * [SSOR前処理](https://en.wikipedia.org/wiki/Symmetric_successive_over-relaxation)
    * [ICCG法](http://www.slis.tsukuba.ac.jp/~fujisawa.makoto.fu/cgi-bin/wiki/index.php?%C1%B0%BD%E8%CD%FD%C9%D5%A4%AD%B6%A6%CC%F2%B8%FB%C7%DB%CB%A1) (不完全コレスキー分解を用いる)

### 直接法
* ガウスの消去法

## メモ・教材など
* 並列化の概念についての説明
    * [東京大学情報基盤センター スーパーコンピューティング部門 講習会ビデオ](https://www.cc.u-tokyo.ac.jp/support/kosyu/materials.html)

* 中島先生の数値解析のスライド
    * <http://nkl.cc.u-tokyo.ac.jp/13n/SolverIterative.pdf>
    * <http://nkl.cc.u-tokyo.ac.jp/13n/SolverPrecond.pdf>

* 条件数の計算機
    * <http://comnuan.com/cmnn0100c/>
