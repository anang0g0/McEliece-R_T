# McEliece-R_T(謎ｗ)

https://hal.inria.fr/file/index/docid/607772/filename/69.pdf

https://arxiv.org/pdf/1108.2462.pdf

https://www.zora.uzh.ch/id/eprint/128304/1/Baldi_et_al2014.pdf

https://www.math.unl.edu/~ckelley2/bgklmr2017.pdf

Reed-Solomonでも暗号に使えるなんて画期的です。

# 20201226

書き忘れたけど、mainブランチはほとんどいじってません。

開発用のmasterブランチで高速化や大きなパラメータを扱っているので、最新のが欲しい人は以下のコマンドでmasterブランチをクローンしてください。

git clone -b master https://github.com/anang0g0/McEliece-R_T

# 20201224

ちなみに今回の暗号はオリジナルじゃないです。

BBCSR方式で、Tの階数zが0、列ベクトルの重みmが1のシンプルなもの。

一般化するとz=0,m=2の場合、異なる置換P,P'に対して、R=PX+P'Yで表される。

今回の場合、Q=0+Tなので、一応定義に従っている。

T = n × n permutation matrix,

R = n × n rank 1 matrix, R = α^T β,

Q = n × n invertible matrix, Q = R + T,

# 20201223

q元Niederreiter暗号を、現実的なパラメータで実装できた。
セキュリティパラメータは大体2^80くらいである。

# 20201221

今まで頻繁に止まっていたバグを、GCDの関数を変えることで解決。

油断は禁物。
こんな簡単なことで本当に強度が上がるんだろうか？
ちょっと懐疑的。

# 20201026

コード管理のために、品質重視のブランチと新規性優先のブランチに分けました。

前者のブランチはmain、後者はmasterです。
stableとdevの違いのようなものです。

# 20201025

半信半疑だったけどちゃんと動いてる。奇跡だｗ

# 20201016

根本的な解決にはなっていないけど、エラーが出るような鍵は選ぶときに外すようにした。

新しい開発メンバーにkuboonさんが加わりました。

# 20201014

気合が足りなくて、バグを追い切れない。
ゆっくりやろう。

# 20201012

久々にバグを見つけてしまいました。

osqrt(多項式の平方根を計算する)で、エラー。
直しますー。
