# McEliece-R_T(謎ｗ)

https://hal.inria.fr/file/index/docid/607772/filename/69.pdf

https://arxiv.org/pdf/1108.2462.pdf

https://www.zora.uzh.ch/id/eprint/128304/1/Baldi_et_al2014.pdf

https://www.math.unl.edu/~ckelley2/bgklmr2017.pdf

Reed-Solomonでも暗号に使えるなんて画期的です。

# 20201103

実践仕様のパラメータで動くようにしました。
自動化はまだです。

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
