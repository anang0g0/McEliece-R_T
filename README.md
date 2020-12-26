# McEliece-R_T(謎ｗ)

文献1：https://hal.inria.fr/file/index/docid/607772/filename/69.pdf

文献2：https://arxiv.org/pdf/1108.2462.pdf

Reed-Solomonでも暗号に使えるなんて画期的です。

文献を基に実装しましたが、まだ行列Qの構成がわかりません。

また、この方法はまだ十分な時間が経過していないので、攻撃法がこれから見つかるかもしれません。

クローンの方法：
git clone -b master https://github.com/anang0g0/McEliece-R_T

# 20201226

しばしばエラーになるけど並列化して高速化した。

並列化しないと微分の計算でかなりイライラする。
失敗すると不安になるけど、成功確率が50%くらいなのでだめなのかも。

並列化はよくわからない。

# 20201225

有限体上の逆行列を返す関数、invmatを作りました。

# 20201103

実践仕様のパラメータ（[256,128,129]）で動くようにしました。(op2.c)

自動化はまだです。

make dual でビルド。

# 20201102

何だか動いたり動かなかったりで全体を把握しきれてません。
体調不良で頭もおかしいです。

op2.cだけ動作確認しました。
コンパイルオプションmake dualで動くはずです。

# 20201026

元ネタは上の論文にある通りなんですが、ｑ元符号でも暗号に使えるというものです。

特別秘密でもないのですがなぜかいっぱいクローンされていますね。

あと、バグだと思っていたのは仕様のようです。
パターソンアルゴリズムでバイナリゴッパ符号を復号するためにはシンドローム多項式とゴッパ多項式が常に互いに疎でなければならず、
それが保証されるのが既約の場合という感じです。

この問題は小さな例では問題になるのですが、多項式の次数が上がるにつれ、このような特殊な場合が減っていく感じです。

なので既約でない場合でも復号できる関数があるのではないかと思うのですが、このような多項式が存在するかどうかはパターソンアルゴリズムを理解するうえで
ユークリッドアルゴリズムの練習問題にはなるのかなと思います。

# 20201016

根本的な解決にはなっていないけど、エラーが出るような鍵は選ぶときに外すようにした。

# 20201014

気合が足りなくて、バグを追い切れない。
ゆっくりやろう。

# 20201012

久々にバグを見つけてしまいました。

osqrt(多項式の平方根を計算する)で、エラー。
直しますー。
