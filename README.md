# McEliece-R_T(謎ｗ)

https://hal.inria.fr/file/index/docid/607772/filename/69.pdf

https://arxiv.org/pdf/1108.2462.pdf

https://www.zora.uzh.ch/id/eprint/128304/1/Baldi_et_al2014.pdf

https://www.math.unl.edu/~ckelley2/bgklmr2017.pdf

https://user.math.uzh.ch/rosenthal/masterthesis/11720935/Weger_2017.pdf

http://algo.epfl.ch/_media/en/projects/lorenz_thesis.pdf （攻撃法）

https://core.ac.uk/download/pdf/52650821.pdf

https://eprint.iacr.org/2007/382.pdf


# 20210419

このリポジトリを使ってオリジナルのテーマをやってみるつもり。
Alternant Codeってなんだろう？

1変数Goppa符号を、関数の列f_iで作ったらどうなるんだろう？

例えば、

1,1,1,1              |   |f_0   |

a_1,a_2,a_3,a_4      |   |f_1   |

a_1^2,a_2^2,a_3^2    |*  |f_2   |

...                  |   |...   |

a_1^{t-1},..,a_n^{t-1}|  |f_{t-1}|

=Code(L,f_i)

# 20210105

BBCRSはだめらしい。ついでに行列のCSPもだめらしい。何とかならんのか。

符号を使った忘却転送（紛失通信）のデモを作る予定。

# 20201231

なぜ符号の構造を隠すことが難しいのかを考えています。

sidelnikovとshestakovの攻撃の一般化について論文を読む必要があり、
それを自分の考えた暗号化にも適用できるかを研究するつもりです。

例えばランダム行列のもつ群としての性質が符号とは異なるなら、どのような意味で異なるのかという感じです。

もしこの問題で新しい発見があれば、既存の暗号理論が吹っ飛ぶようなことになるかもしれません。

ランダム行列の平均的最小距離を比較して明確な差があればそれも発見になるかもしれないです。
またその違いがある程度大きくならないと識別できないのか、それとも大きさに関係なく識別できるのかにも興味があります。

実験の結果が説明できるような理論にしないとプロにはなれないのですが、そこは素人のすることなので大目に見てくださいｗ

もし自分の方法でランダム行列との識別が困難になるなら、それは成功していると言えます。

# 20201226

書き忘れたけど、mainブランチはほとんどいじってません。

開発用のmasterブランチで並列化や大きなパラメータを扱っているので、最新のが欲しい人は以下のコマンドでmasterブランチをクローンしてください。

git clone -b master https://github.com/anang0g0/McEliece-R_T

q元符号はsidelnikov-shestacov攻撃に弱いから今まで使えなかったけど、これでで安全になったかも。
ｑ現符号なのでバイナリで使っていたpattersonアルゴリズムは使いません。

削除したほうがいいかもしれないですね。

# 20201224

ちなみに今回の暗号はオリジナルじゃないです。

BBCSR方式で、Rの階数zが0、列ベクトルTの重みmが1（置換行列）のシンプルなもの。

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
