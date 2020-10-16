#ifndef DEG
#include "struct.h"
#endif

//有限体の元の逆数
unsigned short
oinv (unsigned short a)
{
  int i;

  for (i = 0; i < N; i++)
    {
      if (gf[mlt (fg[a], i)] == 1)
        return (unsigned short) i;
    }

  printf ("no return \n");
  exit (1);
}


//aに何をかけたらbになるか
unsigned short
equ (unsigned short a, unsigned short b)
{
  int i;


  for (i = 0; i < N; i++)
    {
      if (gf[mlt (fg[a], fg[i])] == b)
        break;
    }
  return i;
}


//OP型からベクトル型への変換
vec
o2v (OP f)
{
  vec a = { 0 };
  int i;

#pragma omp parallel for
  for (i = 0; i < DEG; i++)
    {
      if (f.t[i].a > 0 && f.t[i].n < DEG)
        a.x[f.t[i].n] = f.t[i].a;
    }


  return a;
}


//ベクトル型からOP型への変換
OP
v2o (vec a)
{
  int i, j = 0;
  OP f = { 0 };

  //#pragma omp parallel for
  for (i = 0; i < DEG; i++)
    {
      if (a.x[i] > 0)
        {
          f.t[j].n = i;
          f.t[j++].a = a.x[i];
        }
    }


  return f;
}

//OP型を正規化する
OP
conv (OP f)
{
  vec v = { 0 };
  OP g = { 0 };

  v = o2v (f);
  g = v2o (v);

  return g;
}

//多項式を表示する（OP型）
void
oprintpol (OP f)
{
  int i, n;

  f = conv (f);
  n = distance (f);
  printf ("n=%d\n", n);
  printf ("terms=%d\n", terms (f));
  printf ("deg=%d\n", odeg (f));

  for (i = n; i > -1; i--)
    {
      if (f.t[i].a > 0)
        printf ("%ux^%u+", f.t[i].a, f.t[i].n);
    }
}

void
op_print_raw (const OP f)
{
  puts ("op_print_raw:");
  for (int i = 0; i < DEG; i++)
    {
      if (f.t[i].a > 0)
        printf ("[%d] %ux^%u\n", i, f.t[i].a, f.t[i].n);
    }
}

bool
op_verify (const OP f)
{
  bool end = false;
  unsigned short n_max = 0;
  for (int i = 0; i < DEG; i++)
    {
      if (end && (f.t[i].n != 0 || f.t[i].a != 0))
        {
          op_print_raw (f);
          printf ("found data after end: i=%d\n", i);
          print_trace ();
          fflush (stdout);
          return false;
        }
      if (f.t[i].a == 0)
        {
          end = true;
          continue;
        }
      if (f.t[i].n + 1 <= n_max)
        {
          op_print_raw (f);
          printf ("found invalid order: i=%d\n", i);
          print_trace ();
          fflush (stdout);
          return false;
        }
      n_max = f.t[i].n + 1;
    }
  return true;
}


OP
norm (OP f)
{
  OP h = { 0 };
  int i;


  for (i = 0; i < 512; i++)
    {
      if (f.t[i].a > 0)
        {
          //      h.t[f.t[i].n].n=f.t[i].n;
          h.t[f.t[i].n].a = f.t[i].a;
        }
    }

  // exit(1);

  return h;
}

//20200816:正規化したいところだがうまく行かない
//多項式の足し算
OP
oadd (OP f, OP g)
{
  assert (op_verify (f));
  assert (op_verify (g));

  vec a = { 0 }
  , b = {
    0
  }
  , c = {
    0
  };
  int i, j, k, l = 0;
  OP h = { 0 }, f2 = { 0 }, g2 = { 0 };

  a = o2v (f);
  b = o2v (g);

  if (deg (o2v (f)) >= deg (o2v (g)))
    {
      k = deg (o2v (f)) + 1;
    }
  else
    {
      k = deg (o2v (g)) + 1;
    }

  for (i = 0; i < k; i++)
    {
      c.x[i] = a.x[i] ^ b.x[i];
    }
  h = v2o (c);
  assert (op_verify (h));
  return h;
}


//項の順序を降順に揃える
OP
sort (OP f)
{
  oterm o = { 0 };
  int i, j, k;


  k = terms (f);
  for (i = 0; i < k + 1; ++i)
    {
      for (j = i + 1; j < k + 1; ++j)
        {
          if (f.t[i].n > f.t[j].n)
            {
              o = f.t[i];
              f.t[i] = f.t[j];
              f.t[j] = o;
            }
        }
    }

  return f;
}


//リーディングタームを抽出(OP型）
oterm
oLT (OP f)
{
  int i, k, j;
  oterm s = { 0 };


  k = terms (f);
  s = f.t[0];
  for (i = 0; i < k + 1; i++)
    {
      //printf("a=%d %d\n",f.t[i].a,f.t[i].n);
      if (f.t[i].a > 0)
        {
          printf ("in LT=%d %d\n", s.a, s.n);
          for (j = i; j < k + 1; j++)
            {
              if (s.n < f.t[j].n)
                {
                  s.n = f.t[j].n;
                  s.a = f.t[j].a;
                }

              //  else{
              // t=s;
              // }

            }
        }
    }
  //  exit(1);


  return s;
}

//多項式を足し算する（OP型）
OP
add (OP f, OP g)
{
//  vec a={0},b={0},c={0};
  unsigned long long int i, j, n1 = 0, n2 = 0, m1 = 0, count = 0;
  OP h = { 0 };
  oterm o1 = { 0 }, o2 = {
    0
  };


  n1 = terms (f);
  printf ("n1=%d\n", n1);
  n2 = terms (g);
  printf ("n2=%d\n", n2);
  if (n1 > n2)
    {

    }

  oprintpol (f);
  printf (" fff==============\n");
  oprintpol (g);
  printf (" ggg==============\n");
  o1 = oLT (f);
  o2 = oLT (g);
  printf ("LTadd==%d %d\n", o1.n, o2.n);
  m1 = n1 + n2;
  printf ("m1=%d\n", m1);
  // exit(1);

  for (i = 0; i < n1 + 1; i++)
    {
      for (j = 0; j < n2 + 1; j++)
        {
          if (f.t[i].n == g.t[j].n && g.t[j].a > 0 && f.t[i].a > 0)
            {
              o1 = oLT (f);
              o2 = oLT (g);
              printf ("LT==%d %d\n", o1.n, o2.n);
              printf ("f.n==%d %d %d %d\n", f.t[i].n, g.t[j].n, i, j);
              f.t[i].a = 0;
              g.t[j].a = 0;
            }
        }
    }
  for (i = 0; i < n2 + 1; i++)
    {
      if (g.t[i].a > 0)
        {
          h.t[count++] = g.t[i];
          g.t[i].a = 0;
        }
    }
  for (i = 0; i < n1 + 1; i++)
    {
      if (f.t[i].a > 0)
        {
          h.t[count++] = f.t[i];
          f.t[i].a = 0;
        }

    }

  h = sort (h);
  /*
     for (i=0; i<count; ++i) {
     for (j=i+1; j<count; ++j) {
     if (h.t[i].n > h.t[j].n) {
     oo =  h.t[i];
     h.t[i] = h.t[j];
     h.t[j] = oo;
     }
     }
     }
   */
  h = conv (h);
  if (odeg (h) > 0)
    oprintpol (h);
  printf (" addh==============\n");
  //   exit(1);

  return h;
}


//多項式を項ずつ掛ける
OP
oterml (OP f, oterm t)
{
  assert (op_verify (f));
  int i, k;
  OP h = { 0 };
  vec test;
  unsigned short n;

  k = distance (f);
  for (i = 0; i < k + 1; i++)
    {
      h.t[i].n = f.t[i].n + t.n;
      h.t[i].a = gf[mlt (fg[f.t[i].a], fg[t.a])];
    }

  assert (op_verify (h));
  return h;
}


//多項式の掛け算
OP
omul (OP f, OP g)
{
  assert (op_verify (f));
  assert (op_verify (g));
  int i, count = 0, k;
  oterm t = { 0 };
  OP h = { 0 }, e = {
    0
  }, r = {
    0
  };
  vec c = { 0 };

  if (odeg ((f)) > odeg ((g)))
    {
      k = odeg ((f));
    }
  else
    {
      k = odeg ((g));
    }

  for (i = 0; i < k + 1; i++)
    {
      t = g.t[i];
      e = oterml (f, t);
      h = oadd (h, e);
    }
  assert (op_verify (h));
  return h;
}


//リーディグタームを抽出(default)
oterm
LT (OP f)
{
  int i, k;
  oterm t = { 0 };


  k = deg (o2v (f));
  for (i = 0; i < DEG; i++)
    {
      //printf("a=%d %d\n",f.t[i].a,f.t[i].n);
      if (f.t[i].a > 0)
        {
          t.n = f.t[i].n;
          t.a = f.t[i].a;
        }
    }

  return t;
}


//多項式の最後の項を抽出
oterm
LT2 (OP f)
{
  int i, k;
  oterm t = { 0 };

  t.n = f.t[0].n;
  t.a = f.t[0].a;

  return t;
}


//多項式を単行式で割る
oterm
LTdiv (OP f, oterm t)
{
  oterm tt = { 0 }
  , s = {
    0
  };

  tt = LT (f);
  if (tt.n < t.n)
    {
      s.n = 0;
      s.a = 0;
    }
  else if (tt.n == t.n)
    {
      s.n = 0;
      s.a = equ (t.a, tt.a);
    }
  else if (tt.n > t.n)
    {
      s.n = tt.n - t.n;
      s.a = equ (t.a, tt.a);
      //printf("%u\n",s.a);
    }
  else if (t.n == 0 && t.a > 0)
    {
      s.a = gf[mlt (fg[tt.a], oinv (t.a))];
      s.n = tt.n;
    }

  return s;
}


//モニック多項式にする
OP
coeff (OP f, unsigned short d)
{
  int i, j, k;
  vec a, b;


  f = conv (f);
  k = odeg ((f)) + 1;
  for (i = 0; i < k; i++)
    f.t[i].a = gf[mlt (fg[f.t[i].a], oinv (d))];


  return f;

}


//多項式の剰余を取る
OP
omod (OP f, OP g)
{
  int i = 0, j, n, k;
  OP h = { 0 }, e = {
    0
  };
  oterm a, b = { 0 }, c = {
    0
  };


  n = LT (g).n;

  assert (("baka^\n", LT (f).n != 0));

  /*
     if (LT (f).n == 0)
     {
     printf ("baka^\n");
     //return f;
     exit (1);
     }
   */
  assert (("baka(A)\n", LT (g).n != 0));
  /*
     if (LT (g).n == 0)
     {
     printf ("baka('A`)\n");
     //return g;
     exit (1);
     }
   */
  if (LT (f).n < LT (g).n)
    {
      //    exit(1);
      return f;
    }

  //printf ("in omod\n");
  //exit(1);

  k = LT (g).n;
  b = LT (g);



  //printf ("b=========%dx^%d\n", b.a, b.n);
  //printpol (o2v (g));
  assert (("double baka\n", b.a != 0 && b.n != 0));
  /*
     if (b.a == 0 && b.n == 0)
     {
     printf ("double baka\n");
     exit (1);
     }
   */
  //  printf ("\nin omod g=============%d\n", odeg ((g)));
  while (LT (f).n > 0 && LT (g).n > 0)
    {
      // printf ("in!\n");
      //    exit(1);

      c = LTdiv (f, b);
      //   printf("c========%dx^%d\n",c.a,c.n);
      //    exit(1);

      //printpol (o2v (g));
      //printf ("\ng=================\n");

      h = oterml (g, c);
      //printpol (o2v (h));
      //printf ("\n");
      //printf ("modh===================\n");

      //printpol (o2v (f));
      //printf ("\nmodf===================\n");
      //     exit(1);

      f = oadd (f, h);
      f = conv (f);
      if (odeg ((f)) > 0)
        //printpol (o2v (f));
        //printf ("\nff1=====================\n");
        g = conv (g);
      if (odeg ((f)) == 0 || odeg ((g)) == 0)
        {
          printf ("blake1\n");
          break;
        }

      if (c.n == 0 || b.n == 0)
        break;
    }


  return f;
}


//多項式の商を取る
OP
odiv (OP f, OP g)
{
  assert (op_verify (f));
  assert (op_verify (g));
  int i = 0, j, n, k;
  OP h = { 0 }, e = { 0 }, tt = { 0 };
  oterm a, b = { 0 }, c = { 0 };

  if (LT (f).n == 0 && LT (g).a == 0)
    {
      printf ("baka^\n");
      //return f;
      exit (1);
    }
  if (LT (g).a == 0)
    {
      print_trace ();
      exit (1);
    }
  if (LT (g).n == 0 && LT (g).a > 1)
    return coeff (f, LT (g).a);

  k = odeg (g);
  b = LT (g);
  if (b.a == 1 && b.n == 0)
    return f;
  if (b.a == 0 && b.n == 0)
    {
      printf ("baka in odiv\n");
      exit (1);
    }
  if (odeg ((f)) < odeg ((g)))
    {
      return f;
      //  a=LT(f);
    }

  i = 0;
  while (LT (f).n > 0 && LT (g).n > 0)
    {
      c = LTdiv (f, b);
      assert (c.n < DEG);
      tt.t[i] = c;
      i++;

      h = oterml (g, c);

      f = oadd (f, h);
      if (odeg ((f)) == 0 || odeg ((g)) == 0)
        {
          //printf ("blake2\n");
          break;
        }

      if (c.n == 0)
        break;
    }

  // tt は逆順に入ってるので入れ替える
  OP ret = { 0 };
  int tt_terms = terms (tt);
  for (i = 0; i < tt_terms; i++)
    {
      ret.t[i] = tt.t[tt_terms - i - 1];
    }

  assert (op_verify (ret));
  return ret;
}


//多項式のべき乗
OP
opow (OP f, int n)
{
  int i;
  OP g = { 0 };


  g = f;
  //memcpy(g.t,f.t,sizeof(f.t));

  for (i = 1; i < n; i++)
    g = omul (g, f);


  return g;
}


//多項式のべき乗余
OP
opowmod (OP f, OP mod, int n)
{
  OP g;
  int i;

  g = omod (opow (f, n), mod);


  return g;
}


//多項式の代入値
unsigned short
trace (OP f, unsigned short x)
{
  int i, d;
  unsigned short u = 0;


  d = deg (o2v (f));

  for (i = 0; i < d + 1; i++)
    {
      u ^= gf[mlt (fg[f.t[i].a], mltn (f.t[i].n, fg[x]))];
    }


  return u;
}

// invert of polynomial
OP
inv (OP a, OP n)
{
  OP d = { 0 }, x = {
    0
  }, s = {
    0
  }, q = {
    0
  }, r = {
    0
  }, t = {
    0
  }, u = {
    0
  }, v = {
    0
  }, w = {
    0
  }, tt = {
    0
  }, gcd = {
    0
  };
  oterm b = { 0 };
  vec vv = { 0 }, xx = {
    0
  };


  if (odeg ((a)) > odeg ((n)))
    {
      printf ("baka_i\n");
      exit (1);
    }
  if (LT (a).a == 0)
    {
      printf (" a ga 0\n");
      exit (1);
    }


  tt = n;

  d = n;
  x.t[0].a = 0;
  x.t[0].n = 0;
  s.t[0].a = 1;
  s.t[0].n = 0;
  while (odeg ((a)) > 0)
    {
      if (odeg ((a)) > 0)
        r = omod (d, a);
      if (LT (a).a == 0)
        break;
      if (LT (a).a > 0)
        q = odiv (d, a);

      d = a;
      a = r;
      t = oadd (x, omul (q, s));
      ////printpol (o2v (a));
      //printf ("\nin roop a==================%d\n", odeg ((a)));
      //printf ("\n");

      x = s;
      s = t;
    }
  // exit(1);
  //  if(LT(a).a>0){
  d = a;
  a = r;
  ////printpol (o2v (a));
  //printf ("\nin roop a|==================%d\n", odeg ((a)));
  //printf ("\n");

  x = s;
  s = t;

  ////printpol (o2v (d));
  //printf ("\nout1================\n");
  gcd = d;                      // $\gcd(a, n)$
  ////printpol (o2v (gcd));
  //printf ("\n");
  ////printpol (o2v (n));
  //printf ("\n");
  //printf ("out2===============\n");

  printf ("before odiv\n");
  //w=tt;

  b = LT (w);
  ////printpol (o2v (w));
  //printf ("\nw=======%d %d\n", b.a, b.n);
  //w=tt;
  v = oadd (x, n);
  ////printpol (o2v (v));
  //printf ("\n");
  /*
     if (LT (v).a == 0)
     {
     printf ("v=============0\n");
     }
     printf ("d==============\n");
   */
  //  } //end of a>0
  w = tt;
  ////printpol (o2v (v));
  //printf ("\n");
  //printf ("ss==============\n");
  //       exit(1);
  // if(odeg((w))>0)
  if (LT (v).n > 0 && LT (w).n > 0)
    {
      u = omod (v, w);
    }
  else
    {
      //printpol (o2v (v));
      printf (" v===========\n");
      //printpol (o2v (x));
      printf (" x==0?\n");
      //printpol (o2v (n));
      printf (" n==0?\n");

      exit (1);
    }
  //caution !!
  if (LT (u).a > 0 && LT (d).a > 0)
    {
      u = odiv (u, d);
    }

  if (LT (u).a == 0 || LT (d).a == 0)
    {
      printf ("inv div u or d==0\n");
      // exit(1);
    }
  //u=coeff(u,d.t[0].a);
  ////printpol (o2v (u));
  //printf ("\nu==================\n");
  if (LT (u).a == 0)
    {
      printf ("no return at u==0\n");
      exit (1);
    }


  return u;
}


//error locater for decode
OP
vx (OP f, OP g)
{
  OP h = { 0 }
  , ww = {
    0
  };
  OP v[K] = { 0 }
  , vv = {
    0
  };
  oterm a, b;
  int i, j;

  v[0].t[0].a = 0;
  v[0].t[1].n = 0;
  v[1].t[0].a = 1;
  v[1].t[1].n = 0;

  i = 0;

  while (1)
    {
      if (odeg ((g)) == 0)
        break;
      h = omod (f, g);
      if (LT (g).a == 0)
        break;
      ww = odiv (f, g);
      v[i + 2] = oadd (v[i], omul (ww, v[i + 1]));
      f = g;
      g = h;

      vv = v[i + 2];

      if (odeg ((vv)) == T)
        break;
      i++;
    }

  return vv;
}

//最終の項までの距離
int
distance (OP f)
{
  int i, j = 0;

  for (i = 0; i < DEG; i++)
    {
      if (f.t[i].a > 0)
        j = i;
    }

  return j;
}


//項の数
int
terms (OP f)
{
  int i, count = 0;

  for (i = 0; i < DEG; i++)
    if (f.t[i].a > 0)
      count++;

  return count;
}


//多項式の次数(degのOP型)
int
odeg (OP f)
{
  int i, j = 0, k;


  //k=terms(f);
  for (i = 0; i < 512; i++)
    {
      if (j < f.t[i].n && f.t[i].a > 0)
        j = f.t[i].n;
    }

  return j;
}

//０多項式かどうかのチェック
unsigned char
chk (OP f)
{
  int i, flg = 0;
  vec x = { 0 };

  x = o2v (f);
  for (i = 0; i < DEG; i++)
    {
      if (x.x[i] > 0)
        {
          flg = 1;
          return 1;
        }
    }
  if (flg == 0)
    return 0;

  exit (1);
}


//decode用の多項式の最大公約数
OP
ogcd (OP f, OP g)
{
  OP h;
  //oterm a, b;
  int i = 0;


  //oprintpol((f));
  //oprintpol((g));
  //  exit(1);

  for (i = 0; i < T; i++)
    {
      if (odeg ((g)) == 0)
        break;
      h = omod (f, g);
      if (odeg ((h)) == T - 1)
        {
          //printpol (o2v (h));
          printf (" in ogcd=============\n");
          //wait();
          //break;
          return h;
        }
      f = g;
      g = h;
    }
  // exit(1);


  return h;
}


//拡張ユークリッドアルゴリズム
EX
xgcd (OP f, OP g)
{
  OP h = { 0 }
  , ww = {
    0
  }
  , *v, *u;
  oterm a, b;
  int i = 0, j, k, flg = 0;
  EX e = { 0 };


  v = (OP *) malloc (sizeof (OP) * (DEG));
  u = (OP *) malloc (sizeof (OP) * (DEG));
  memset (v, 0, sizeof (OP) * DEG);
  memset (u, 0, sizeof (OP) * DEG);


  u[0].t[0].a = 1;
  u[0].t[0].n = 0;
  u[1].t[0].a = 0;
  u[1].t[0].n = 0;
  u[2].t[0].a = 1;
  u[2].t[0].n = 0;

  v[0].t[0].a = 0;
  v[0].t[0].n = 0;
  v[1].t[0].a = 1;
  v[1].t[0].n = 0;


  //printpol (o2v (f));
  printf (" f===============\n");
  //printpol (o2v (g));
  printf (" s===============\n");
  // exit(1);


  k = 0;
  i = 0;
  while (1)
    {
      if (LT (g).n == 0)
        {
          printf ("v[%d]=%d skipped deg(g)==0!\n", i, odeg ((v[i])));
          printf (" g========\n");
          exit (1);
        }

      if (LT (g).n > 0)
        h = omod (f, g);

      if (LT (g).a > 0)
        ww = odiv (f, g);

      v[i + 2] = oadd (v[i], omul (ww, v[i + 1]));
      u[i + 2] = oadd (u[i], omul (ww, u[i + 1]));
      //printf ("i+1=%d %d %d g=%d\n", i + 1, odeg ((v[i])), T - 1, odeg ((g)));
      f = g;
      g = h;

      //if(
      if (odeg ((f)) == T - 1 || odeg ((v[i])) == T - 1)
        {
          break;
        }
      i++;
    }

  //v[i]=odiv(v[i],h);
  //u[i]=odiv(u[i],h);
  // h.t[0].a=1;
  //h.t[0].n=0;
  printf ("i=%d\n", i);
  //printpol (o2v (v[i]));
  printf (" v=============\n");
  //printpol (o2v (u[i]));
  printf (" u=============\n");
  //printpol (o2v (f));
  printf (" f=============\n");
  //exit(1);


  e.d = f;
  e.v = v[i];
  e.u = u[i];

  free (v);
  free (u);


  return e;
}


//拡張ユークリッドアルゴリズム(Tで止まらない)
EX
gcd (OP f, OP g)
{
  OP h = { 0 }
  , ww = {
    0
  }
  , v[3] = {
    0
  }
  , u[3] = {
    0
  };
  oterm a, b;
  int i = 0, j, k;
  EX e = { 0 };

  /*
     v = malloc (sizeof (OP) * (DEG));
     u = malloc (sizeof (OP) * (DEG));
     memset (v, 0, sizeof (OP)*DEG);
     memset (u, 0, sizeof (OP)*DEG);
   */

  u[0].t[0].a = 1;
  u[0].t[0].n = 0;
  u[1].t[0].a = 0;
  u[1].t[0].n = 0;
  u[2].t[0].a = 1;
  u[2].t[0].n = 0;

  v[0].t[0].a = 0;
  v[0].t[0].n = 0;
  v[1].t[0].a = 1;
  v[1].t[0].n = 0;


  ////printpol(o2v(f));
  ////printpol(o2v(g));
  //  exit(1);


  k = 0;
  //i=1;
  while (odeg ((g)) > 0)
    {
      if (LT (g).a == 0)
        break;
      if (odeg ((g)) > 0)
        h = omod (f, g);
      if (LT (g).a == 0)
        break;
      if (LT (g).a > 0)
        ww = odiv (f, g);

      v[2] = oadd (v[0], omul (ww, v[1]));
      u[2] = oadd (u[0], omul (ww, u[1]));
      printf ("i+1=%d\n", i + 1);
      f = g;
      g = h;

    }

  //v[i]=odiv(v[i],h);
  //u[i]=odiv(u[i],h);
  // h.t[0].a=1;
  //h.t[0].n=0;
  printf ("i=%d\n", i);
  //printpol (o2v (v[i]));
  printf (" v=============\n");
  //printpol (o2v (u[i]));
  printf (" u=============\n");
  //printpol (o2v (h));
  printf (" h=============\n");
  //exit(1);

  e.d = h;
  e.v = v[1];
  e.u = u[1];

  //free(u);
  //free(v);

  return e;
}

OP
init_pol (OP f)
{
  int i;

  for (i = 0; i < DEG; i++)
    {
      f.t[i].a = 0;
      f.t[i].n = 0;
    }

  return f;
}
