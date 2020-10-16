//date 20200331:xgcdの終了条件が２つになってしまう。ogcdとxgcdで使い分ける。
//date : 20200326 鍵生成が det と deta で高い確率で一致する。detaは並列処理。
//date 20200229 : pattarson algorithm implementation ver 1.0
// xgcd & osqrtを追加した
//date      :  20160310,20191218,20191220,20191221,20191223,20191224,20191225,20191229,20191230
//auther    : the queer who thinking about cryptographic future
//code name :  一変数多項式演算ライブラリのつもり
//status    : now in debugging (ver 0.8)
// 0ベクトルが出ないように生成多項式のトレースチェックを入れた。
//date      :  20160310
//auther    : the queer who thinking about cryptographic future
//code name : OVP - One Variable Polynomial library with OpenMP friendly
//status    : now in debugging (ver 0.1)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include <stdbool.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <assert.h>
#include <execinfo.h>


//#include "omp.h" //clang-12
#include <omp.h>                //clang-10

#include "debug.c"
//#include "8192.h"
#include "global.h"
#include "struct.h"

#include "chash.c"
#include "lu.c"
#include "sha3.c"
#include "inv_mat.c"
#include "op.c"

extern unsigned long xor128 (void);
extern int mlt (int x, int y);
extern int mltn (int n, int a);
extern void makeS ();

//#pragma omp threadprivate(mat)
//シンドロームのコピー
unsigned short sy[K] = { 0 };

//Goppa多項式
static unsigned short g[K + 1] = { 0 }; //{1,0,6,0,0,0,0,8,1};
unsigned short zz[N] = { 0 };


/*
static unsigned short g[K+1]={1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,
		       0,0,0,0,0,0,0,0,0,0,
		       		       0,0,0,0,0,0,0,0,0,0,
		       		       0,0,0,0,0,0,0,0,0,0,
		       		       0,0,0,0,0,0,0,0,0,0,
		       0,1};
*/

//ランダム多項式の生成
static void
ginit (void)
{
  int j, count = 0;
  unsigned short gg[K + 1] = { 0 };
  //unsigned short g[K+1]={0};

  printf ("in ginit\n");


  g[K] = 1;                     //xor128();
  g[0] = rand () % N;
  while (count < ((K / 2) - 1))
    {
      printf ("@\n");
      j = rand () % (K - 1);
      if (j < K && j > 0 && g[j] == 0)
        {
          g[j] = rand () % N;
          count++;
        }
    }


  for (j = 0; j < K + 1; j++)
    gg[j] = g[K - j];

  memcpy (g, gg, sizeof (g));


}



//ランダム置換の生成（Niederreoter 暗号における置換）
void
random_permutation (unsigned short *a)
{
  int i, j, x;

  for (i = 0; i < N; i++)
    {
      a[i] = i;
    }
  for (i = 0; i < N - 2; i++)
    {
      j = (rand () % (N - 1 - i)) + i + 1;

      x = a[j];
      a[j] = a[i];
      a[i] = x;
    }
  if (a[N - 1] == N - 1)
    {
      a[N - 1] = a[N - 2];
      a[N - 2] = N - 1;
    }


}


//配列から置換行列への変換
void
P2Mat (unsigned short P[N])
{
  int i;

  for (i = 0; i < N; i++)
    A[i][P[i]] = 1;
}


unsigned short
b2B (unsigned short b[E])
{
  int i;
  unsigned short a = 0;

  for (i = E - 1; i > -1; i--)
    a ^= (b[E - i - 1] << i);

  return a;
}


//多項式の次数(default)
int
deg (vec a)
{
  int i, n = 0;

  //#pragma omp parallel for
  for (i = 0; i < DEG; i++)
    {
      if (a.x[i] > 0)
        n = i;
    }


  return n;
}



//整数からベクトル型への変換
vec
i2v (unsigned int n)
{
  vec v;
  int i;


  for (i = 0; i < 32; i++)
    {
      v.x[i] = n % 2;
      n = (n >> 1);
    }


  return v;
}


//ベクトル型から整数への変換
unsigned int
v2i (vec v)
{
  unsigned int d = 0, i;


  for (i = 0; i < 32; i++)
    {
      d = (d << 1);
      d ^= v.x[i];
    }


  return d;
}


//配列からベクトル表現の多項式へ変換する
vec
Setvec (int n)
{
  int i;
  vec v = { 0 };


  for (i = 0; i < n; i++)
    {
      v.x[n - 1 - i] = c[i];
    }


  return v;
}


void
printvec (vec v)
{
  int i, j;

  for (i = 0; i < deg (v) + 1; i++)
    {
      if (v.x[i] > 0)
        printf ("g[%d]=%d;\n", i, v.x[i]);
    }

}


//多項式を表示する(default)
void
printpol (vec a)
{
  int i, n;

  n = deg (a);

  //printf ("baka\n");
  assert (("baka\n", n >= 0));



  for (i = n; i > -1; i--)
    {
      if (a.x[i] > 0)
        {
          printf ("%u", a.x[i]);
          //if (i > 0)
          printf ("x^%d", i);
          //if (i > 0)
          printf ("+");
        }
    }
  //  printf("\n");

  return;
}



vec
vadd (vec a, vec b)
{
  vec c = { 0 };
  int i, k;


  if (deg (a) >= deg (b))
    {
      k = deg (a) + 1;
    }
  else
    {

      k = deg (b) + 1;

    }
  for (i = 0; i < k; i++)
    {
      c.x[i] = a.x[i] ^ b.x[i];
    }
  // c.x[i] = a.x[i] ^ b.x[i];
  //  h = v2o (c);

  return c;
}



//停止コマンド
void
wait (void)
{
  int n;                        // 読み込む変数はローカルに取るべし
  printf (" (enter number and hit return) ");   // 何か表示させたほうが良いだろう
  fflush (stdout);              // just in case
  scanf ("%d", &n);             // fgets(line, LINESIZE, stdin); という手も
}



//整数のべき乗
unsigned int
ipow (unsigned int q, unsigned int u)
{
  unsigned int i, m = 1;

  for (i = 0; i < u; i++)
    m *= q;

  printf ("in ipow====%d\n", m);

  return m;
}



//多項式の形式的微分
OP
bibun (vec a)
{
  OP w[T * 2] = { 0 };
  OP l = { 0 }
  , t = {
    0
  };
  int i, j, n, id;
  vec tmp = { 0 };



  n = deg (a);
  printf ("n=%d\n", n);
  if (n == 0)
    {
      printf ("baka8\n");
      //  exit(1);
    }

  //
  for (i = 0; i < T; i++)
    {
      w[i].t[0].a = a.x[i];
      w[i].t[0].n = 0;
      w[i].t[1].a = 1;
      w[i].t[1].n = 1;
      ////printpol(o2v(w[i]));
    }
  //  exit(1);

  tmp.x[0] = 1;
  //

  //#pragma omp parallel for private(i,j)
  for (i = 0; i < T; i++)
    {
      t = v2o (tmp);
      //
      for (j = 0; j < T; j++)
        {
          // #pragma omp critical
          if (i != j)
            {
              t = omul (t, w[j]);
            }
        }

      ////printpol(o2v(t));

      if (odeg ((t)) == 0)
        {
          printf ("baka9\n");
          // exit(1);
        }
      l = oadd (l, t);
    }


  return l;
}


//chen探索
vec
chen (OP f)
{
  vec e = { 0 };
  int i, count = 0, n, x = 0;
  unsigned short z;


  n = odeg ((f));
//exit(1);
//#pragma omp parallel for private(i)
  for (x = 0; x < N; x++)
    {
      z = 0;
#pragma omp parallel for reduction (^:z)
      for (i = 0; i < n + 1; i++)
        {
          if (f.t[i].a > 0)
            z ^= gf[mlt (mltn (f.t[i].n, fg[x]), fg[f.t[i].a])];
        }
      if (z == 0)
        {
          e.x[count] = x;
          count++;
        }
      //    printf("%d\n",x);
    }


  return e;
}


//ユークリッドアルゴリズムによる復号関数
OP
decode (OP f, OP s)
{
  int i, j, k, count = 0;
  OP r = { 0 }, w = {
    0
  }, e = {
    0
  }, l = {
    0
  };
  oterm t1, t2, d1, a, b;
  vec x = { 0 };
  unsigned short d = 0;
  OP h = { 0 };
  EX hh = { 0 };


  printf ("in decode\n");
  //printpol (o2v (s));
  printf ("\nsyn===========\n");
  r = vx (f, s);
  //h=ogcd(f,s);

  if (odeg ((r)) == 0)
    {
      printf ("baka12\n");
      exit (1);
    }
  k = 0;
  // exit(1);
  x = chen (r);
  // exit(1);



  for (i = 0; i < T; i++)
    {
      // printf ("x[%d]=1\n", x.x[i]);
      if (x.x[i] == 0)
        k++;
      if (k > 1)
        {
          printf ("baka0\n");
          printvec (o2v (f));
          for (i = 0; i < N; i++)
            printf ("%d,", zz[i]);
          exit (1);
          //return f;
        }
    }
  //exit(1);

  //  printf("\n");

  printf ("あっ、でる！\n");
  //  exit(1);

  if (odeg ((r)) < T)
    {
      //printpol (o2v (r));
      printf ("baka5 deg(r)<T\n");
      exit (1);
    }


  w = bibun (x);
  //exit(1);
  //  w=oterml(w,d1);
  //printpol (o2v (w));
  printf ("@@@@@@@@@\n");
  //exit(1);

  hh = xgcd (f, s);
  //printpol (o2v (hh.d));
  //wait();


  //  exit(1);
  t1 = LT (r);

  t2.a = t1.a;
  t2.n = 0;

  if (odeg ((w)) == 0)
    {
      //printpol (o2v (w));
    }
  l = oterml (w, t2);

  j = deg (x) + 1;
  printf ("%d\n", j);

  //    exit(1);

  for (i = 0; i < j; i++)
    {
      if (x.x[i] > 0)
        {
          e.t[i].a =
            gf[mlt (fg[trace (hh.d, x.x[i])], oinv (trace (l, x.x[i])))];
          //e.t[i].a = gf[mlt (fg[trace (h, x.x[i])], oinv (trace (l, x.x[i])))];
          e.t[i].n = x.x[i];
        }
    }
  //printpol (o2v (f));
  printf (" f============\n");
  //printpol (o2v (l));
  printf (" l============\n");
  //  exit(1);


  for (i = 0; i < T; i++)
    if (gf[trace (h, x.x[i])] == 0)
      printf ("h=0");
  //printf("\n");
  for (i = 0; i < T; i++)
    if (gf[oinv (trace (l, x.x[i]))] == 0)
      printf ("l=0\n");
  //  printf("\n");


  return e;
}



//配列の値を係数として多項式に設定する
OP
setpol (unsigned short f[], int n)
{
  OP g;
  vec a;
  int i;


  memset (c, 0, sizeof (c));
  memcpy (c, f, 2 * n);
  a = Setvec (n);

  g = v2o (a);


  return g;
}



unsigned short tr[N] = { 0 };
unsigned short ta[N] = { 0 };

void
det2 (int i, unsigned short g[])
{
  OP f[16] = { 0 }, h[16] = {
    0
  }, w, u[16] = {
    0
  };
  unsigned short cc[K + 1] = { 0 }, d[2] = {
    0
  };
  int j, a, b, k, t1, l = 0, flg = 0, id;
  oterm t[16] = { 0 };
  vec e[16] = { 0 };
  OP ww[16] = { 0 };


  memcpy (cc, g, sizeof (cc));

  k = cc[K];
  w = setpol (g, K + 1);

  //  omp_set_num_threads(8);
  id = omp_get_thread_num ();

  // h[id] = x+i
  if (i == 0)
    {
      h[id].t[0].a = 1;
      h[id].t[0].n = 1;
    }
  else
    {
      h[id].t[0].a = i;
      h[id].t[1].a = 1;
      h[id].t[1].n = 1;
    }

  t[id].n = 0;

  f[id] = setpol (cc, K + 1);

  cc[K] = k ^ ta[i];
  //tr[i];
  f[id] = setpol (cc, K + 1);

  //f.t[0].a=k^ta[i]; //cc[K];

  ww[id] = odiv (f[id], h[id]);

  //b = oinv (a);
  t[id].a = gf[tr[i]];
  u[id] = oterml (ww[id], t[id]);
  e[id] = o2v (u[id]);

  memcpy (mat[i], e[id].x, sizeof (mat[i]));


}





void
f1 (unsigned short g[])
{
  int i;

  for (i = 0; i < 836; i++)
    {
      det2 (i, g);
    }


}

void
f2 (unsigned short g[])
{
  int i;

  for (i = 836; i < 1672; i++)
    {
      det2 (i, g);
    }


}

void
f3 (unsigned short g[])
{
  int i;

  for (i = 1672; i < 2508; i++)
    {
      det2 (i, g);
    }


}

void
f4 (unsigned short g[])
{
  int i;

  for (i = 2508; i < 3344; i++)
    {
      det2 (i, g);
    }


}

void
f5 (unsigned short g[])
{
  int i;

  for (i = 3344; i < 4180; i++)
    {
      det2 (i, g);
    }


}

void
f6 (unsigned short g[])
{
  int i;

  for (i = 4180; i < 5016; i++)
    {
      det2 (i, g);
    }


}

void
f7 (unsigned short g[])
{
  int i;

  for (i = 5016; i < 5852; i++)
    {
      det2 (i, g);
    }


}

void
f8 (unsigned short g[])
{
  int i;

  for (i = 5852; i < 6688; i++)
    {
      det2 (i, g);
    }


}


int
detb (unsigned short g[])
{
  int i, j, flg = 0;
  /*
     for(i=0;i<N;i++){
     for(j=0;j<K;j++)
     mat[i][j]=0;
     }
   */
#pragma omp parallel num_threads(8)
  {
#pragma omp sections
    {
      //if(omp_get_thread_num() == 0){
#pragma omp section
      f1 (g);

      //if(omp_get_thread_num() == 1){
#pragma omp section
      f2 (g);
      //}
      //if(omp_get_thread_num() == 2){
#pragma omp section
      f3 (g);
      //}
      //if(omp_get_thread_num() == 3){
#pragma omp section
      f4 (g);
      //}
      //if(omp_get_thread_num() == 4){
#pragma omp section
      f5 (g);
      //}
      //if(omp_get_thread_num() == 5){
#pragma omp section
      f6 (g);
      //}
      //if(omp_get_thread_num() == 6){
#pragma omp section
      f7 (g);
      //}
      //if(omp_get_thread_num() == 7){
#pragma omp section
      f8 (g);
      //}
    }
  }
  printf ("enf of detb\n");
  for (j = 0; j < N; j++)
    {
      flg = 0;
      for (i = 0; i < K; i++)
        {
          //printf("%d,",mat[i][j]);
          if (mat[j][i] > 0)
            flg = 1;
          //      printf("\n");
        }
      if (flg == 0)
        {
          printf ("0 is %d\n", j);
          //exit(1);
          return -1;
        }
    }

  return 0;
}



//パリティチェック行列を生成する
int
deta (unsigned short g[])
{
  int i, j, a, b, k, t1, l = 0, flg = 0, id;


  //
  //
#pragma omp parallel num_threads(8)
  {
#pragma omp for schedule(static)
    for (i = 0; i < D; i++)
      {
        det2 (i, g);
      }
  }
  for (j = 0; j < D; j++)
    {
      flg = 0;
      for (i = 0; i < K; i++)
        {
          printf ("%d,", mat[i][j]);
          if (mat[j][i] > 0)
            flg = 1;
        }
      printf ("\n");
      if (flg == 0)
        {
          printf ("0 is %d\n", j);
          //exit(1);
          return -1;
        }
    }
  printf ("end2\n");
  // exit(1);
  return 0;
}



unsigned short *base;

//パリティチェック行列を生成する
void
det (unsigned short g[])
{
  OP f, h = { 0 }, w, u;
  unsigned short cc[K + 1] = { 0 }, d[2] = {
    0
  }, pp[20][K] = {
    0
  };
  int i, j, a, b, k, t1, l = 0, flg = 0, count = 0;
  oterm t = { 0 };
  vec e;


  memcpy (cc, g, sizeof (cc));
  /*
     for (i = 0; i < K + 1; i++)
     {
     cc[i] = g[i];
     printf ("%d,", g[i]);
     }
   */
  //printf ("\n");
  //  exit(1);
  //    cc[i]=g[i];
  k = cc[K];
  w = setpol (g, K + 1);

  OP ww = { 0 };

  h.t[0].n = 0;
  h.t[1].a = 1;
  h.t[1].n = 1;
  t.n = 0;
  t1 = 2 * T;
  // #pragma omp parallel for
  /*
     unsigned short tr[N];
     unsigned short ta[N];
     for(i=0;i<N;i++){
     ta[i] = trace (w, i);
     if(ta[i]==0){
     printf("%d %d\n",i,ta[i]);
     exit(1);
     }   
     tr[i] = oinv (ta[i]);    
     }
   */

  //
  f = setpol (cc, K + 1);

  for (i = 0; i < D; i++)
    {

      a = trace (w, i);
      // cc[K] = k;

      cc[K] = k ^ a;
      //tr[i];
      f = setpol (cc, K + 1);

      //f.t[0].a=k^ta[i]; //cc[K];
      h.t[0].a = i;

      ww = odiv (f, h);

      b = oinv (a);
      t.a = gf[b];


      u = oterml (ww, t);
      e = o2v (u);

      // #pragma omp parallel for
      //for (j = 0; j < K; j++)
      //mat[i][j]= e.x[j];
      memcpy (mat[i], e.x, sizeof (mat[i]));
    }


  for (j = 0; j < D; j++)
    {
      flg = 0;
      for (i = 0; i < K; i++)
        {
          printf ("%d,", mat[i][j]);
          if (mat[j][i] > 0)
            flg = 1;
        }
      printf ("\n");
      if (flg == 0)
        {
          printf ("0 is %d\n", j);
          //exit(1);
        }
    }
  //exit(1);

}


void
detc (unsigned short g[])
{



}


//バイナリ型パリティチェック行列を生成する
void
bdet ()
{
  int i, j, k, l;
  unsigned char dd[E * K] = { 0 };
  FILE *ff;


  ff = fopen ("Hb.key", "wb");


  for (i = 0; i < N; i++)
    {
      for (j = 0; j < K; j++)
        {
          l = mat[i][j];
#pragma omp parallel for
          for (k = 0; k < E; k++)
            {
              BH[j * E + k][i] = l % 2;
              l = (l >> 1);
            }
        }
    }
  for (i = 0; i < N; i++)
    {
#pragma omp parallel for
      for (j = 0; j < E * K; j++)
        {
          //  printf("%d,",BH[j][i]);
          dd[j] = BH[j][i];
        }
      fwrite (dd, 1, E * K, ff);
      //printf("\n");
    }
  fclose (ff);

}




//Niederreiter暗号の公開鍵を作る
void
pubkeygen ()
{
  int i, j, k, l;
  FILE *fp;
  unsigned char dd[E * K] = { 0 };


  fp = fopen ("pub.key", "wb");
#pragma omp parallel for private(j,k)
  for (i = 0; i < E * K; i++)
    {
      for (j = 0; j < N; j++)
        {
          //tmp[i][j]=0;
          l = 0;
          //#pragma omp parallel for reduction (^:l)
          for (k = 0; k < E * K; k++)
            {
              tmp[i][j] ^= cl[i][k] & BH[k][j];
              //l^= cl[i][k] & BH[k][j];
            }
          //tmp[i][j]=l;
        }
      //}
    }
  P2Mat (P);

#pragma omp parallel for
  for (i = 0; i < E * K; i++)
    {
      //  for(j=0;j<N;j++){
#pragma omp parallel for
      for (k = 0; k < N; k++)
        pub[i][k] = tmp[i][P[k]];       //&A[k][j];
      //    }
    }

#pragma omp parallel for private(j)
  for (i = 0; i < N; i++)
    {
      for (j = 0; j < E * K; j++)
        {
          dd[j] = pub[j][i];
        }
      fwrite (dd, 1, E * K, fp);
    }
  fclose (fp);

}


//秘密置換を生成する
void
Pgen ()
{
  unsigned int i, j;
  FILE *fp;


  fp = fopen ("P.key", "wb");
  random_permutation (P);
  for (i = 0; i < N; i++)
    inv_P[P[i]] = i;
  fwrite (P, 2, N, fp);
  fclose (fp);

  //for (i = 0; i < N; i++)
  //printf ("%d,", inv_P[i]);
  //printf ("\n");


  fp = fopen ("inv_P.key", "wb");
  fwrite (inv_P, 2, N, fp);
  fclose (fp);

}


//鍵生成
void
key2 (unsigned short g[])
{
  FILE *fp;
  unsigned short dd[K] = { 0 };
  int i, j, k;

  printf ("鍵を生成中です。４分程かかります。\n");
  fp = fopen ("H.key", "wb");
  i = 0;
  do
    {
      i = deta (g);
    }
  while (i == -1);

  //exit(1);
  for (i = 0; i < N; i++)
    {
      for (j = 0; j < K; j++)
        dd[j] = mat[i][j];
      fwrite (dd, 2, K, fp);

    }
  fclose (fp);
  fp = fopen ("sk.key", "wb");
  fwrite (g, 2, K + 1, fp);
  fclose (fp);

}


//すべての鍵を生成する
void
keygen (unsigned short *g)
{
  int i;
  FILE *fp;


  key2 (g);
  printf ("end of ky2\n");
  makeS ();
  printf ("end of S\n");
  bdet ();
  printf ("end of bdet\n");
  Pgen ();
  printf ("end of Pgen\n");
  pubkeygen ();

}


//ハッシュ１６進表示
static void
byte_to_hex (uint8_t b, char s[23])
{
  unsigned i = 1;
  s[0] = s[1] = '0';
  s[2] = '\0';
  while (b)
    {
      unsigned t = b & 0x0f;
      if (t < 10)
        {
          s[i] = '0' + t;
        }
      else
        {
          s[i] = 'a' + t - 10;
        }
      i--;
      b >>= 4;
    }
}




//有限体の元の平方を計算する
int
isqrt (unsigned short u)
{
  int i, j, k;


  for (i = 0; i < N; i++)
    {
      if (gf[mlt (i, i)] == u)
        return i;
    }

  printf ("来ちゃいけないところに来ました\n");
  exit (1);
}


//多項式の平方を計算する
OP
osqrt (OP f, OP w)
{
  int i, j, k, jj, n;
  OP even = { 0 }, odd = {
    0
  }, h = {
    0
  }, r = {
    0
  }, ww = {
    0
  }, s = {
    0
  }, tmp = {
    0
  }, t = {
    0
  };
  oterm o = { 0 };
  vec v = { 0 };

  j = 0;
  jj = 0;
  k = deg (o2v (f));
  for (i = 0; i < k + 1; i++)
    {
      if (f.t[i].n % 2 == 0 && f.t[i].a > 0)
        {
          even.t[j].n = f.t[i].n / 2;
          even.t[j++].a = gf[isqrt (f.t[i].a)];
          printf ("a=%d %d\n", f.t[i].a, i);
        }
      if (f.t[i].n % 2 == 1 && f.t[i].a > 0)
        {
          odd.t[jj].n = (f.t[i].n - 1) / 2;
          odd.t[jj++].a = gf[isqrt (f.t[i].a)];
          printf (" odd %d\n", i);
        }
    }


  k = deg (o2v (w));
  //printf ("%d\n", k);
  //exit(1);
  j = 0;
  jj = 0;
  for (i = 0; i < k + 1; i++)
    {
      if (w.t[i].n % 2 == 0 && w.t[i].a > 0)
        {
          h.t[j].a = gf[isqrt (w.t[i].a)];
          h.t[j++].n = w.t[i].n / 2;
          printf ("h==%d %d\n", (w.t[i].n / 2), i);
        }
      if (w.t[i].n % 2 == 1 && w.t[i].a > 0)
        {
          r.t[jj].a = gf[isqrt (w.t[i].a)];
          r.t[jj++].n = (w.t[i].n - 1) / 2;
          printf ("r=====%d %d\n", (w.t[i].n - 1) / 2, i);
        }
    }
  printpol (o2v (r));
  printf (" sqrt(g1)=======\n");

  //  exit(1);
  if (LT (r).n > 0)
    {
      s = inv (r, w);
    }
  else if (LT (r).n == 0)
    {
      printpol (o2v (r));
      printf (" deg(r)======0!\n");
      printpol (o2v (w));
      printf (" goppa======0\n");
      printpol (o2v (f));
      printf (" syn======0\n");
      if (LT (r).a > 0)
        {
          s.t[0].a = gf[isqrt (LT (r).a)];
          s.t[0].n = 0;
          return omod (oadd (even, omul (coeff (h, s.t[0].a), odd)), w);
        }
      return even;
      //s = inv (r, w);
      //wait ();
      //exit (1);
    }
  if (deg (o2v (s)) > 0)
    tmp = omod (omul (s, r), w);
  if (odeg ((tmp)) > 0)
    {
      //printpol (o2v (tmp));
      printf (" r is not inv==========\n");
      wait ();
      exit (1);
    }
  if (LT (h).n > 0)
    {
      ww = omod (omul (h, s), w);
    }
  if (LT (h).n == 0)
    {
      printpol (o2v (h));
      printf (" h=========0\n");
      //ww = omod (omul (h, s), w);
      //wait();
      exit (1);
    }

  if (LT (ww).n == 0 && LT (ww).a == 0)
    {
      //printpol(o2v(s));
      printf (" s===========\n");
      //printpol(o2v(w));
      printf (" w==============\n");
      //printpol(o2v(r));
      printf (" r===========\n");
      //printpol(o2v(h));
      printf (" h============\n");
      //printpol (o2v (ww));
      printf (" ww==============\n");
      printf (" wwが０になりました。error\n");
      wait ();
      //return ww;;
      exit (1);
    }

  tmp = omod (omul (ww, ww), w);
  if (LT (tmp).n == 1)
    {
      printpol (o2v (ww));
      printf (" ww succsess!===========\n");
    }
  else
    {
      //printpol (o2v (tmp));
      printf (" mod w^2==========\n");
      //printpol (o2v (ww));
      printf (" ww^2 failed!========\n");
      printpol (o2v (s));
      printf (" g1^-1==============\n");
      //printpol (o2v (w));
      printf (" w==============\n");
      //printpol (o2v (h));
      printf (" g0===========\n");
      //printpol (o2v (r));
      printf (" r===========\n");
      printf ("この鍵では逆元が計算できません。error");
      wait ();
      //return ww;
      exit (1);
    }


  //    exit(1);
  printpol (o2v (s));
  printf (" g1^-1=========\n");
  printpol (o2v (h));
  printf (" g0=========\n");
  //exit(1);
  printpol (o2v (ww));
  printf (" ww==========\n");
  //  exit(1);
  h = ww;
  if (odeg (omod (omul (h, ww), w)) == 1)
    {
      ww = h;
      h = omod (oadd (even, omul (ww, odd)), w);
      return h;
    }
  else if (LT (ww).a == 0)
    {
      printf ("vaka\n");
      exit (1);
    }

  // //printpol(o2v(ww));
  printf (" 来ちゃだめなところに来ました\n");

  exit (1);
}


//パターソンアルゴリズムでバイナリGoppa符号を復号する
vec
pattarson (OP w, OP f)
{
  OP g1 = { 0 }
  , ll = {
    0
  };
  int i, j, k, l, c, n, count = 0;
  int flg, o1 = 0;
  OP tt = { 0 }, ff = {
    0
  }, h = {
    0
  };
  EX hh = { 0 };
  vec v;
  oterm rr = { 0 };
  OP r2 = { 0 }, b2 = {
    0
  };
  //unsigned short g[K+1]={2,2,12,1,2,8,4,13,5,10,8,2,15,10,7,3,5};


  tt.t[0].n = 1;
  tt.t[0].a = 1;


  ff = inv (f, w);
//  //printpol (o2v (ff));
  b2 = omod (omul (ff, f), w);
  if (odeg ((b2)) > 0)
    {
      printf ("逆元が計算できません。error\n");
      wait ();
      exit (1);
    }
  printf ("locater==========\n");
  //exit(1);
  r2 = oadd (ff, tt);
  printpol (o2v (r2));
  printf (" h+x==============\n");
  //wait();
  //  exit(1);
  g1 = osqrt (r2, w);
  printpol (o2v (g1));
  printf (" sqrt(h+x)=======\n");
  b2 = omod (omul (g1, g1), w);
  printpol (o2v (b2));
  printf (" g1*g1%w===========\n");
  if (deg (o2v (b2)) == 0)
    exit (1);
  if (LT2 (b2).a != LT2 (r2).a)
    {
      printpol (o2v (w));
      printf (" w============\n");
      printpol (o2v (r2));
      printf (" r2============\n");
      printpol (o2v (g1));
      printf (" g1============\n");
      printf (" g1は平方ではありません。error");
      wait ();
      exit (1);
    }
  printpol (o2v (w));
  printf (" ========goppa\n");
  printpol (o2v (g1));
  printf (" g1!=========\n");
  if (LT (g1).n == 0 && LT (g1).a == 0)
    {
      //printpol (o2v (w));
      printf (" badkey=========\n");
      //printpol (o2v (g1));
      printf (" g1============\n");
      printf ("平方が０になりました。 error\n");
      wait ();
      exit (1);
    }
  //exit(1);
  hh = xgcd (w, g1);
  flg = 0;


  ff = omod (omul (hh.v, g1), w);
  printpol (o2v (ff));
  printf (" beta!=========\n");
  printpol (o2v (w));
  printf (" goppa=========\n");
  printpol (o2v (f));
  printf (" syn=========\n");
  if (odeg ((ff)) != K / 2)
    {
      flg = 1;
      printpol (o2v (ff));
      printf (" locater function failed!! error\n");
      printf ("cannot correct(bad key) error============\n");
      wait ();
      //return -1;
    }


  printpol (o2v (hh.v));
  printf (" alpha!=========\n");
  //exit(1);
  if (odeg ((ff)) == K / 2)
    {
      ll = oadd (omul (ff, ff), omul (tt, omul (hh.v, hh.v)));
    }
  else if (odeg ((ff)) == 1)
    {
      ll = oadd (omul (ff, ff), omul (tt, omul (hh.v, hh.v)));  //ff;
      wait ();
    }
  else
    {
      printf ("locate degree is !=K/2 %d\n", odeg ((ff)));
      exit (1);
    }
  if (odeg ((ll)) == 0)
    {
      printf (" locater degree is 0\n");
      exit (1);
    }
  printf ("あっ、でる・・・！\n");
  count = 0;
  //printpol (o2v (ll));
  //printf (" ll=================\n");
  //printpol (o2v (f));
  //printf (" syn=================\n");
  v = chen (ll);
  /*
     if (count < 2 * T - 1)
     {
     printf (" decode failed error===========\n");
     for (i = 0; i < N; i++)
     {
     if (trace (ll, i) == 0)
     printf ("x=%d\n", i);
     }
     }
   */

  return v;
}


//512bitの秘密鍵を暗号化
void
encrypt (char buf[], unsigned char sk[64])
{
  const uint8_t *hash = { 0 };
  sha3_context c = { 0 };
  int image_size = 512, i;
  FILE *fp;
//  unsigned short dd=0;



  printf ("plain text=");
  for (i = 0; i < 64; i++)
    printf ("%u,", sk[i]);
  printf ("\n");
  //  puts(buf);
  //printf("\n");
  //exit(1);

  //scanf("%s",buf);
  sha3_Init256 (&c);
  sha3_Update (&c, (char *) buf, strlen (buf));
  hash = sha3_Finalize (&c);

  //j=0;

  for (i = 0; i < 64; i++)
    {
      printf ("%d", hash[i]);
      //char s[3];
      //byte_to_hex(hash[i],s);

      sk[i] ^= hash[i];

    }
  printf ("\nencrypt sk=");
  for (i = 0; i < 64; i++)
    printf ("%d,", sk[i]);
  printf ("\n");

  fp = fopen ("enc.sk", "wb");
  fwrite (sy, 2, K, fp);
  fwrite (sk, 1, 64, fp);
  fclose (fp);

}


void
decrypt (OP w)
{
  FILE *fp;
  int i, j;
  unsigned char sk[64] = { 0 }, err[N] = {
    0
  };
  unsigned short buf[K] = { 0 }, tmp[K] = {
    0
  };
  OP f = { 0 };
  vec v = { 0 };
  const uint8_t *hash = { 0 };
  sha3_context c = { 0 };
  int image_size = 512;


  j = 0;
  fp = fopen ("enc.sk", "rb");

  fread (tmp, 2, K, fp);
  fread (sk, 1, 64, fp);
  fclose (fp);

  for (i = 0; i < K; i++)
    buf[i] = tmp[K - i - 1];

  printf ("in decrypt\n");
  f = setpol (buf, K);
  v = pattarson (w, f);
  //v=o2v(h);
  //  exit(1);


  j = 0;
  if (v.x[1] > 0 && v.x[0] == 0)
    {
      err[0] = 1;
      j++;
    }

  printf ("j=%d\n", j);
  printf ("after j\n");
  for (i = j; i < 2 * T; i++)
    {
      if (v.x[i] > 0 && v.x[i] < N)
        {
          err[v.x[i]] = 1;
        }
    }

  char buf0[8192] = { 0 }, buf1[10] = {
    0
  };

  //#pragma omp parallel for
  for (i = 0; i < N; i++)
    {
      snprintf (buf1, 10, "%d", err[i]);
      strcat (buf0, buf1);
    }
  //puts (buf0);
  printf ("vector=%d\n", strlen (buf0));
  //exit(1);
  printf ("cipher sk2=");
  for (i = 0; i < 64; i++)
    printf ("%u,", sk[i]);
  printf ("\n");

  sha3_Init256 (&c);
  sha3_Update (&c, (char *) buf0, strlen (buf0));
  hash = sha3_Finalize (&c);

  j = 0;
  printf ("hash=");
  for (i = 0; i < 64; i++)
    {
      printf ("%d", hash[i]);
      //char s[3];
      //byte_to_hex(hash[i],s);

      sk[i] ^= hash[i];
    }
  printf ("\ndecript sk=");
  for (i = 0; i < 64; i++)
    printf ("%u,", sk[i]);
  printf ("\n");
  //  exit(1);

  return;
}



OP
synd (unsigned short zz[])
{
  unsigned short syn[K] = { 0 }, s = 0;
  int i, j, t1;
  OP f = { 0 };


  printf ("in synd\n");

#pragma omp parallel for        //num_threads(8)
  for (i = 0; i < K; i++)
    {
      syn[i] = 0;
      s = 0;
      //#pragma omp parallel num_threads(8)
      for (j = 0; j < N; j++)
        {
          syn[i] ^= gf[mlt (fg[zz[j]], fg[mat[j][i]])];
        }
      sy[i] = syn[i];           //=s;
      //printf ("syn%d,", syn[i]);
    }
  //printf ("\n");

  for (int j = 0; j < K / 2; j++)
    {
      t1 = syn[j];
      syn[j] = syn[K - j - 1];
      syn[K - j - 1] = t1;
    }


  f = setpol (syn, K);
  //printpol (o2v (f));
  //printf (" syn=============\n");
  //  exit(1);

  return f;
}


//64バイト秘密鍵の暗号化と復号のテスト
void
test (OP w, unsigned short zz[])
{
  int i;
  vec v = { 0 };
  const uint8_t *hash;
  sha3_context c;
  //int image_size=512;
  OP f = { 0 };
  FILE *fp;


  fp = fopen ("aes.key", "rb");
  /*
     static char base64[] = {
     'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H',
     'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P',
     'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X',
     'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f',
     'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n',
     'o', 'p', 'q', 'r', 's', 't', 'u', 'v',
     'w', 'x', 'y', 'z', '0', '1', '2', '3',
     '4', '5', '6', '7', '8', '9', '+', '/',
     };
   */
  char buf[8192] = { 0 }, buf1[10] = {
    0
  };
  unsigned char sk[64] = { 0 };
  // unsigned short s[K]={0};
  //fread(sk,1,32,fp);
  for (i = 0; i < 64; i++)
    sk[i] = i + 1;

  for (i = 0; i < N; i++)
    {
      snprintf (buf1, 10, "%d", zz[i]);
      strcat (buf, buf1);
    }
  //puts (buf);
  printf ("vector=%u\n", strlen (buf));
  //exit(1);


  printf ("sk0=");
  for (i = 0; i < 64; i++)
    printf ("%u,", sk[i]);
  printf ("\n");
  //exit(1);

  f = synd (zz);
  v = o2v (f);
  //printf("v=");
  for (i = 0; i < K; i++)
    {
      sy[i] = v.x[i];
      printf ("%d,", sy[i]);
    }
  printf ("\n");

  encrypt (buf, sk);
  decrypt (w);


  sha3_Init256 (&c);
  sha3_Update (&c, (char *) buf, strlen (buf));
  hash = sha3_Finalize (&c);


}

/*
void trap(OP w,OP f){
  OP hh={0},d,tt,ff,tmp,r2;
  
  hh = gcd (w, f);
      if (odeg ((hh.d)) > 0)
	{
	  printf (" s,wは互いに素じゃありません。\n");
	  wait ();
	  goto label;
	}
      tt.t[0].n = 1;
      tt.t[0].a = 1;
      ff = inv (f, w);
      tmp = omod (omul (ff, f), w);
      if (odeg ((tmp)) > 0)
	{
	  //printpol (o2v (tmp));
	  printf (" inv(h+x)============\n");
	  //printpol (o2v (w));
	  printf (" w============\n");
	  printf ("この多項式では逆元計算ができません。");
	  printf ("count=%d\n", k);
	  wait ();
	  //  goto label;
	}
      r2 = oadd (ff, tt);
      //printpol (o2v (r2));
      printf (" h+x==============\n");
      //  exit(1);
      g1 = osqrt (r2, w);
      //printpol (o2v (g1));
      printf (" g1!=========\n");
      r1 = omod (omul (g1, g1), w);
      //printpol (o2v (r1));
      printf (" g1^2 mod w===========\n");
      //printpol (o2v (r2));
      printf (" r2===========same?\n");
      //scanf("%d",&n);
      if (odeg ((r1)) != odeg ((r2)) && odeg ((g1)) > 0)
	{
	  //printpol (o2v (w));
	  printf (" w===========\n");
	  //printpol (o2v (r2));
	  printf (" r2===========\n");
	  printf ("平方根の計算に失敗しました。\n");
	  printf ("count=%d\n", k);
	  wait ();
	  goto label;
	  //exit(1);
	}
      if (odeg ((g1)) == 0)
	{
	  //printpol (o2v (g1));
	  printf (" sqrt(h+x)==============\n");
	  //printpol (o2v (w));
	  printf (" badkey=========\n");
	  printf ("平方根が０になりました。\n");
	  printf ("count=%d\n", k);
	  wait ();
	  //exit(1);
	  goto label;
	}
      hh = xgcd (w, g1);
      ff = omod (omul (hh.v, g1), w);
      //printpol (o2v (ff));
      printf (" beta!=========\n");
      if (odeg ((ff)) != K / 2)
	{
	  printf ("deg(l)!=K/2=%d %d %d\n", odeg ((ff)), K, k);
	  //exit(1);
	  goto label;
	}
}
*/


void
readkey ()
{

  /*
     //鍵をファイルに書き込むためにはkey2を有効にしてください。
     key2 (g);

     fp=fopen("sk.key","rb");
     fread(g,2,K+1,fp);
     fclose(fp);

     //固定した鍵を使いたい場合はファイルから読み込むようにしてください。  

     fq = fopen ("H.key", "rb");
     fread (dd, 2, K * N, fq);
     //#pragma omp parallel for
     for (i = 0; i < N; i++)
     {
     for (j = 0; j < K; j++)
     mat[i][j] = dd[K * i + j];
     }
     fclose (fq);
   */

}


//言わずもがな
int
main (void)
{
  int i, j, k, l;
  int count = 0;
  FILE *fp, *fq;
  unsigned short z1[N] = { 0 };
  int flg, o1 = 0;
  OP f = { 0 }, r = {
    0
  }, w = {
    0
  }, ff = {
    0
  }, tt = {
    0
  };
  EX hh = { 0 };
  vec v;
  //  unsigned short dd[K*N] = {0},gg[K+1]={0};
  time_t t;
  OP r1 = { 0 }, r2 = {
    0
  };
  OP g1 = { 0 }, tmp = {
    0
  };
  int fail = 0;


//公開鍵matはグローバル変数でスタック領域に取る。
//ヒープ領域は使うときはここを有効にしてください。
/*
  mat = malloc (N * sizeof (unsigned short *));
  base=malloc(sizeof(unsigned short)*K*N);
  #pragma omp parallel for
  for(i=0;i<N;i++){
  //連続したメモリ空間にデータを配置
  mat[i]=base+i*K;
  memset(mat[i],0,K);
  }
*/

#ifdef SRAND
  srand (SRAND);
#else
  const unsigned int seed = clock () + time (&t);
  printf ("srand(%u)\n", seed);
  srand (seed);
#endif

label:

  //makeS();
  //exit(1);
  do
    {
      fail = 0;
      k = 0;
      flg = 0;
      memset (g, 0, sizeof (g));
      memset (ta, 0, sizeof (ta));
      ginit ();

      for (i = 0; i < K + 1; i++)
        {
          if (g[K - 1] > 0)
            flg = 1;
          if (i % 2 == 1 && g[i] > 0 && i < K)
            k++;
        }
      if ((k > 0 && flg == 0) || (k > 1 && flg == 1))
        {
          w = setpol (g, K + 1);
          oprintpol (w);
        }

      //w = setpol (g, K + 1);
      //oprintpol (w);

      //多項式の値が0でないことを確認
      for (i = 0; i < D; i++)
        {
          ta[i] = trace (w, i);
          if (ta[i] == 0)
            {
              printf ("trace 0 @ %d\n", i);
              fail = 1;
              break;
            }
        }

    }
  while (fail);

#pragma omp parallel for
  for (i = 0; i < N; i++)
    tr[i] = oinv (ta[i]);


  memset (mat, 0, sizeof (mat));

  printf ("\nすげ、オレもうイキそ・・・\n");
  //keygen(g);



  //パリティチェックを生成する。
  //パリティチェックに0の列があったら、なくなるまでやり直す。
  do
    {
      i = deta (g);
    }
  while (i < 0);



lab:

  matmul ();
  matinv ();

//decode開始
  k = 0;
  while (1)
    {

      o1 = 0;

      count = 0;

      memset (zz, 0, sizeof (zz));


      j = 0;
      while (j < T)
        {
          l = xor128 () % D;
          //printf("l=%d\n",l);
          if (0 == zz[l] && l > 0)
            {
              zz[l] = l;
              j++;
            }
        }

      for (i = 0; i < D; i++)
        {
          if (zz[i] > 0)
            printf ("l=%d %d\n", i, zz[i]);
        }
      //exit(1);

      f = synd (zz);

      count = 0;
      /*
         count = 0;
         for (i = 0; i < N; i++)
         {
         if (zz[i] > 0)
         count++;
         }
         printf ("%d\n", count);
       */

      r = decode (w, f);

      for (i = 0; i < T; i++)
        {
          if (r.t[i].a > 0 && count > 0)        // == r.t[i].n)
            {
              printf ("e=%d %d %s\n", r.t[i].a, r.t[i].n, "お");
              count++;
            }
          if (count == 0 && r.t[i].a > 0)
            {
              printf ("\ne=%d %d %s\n", r.t[i].a, r.t[i].n, "う");
              count++;
            }

        }
      if (count != T)
        {
          printf ("error pattarn too few %d\n", count);
          exit (1);
        }

      o1 = 0;
      //  #pragma omp parallel for
      for (i = 0; i < D; i++)
        {
          if (zz[i] > 0)
            o1++;
        }
      printf ("err=%dっ！！\n", count);
      //exit(1);
      //goto label;



      //printf("パターソンアルゴリズムを実行します。何か数字を入れてください。\n");
      //wait();

      memset (z1, 0, sizeof (z1));


      j = 0;
      while (j < T * 2)
        {
          l = xor128 () % D;
          //printf ("l=%d\n", l);
          if (0 == z1[l])
            {
              z1[l] = 1;
              j++;
            }
        }

      //encryotion
      //test (w, z1);

      for (i = 0; i < D; i++)
        printf ("%d= %d,\n", i, z1[i]);
      printf ("\n");

      f = synd (z1);


      //バグトラップのためのコード（省略）
      //trap(w,f);
      //バグトラップ（ここまで）

      count = 0;
      //復号化の本体
      v = pattarson (w, f);
      //エラー表示
      for (i = 0; i < T * 2; i++)
        {
          if (i == 0)
            printf ("error position=%d %d う\n", i, v.x[i]);
          if (i > 0 && v.x[i] > 0)
            printf ("error position=%d %d お\n", i, v.x[i]);
          if (i > 0 && v.x[i] == 0)
            {
              printf ("baka %d %d\n", i, v.x[i]);
              exit (1);
            }

          count++;
        }

      printf ("err=%dっ!! \n", count);
      if (count != T * 2)
        printf ("error is too few\n");

      //exit(1);
      //goto lab;
      //wait();

      break;
    }

  return 0;
}
