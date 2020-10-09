///////////////////////////////////////////////////////////////////////////
//
// (over finite field) Gauss-Jordan法による逆行列
//

#define MAXN 4
#define N 16
#define M 32
#define F 32 //2040

//elements of GF16
unsigned char gf[N]={0,1,2,4,8,9,11,15,7,14,5,10,13,3,6,12};
//index of GF16
unsigned char fg[N]={0,1,2,13,3,10,14,8,4,5,11,6,15,12,9,7};
unsigned char c[M][M]={0};

unsigned char a[M][M]; //={{1,2,0,1},{1,1,2,0},{2,0,1,1},{1,2,1,1}}; //入力用の配列
unsigned char cc[M][M]={0};

int mlt(int x, int y){

    if(x==0||y==0)
        return 0;

  return ((x+y-2)%(N-1))+1;
}


int mltn(int n,int x){
  int i,j;

  if(n==0)
    return 1;
  i=x;
    for(j=0;j<n-1;j++)
      i=mlt(i,x);

  return i;
}

int Inv(unsigned char b){
  int i;

      for(i=0;i<N;i++){
        if(gf[mlt(i,b)]==1)
          return i;
      }

}


int det(){
  //double a[M][M]={{2,-2,4,2},{2,-1,6,3},{3,-2,12,12},{-1,3,-4,4}};
double det=1.0,buf;

int n=M;  //配列の次数
int i,j,k;

// lab:
/*
 for(i=0;i<M;i++){
   for(j=0;j<M;j++)
     a[i][j]=rand()%16;
 }
*/
//三角行列を作成
for(i=0;i<n;i++){
 for(j=0;j<n;j++){
  if(i<j){
    //    if(a[i][i]==0)
    //   goto lab;
   buf=a[j][i]/a[i][i];
   for(k=0;k<n;k++){
   a[j][k]-=a[i][k]*buf;
   }
  }
 }
}
 
//対角部分の積
for(i=0;i<n;i++){
 det*=a[i][i];
}
 
printf("%f\n",det); // -> 120.000000
 if(det==0.0){
   printf("baka\n");
   exit(1);
 }
 if(det==1.0)
 return 0;
 if(det!=1.0)
   return -1;
}


void g2(){
  int i,j,k;
unsigned char bb[F][F]={0};

#pragma omp parallel for    
  for(i=0;i<F;i++){
    a[i][i]=rand()%16;
    bb[i][i]=1;
  }
  for(i=0;i<F;i++){
    //#pragma omp parallel for  
    for(j=i+1;j<F;j++){
      a[i][j]=rand()%16;
      // printf("%d,",a[i][j]);
    }
    //printf("\n");
  }
  for(i=0;i<F;i++){
    for(j=0;j<F;j++)
      printf("%d,",a[i][j]);
    printf("\n");
  }

  for(i=0;i<F;i++){
    //#pragma omp parallel for  
    for(j=i+1;j<F;j++){
      bb[j][i]=rand()%16;      
    }

  }
  for(i=0;i<F;i++){
    for(j=0;j<F;j++)
      printf("%d,",bb[i][j]);
    printf("\n");
  }
  //exit(1);
  for(i=0;i<F;i++){
    for(j=0;j<F;j++){
      //#pragma omp parallel for  
      for(k=0;k<F;k++){
	cc[i][j]^=bb[i][k]&a[k][j];
      }
    }
  }
  for(i=0;i<F;i++){
    for(j=0;j<F;j++)
      printf("%d,",cc[i][j]);
    printf("\n");
  }
  //exit(1); 
}


int matinv(){
  //unsigned char a[M][M]; //={{1,2,0,1},{1,1,2,0},{2,0,1,1},{1,2,1,1}}; //入力用の配列
unsigned char inv_a[M][M]; //ここに逆行列が入る
unsigned char buf; //一時的なデータを蓄える
 unsigned char b[M][M]={0},dd[M][M]={0};
 int i,j,k,count; //カウンタ
 int n=M;

 lab:
 
 for(i=0;i<n;i++){
   for(j=0;j<n;j++){
     a[i][j]=rand()%16;
     printf("%d,",a[i][j]);
   }
    printf("\n");
 }
 

 // printf("\n");
 for(i=0;i<n;i++){
   for(j=0;j<n;j++)
     c[i][j]=a[i][j];
 }
//単位行列を作る
for(i=0;i<n;i++){
 for(j=0;j<n;j++){
 inv_a[i][j]=(i==j)?1:0;
 }
}
//掃き出し法
for(i=0;i<n;i++){
  buf=gf[Inv(fg[a[i][i]])];
 for(j=0;j<n;j++){
   a[i][j]=gf[mlt(fg[buf],fg[a[i][j]])];
   inv_a[i][j]=gf[mlt(fg[buf],fg[inv_a[i][j]])];
}
for(j=0;j<n;j++){
 if(i!=j){
  buf=a[j][i];
  for(k=0;k<n;k++){
    a[j][k]^=gf[mlt(fg[a[i][k]],fg[buf])];
    inv_a[j][k]^=gf[mlt(fg[inv_a[i][k]],fg[buf])];
  }
 }
}
}
//逆行列を出力
for(i=0;i<n;i++){
  count=0;
 for(j=0;j<n;j++){
   if(inv_a[i][j]==0)
     count++;
   if(count==n)
     goto lab;
  printf(" %d",inv_a[i][j]);
 }
 printf("\n");
}
//exit(1); 
//検算
 for(i=0;i<n;i++){
   for(j=0;j<n;j++){
     for(k=0;k<n;k++)
       b[i][j]^=gf[mlt(fg[c[k][j]],fg[inv_a[i][k]])];

     printf("%d,",b[i][j]);
     // if(j==i && b[i][j]!=1 && j!=i && b[i][j]>0)
     //goto lab;
       
   }
   printf("\n");
 }
 //exit(1);
 
 return 0;
}


int main(){
  
    int i,j;
    double b[4],k=0;

    srand(clock());

    //g2();
 lab:
    matinv();
    //if(det()!=1.0)
    //goto lab;
    
    return 0;
}
