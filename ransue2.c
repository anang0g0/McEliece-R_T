
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#define N 16
#define NN 2


unsigned char x0[N];//={1,2,3,4,5,6,7,0};
unsigned char x1[N];//={2,3,1,6,5,7,0,4};
unsigned char x2[N]={0};

 void rp(unsigned char* a) {
        int i, j, x;
        time_t t;

	srand(clock()+time(&t));


        for (i = 0; i < N; i++) {
            a[i] = i;
        }
        for (i = 0; i < N - 2; i++) {
            // rand from i+1 to N-1
            j = (rand() % (N - 1 - i)) + i + 1;

            // swap a[i] and a[j]
            x = a[j];
            a[j] = a[i];
            a[i] = x;
        }
        if (a[N - 1] == N - 1) {
            a[N - 1] = a[N - 2];
            a[N - 2] = N - 1;
        }


    }




unsigned int xor(void){
  static unsigned int y = 2463534242;
  y = y ^ (y << 13); y = y ^ (y >> 17);
  return y = y ^ (y << 15);
}

unsigned long long int xor64(void) {
  static unsigned long long int x = 88172645463325252ULL;
  x = x ^ (x << 13); x = x ^ (x >> 7);
  return x = x ^ (x << 17);
}

void data(){
  unsigned long long int i,j=0,k=0;
  unsigned char a[N]={1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  
  unsigned short aa=0;
  
  unsigned char z[N];
  unsigned char w[N];
  FILE *fp;
  int cnt=0;
  unsigned char ww[N]={0};
  unsigned int flg[N]={0},count=0;

  for(i=0;i<N;i++)
    x2[x0[i]]=i;
  fp=fopen("rand.dat","wb");

  //  for(i=0;i<N;i++)
  //a[i]=xor();
  j=0;
  i=0;
  while(j<10000000){
    count=0;
    // printf("=======\n");
    for(i=0;i<N;i++){
      /*
      if(i<x1[i]){
	
	//    printf("%u,%u\n",i,x1[i]);
      }
      */
      // if(i>x1[i]){
      a[i]^=a[x1[i]];
	//printf("*%u,%u\n",i,x1[i]);
	//count++;
      }
      
    //    }
    //    printf("%u\n",count);
    //for(i=0;i<N;i++)
    //a[i]^=w[i];
    aa=0;
    
    for(cnt=0;cnt<NN;cnt++){
      
      for(i=cnt*8;i<8*cnt+8;i++){
	aa=aa<<1;
	aa^=a[i];
      }
    }

    
    fwrite(&aa,NN,1,fp);    
    //printf("%llu,\n",aa);
        
    for(i=0;i<N;i++)
      w[i]=x0[x1[x2[i]]];
    
    for(i=0;i<N;i++)
      x1[i]=w[i];
    
    
    j++;
  }

  fclose(fp);

}



int main(){
  int i;
    rp(x0);
    rp(x1);

    for(i=0;i<N;i++)
      x2[x0[i]]=i;

    for(i=0;i<N;i++)
      printf("%d,",x0[i]);
    printf("\n");
    for(i=0;i<N;i++)
      printf("%d,",x2[i]);
    printf("\n");
    
    //  data();
  
  return 0;
}
