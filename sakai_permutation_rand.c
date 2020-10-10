#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#define N 32
#define BYTE
//BIT_VERSION

void random_permutation(unsigned char* a){
	int i,j,x;
	for(i = 0; i < N; i++){
		a[i] = i;
	}
	for(i = 0; i < N - 2; i++){
		// rand from i+1 to N-1
		j = (rand() % (N-1-i)) + i + 1;

		// swap a[i] and a[j]
		x = a[j];
		a[j] = a[i];
		a[i] = x;
	}
	if(a[N-1] == N-1){
		a[N-1] = a[N-2];
		a[N-2] = N - 1;
	}
}

unsigned int toInt(unsigned char * a){
	unsigned int i = 0, s = 1, ret = 0;
	for(i=0; i < N; i++){
		ret += s * a[i];
		s *= 2;
		ret << 1;
	}
	return ret;
}

int p_rand(){
#ifdef BIT_VERSION
	static unsigned char a[N] = {
		1,1,1,1,1,1,1,1,
		1,1,1,1,1,1,1,1,
		0,0,0,0,0,0,0,0,
		0,0,0,0,0,0,0,0,
	};
#else
	static unsigned char a[N] = {
		1,2,3,5,7,11,13,17,
		19,23,29,31,33,37,39,41,
		143,151,157,159,161,167,171,175,
		212,85,1,0,0,0,0,0,
	};
#endif

	static unsigned char x[N];
	static unsigned char y[N];
	static unsigned char inv_x[N];
	static unsigned char first = 1;

	static unsigned char tmp[N];
	int i;

	if(first){
		first = 0;

		random_permutation(x);
		random_permutation(y);

		for(i=0;i<N;i++){
			inv_x[x[i]]=i; 
		}		
	}
	for(i=0;i<N;i++)
	  printf("%d,",x[i]);
	printf("\n");
	for(i=0;i<N;i++)
	  printf("%d,",inv_x[i]);
	printf("\n");
	exit(1);

	// a ^= a * y
	memcpy(tmp, a, sizeof(tmp)); 	// tmp = a;
	for(i=0;i<N;i++){ 
		a[i] ^= tmp[y[i]];
	}

	// y = x * y * ~x
	for(i=0;i<N;i++) {
		tmp[i]=x[y[inv_x[i]]];
	}
	memcpy(y, tmp, sizeof(tmp));

	return toInt(a);
}

int main(){ 
	FILE *fp = fopen("out.txt","wb");  
	int i,j=0; 

	//srand((unsigned int)time(NULL)); 

	while(j<32){
		unsigned int a = p_rand();

		if(fp){
			fprintf(fp, "%u\n", a);
		}else{
			printf("%u\n", a);
		}

		j++; 
	} 
	return 0;
} 
