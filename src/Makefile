
all:
	clang -Wall -g -pg -O6 -mtune=native -march=native -ffast-math -funroll-loops  -fopenmp op.c

amd:
	clang -w -g -pg -O3 -DNDEBUG -mtune=znver1 -march=znver1 -ffast-math -funroll-loops -fopenmp op.c

test:
	clang -Wall -g -pg -O6 -mtune=znver1 -march=znver1  -ffast-math -funroll-loops -fopenmp op.c

gcc:
	gcc -Wall -g -pg  -O2 -DSRAND=1602886277 -mtune=native -march=native -ffast-math -funroll-loops  -fopenmp op.c

oo:
	clang-10 -Wall -g -pg -O3 -DSRAND=1602886277 -mtune=znver2 -march=znver2 -ffast-math -funroll-loops  -fopenmp o.c


stable:
	gcc -Wall -O3 -g -pg -mtune=native -march=native -ffast-math -funroll-loops -fopenmp oplib.c

clean:
	rm -f a.out
