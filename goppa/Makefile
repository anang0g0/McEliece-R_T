files=op.c main.c
debug:
	clang-10 -w -g -rdynamic -pg -O3 -DSRAND=1602886277 -mtune=znver2 -march=znver2 -ffast-math -funroll-loops -fopenmp main.c

release:
	clang-10 -w -g -pg -O3 -DNDEBUG -mtune=znver2 -march=znver2 -ffast-math -funroll-loops -fopenmp main.c

gcc:
	gcc -Wall -g -pg -O3 -mtune=znver2 -march=znver2 -ffast-math -funroll-loops  -fopenmp main.c

stable:
	gcc -Wall -O3 -g -pg -mtune=native -march=native -ffast-math -funroll-loops -fopenmp oplib.c

format:
	VERSION_CONTROL=none indent -nut $(files)

clean:
	rm -f a.out
