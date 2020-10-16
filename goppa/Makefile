files=op.c main.c
dev:
	clang-10 -w -g -rdynamic -pg -O3 -mtune=znver2 -march=znver2 -ffast-math -funroll-loops  -fopenmp main.c

clang:
	clang-10 -w -g -pg -O3 -mtune=znver2 -march=znver2 -ffast-math -funroll-loops  -fopenmp main.c

gcc:
	gcc -Wall -g -pg -O3 -mtune=znver2 -march=znver2 -ffast-math -funroll-loops  -fopenmp main.c

stable:
	gcc -Wall -O3 -g -pg -mtune=native -march=native -ffast-math -funroll-loops -fopenmp oplib.c

format:
	VERSION_CONTROL=none indent -nut $(files)

clean:
	rm -f a.out
