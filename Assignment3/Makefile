EXECUTABLES = jacobi-omp gs-omp
COMPILER = gcc
FLAGS = -fopenmp

all: $(EXECUTABLES)

jacobi-omp: jacobi-omp.c
	$(COMPILER) $(FLAGS) jacobi-omp.c -o jacobi-omp

gs-omp: gs-omp.c
	$(COMPILER) $(FLAGS) gs-omp.c -o gs-omp

clean:
	rm -rf $(EXECUTABLES)
