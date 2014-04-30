# CC says that the variable CC will be the compiler to use.
CC=mpicc

# CFLAGS will be the options I'll pass to the compiler.
# -Ofast := all optimizations
# -mark=native := activate machine instructions sets
FAST=-Ofast -march=native

# -Wall := all warnings
# -Werror := error on warning
STRICT=-ansi -pedantic -Wall -Werror

dir = bin

NAME = $(dir)/hpc_mpi_cannon

all: update base

base: cannon.c 
	$(CC) -o $(NAME) cannon.c -lm

strict: cannon.c 
	$(CC) $(STRICT) -o $(NAME) cannon.c -lm

fast_s:  cannon.c
	$(CC) $(FAST) -o $(NAME) cannon.c -lm

fast: update fast_s

clean:
	rm $(dir)/*.out

update: 
	svn update
