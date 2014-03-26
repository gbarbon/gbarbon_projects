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

NAME = $(dir)/hpc_mpi_farm

all: update base

base: farm.c 
	$(CC) -o $(NAME) farm.c

strict: farm.c 
	$(CC) $(STRICT) -o $(NAME) farm.c

fast_s:  farm.c
	$(CC) $(FAST) -o $(NAME) farm.c

fast: update fast_s

clean:
	rm $(dir)/*.out

update: 
	svn update
