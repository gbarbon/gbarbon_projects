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

NAME = $(dir)/hpc_mpi_mm

all: update base

base: mm.c 
	$(CC) -o $(NAME) mm.c

strict: mm.c 
	$(CC) $(STRICT) -o $(NAME) mm.c

fast:  mm.c
	$(CC) $(FAST) -o $(NAME) mm.c

clean:
	rm $(dir)/*.out
	
update: 
	svn update
