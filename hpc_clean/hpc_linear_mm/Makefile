# CC says that the variable CC will be the compiler to use.
CC=gcc

# CFLAGS will be the options I'll pass to the compiler.
# -Ofast := all optimizations
# -mark=native := activate machine instructions sets
FAST=-Ofast -march=native

# -Wall := all warnings
# -Werror := error on warning
STRICT=-ansi -pedantic -Wall -Werror

dir = bin

NAME = $(dir)/hpc_linear_mm

all: update base

base: linear_mm.c 
	$(CC) -o $(NAME) linear_mm.c -lm

strict: linear_mm.c 
	$(CC) $(STRICT) -o $(NAME) linear_mm.c -lm

fast_s:  linear_mm.c
	$(CC) $(FAST) -o $(NAME) linear_mm.c -lm

fast: update fast_s

clean:
	rm $(dir)/*.out

update: 
	svn update
