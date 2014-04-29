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

NAME = $(dir)/matcreator

all: update base

base: matcreator.c 
	$(CC) -o $(NAME) matcreator.c

strict: matcreator.c 
	$(CC) $(STRICT) -o $(NAME) matcreator.c

fast_s:  matcreator.c
	$(CC) $(FAST) -o $(NAME) matcreator.c

fast: update fast_s

clean:
	rm $(dir)/*.out

update: 
	svn update
