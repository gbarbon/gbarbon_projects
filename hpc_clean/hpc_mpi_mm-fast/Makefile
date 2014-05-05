# CC says that the variable CC will be the compiler to use.
CC=mpicc

# CFLAGS will be the options I'll pass to the compiler.
# -Ofast := all optimizations
# -mark=native := activate machine instructions sets

# uname: -s: System name, -p: processor type, -m: machine type
UNAME_S := $(shell uname -s)
UNAME_P := $(shell uname -p)
UNAME_M := $(shell uname -m)

OPTIMIZED=-Ofast
ifeq ($(UNAME_S),Linux)
	OPTIMIZED += -march=native
endif
ifeq ($(UNAME_S),Darwin)
	ifeq ($(UNAME_M),x86_64)
		OPTIMIZED += -mtune=core-avx-i
	else
		OPTIMIZED += -msse4.2
	endif
endif

# -Wall := all warnings
# -Werror := error on warning
STRICT=-ansi -pedantic -Wall -Werror

dir = bin

NAME = $(dir)/hpc_mpi_mm-fast

all: update base

base: mmfast.c 
	$(CC) -o $(NAME) mmfast.c -lm -DOPTI='""'

strict: mmfast.c 
	$(CC) $(STRICT) -o $(NAME) mmfast.c -lm -DOPTI='""'

opti_s:  mmfast.c
	$(CC) $(OPTIMIZED) -o $(NAME)-op mmfast.c -lm -DOPTI='"-op"'

opti: update opti_s

clean:
	rm $(dir)/*.out

update: 
	svn update
