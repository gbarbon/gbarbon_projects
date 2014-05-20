# CC says that the variable CC will be the compiler to use.
CC=mpicc

# OPTIMIZED will be the optimization options passed to the compiler.
# -Ofast := all optimizations
# -mark=native := activate machine instructions sets (only for Lynux systems)
# -mtune=core-avx-i := for Intel_Core OSX architectures
# -msse4.2 := for non-Intel_Core OSX architectures

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

NAME = $(dir)/hpc_linear_mm
CFILENAME = linear_mm.c

all: update base

base: 
	$(CC) -o $(NAME) $(CFILENAME) -lm

strict: 
	$(CC) $(STRICT) -o $(NAME) $(CFILENAME) -lm

opti_s:  
	$(CC) $(OPTIMIZED) -o $(NAME)-op $(CFILENAME) -lm -DOPTI='"-op"'

opti: update opti_s

clean:
	rm $(dir)/*.out

update: 
	svn update
