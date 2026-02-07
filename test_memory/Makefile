CC = gcc
CFLAGS = -O3 -Wall -Wextra -std=c99 -g

# Paths
COMMON = ../common
MLDSA65 = ../crypto_sign/ml-dsa-65/clean

# Source files
MLDSA65_SRCS = $(MLDSA65)/sign.c $(MLDSA65)/packing.c $(MLDSA65)/polyvec.c \
               $(MLDSA65)/poly.c $(MLDSA65)/ntt.c $(MLDSA65)/reduce.c \
               $(MLDSA65)/rounding.c $(MLDSA65)/symmetric-shake.c

COMMON_SRCS = $(COMMON)/fips202.c $(COMMON)/randombytes.c

TEST_SRC = test_mldsa65_memory.c
STACK_SRC = test_mldsa65_stack_analysis.c

all: test_mldsa65_memory test_mldsa65_stack_analysis

test_mldsa65_memory: $(TEST_SRC) $(MLDSA65_SRCS) $(COMMON_SRCS)
	$(CC) $(CFLAGS) -I$(COMMON) -I$(MLDSA65) -o $@ $^

test_mldsa65_stack_analysis: $(STACK_SRC) $(MLDSA65_SRCS) $(COMMON_SRCS)
	$(CC) $(CFLAGS) -I$(COMMON) -I$(MLDSA65) -o $@ $^

run: test_mldsa65_memory
	./test_mldsa65_memory

# Run with valgrind massif for heap profiling
massif: test_mldsa65_memory
	valgrind --tool=massif --stacks=yes --massif-out-file=massif.out ./test_mldsa65_memory
	ms_print massif.out

# Run with valgrind memcheck
memcheck: test_mldsa65_memory
	valgrind --tool=memcheck --leak-check=full ./test_mldsa65_memory

clean:
	rm -f test_mldsa65_memory massif.out massif.out.*

.PHONY: all run massif memcheck clean
