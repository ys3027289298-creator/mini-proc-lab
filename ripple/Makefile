CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra
SHELL := /bin/bash

all: ripple

ripple: src/ripple.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: ripple
	./ripple

matrix: ripple
	mkdir -p out
	set -o pipefail; \
	python3 -I tools/run_matrix.py | tee out/fixed_results.txt; \
	python3 -I tools/run_matrix.py --coupled | tee out/coupled_results.txt || true

montage: matrix
	python3 -I tools/montage.py fixed
	python3 -I tools/montage.py coupled

clean:
	rm -rf ripple out
