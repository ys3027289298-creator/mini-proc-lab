CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: lod

lod: src/lod.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: lod
	./lod --selftest

run: lod
	mkdir -p ../artifacts/lod
	./lod ../artifacts/lod | tee ../artifacts/lod/metrics.txt
	python3 tools/assemble.py ../artifacts/lod

clean:
	rm -f lod
