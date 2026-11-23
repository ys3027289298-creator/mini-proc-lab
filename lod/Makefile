CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

OUT ?= ../artifacts/lod

all: lod

lod: src/lod.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: lod
	mkdir -p $(OUT)
	./lod $(OUT)

clean:
	rm -f lod
