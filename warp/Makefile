CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: warp

warp: src/warp.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: warp
	mkdir -p out
	./warp 0.0 out/base.bmp
	./warp 0.3333 out/a1.bmp
	./warp 0.6667 out/a2.bmp
	./warp 1.0 out/a3.bmp

clean:
	rm -f warp
	rm -rf out
