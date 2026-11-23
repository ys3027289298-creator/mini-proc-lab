CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: scatter bluenoise

scatter: src/scatter.c
	$(CC) $(CFLAGS) $< -o $@ -lm

bluenoise: src/bluenoise.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: scatter bluenoise
	./scatter
	./bluenoise selftest

clean:
	rm -f scatter bluenoise
