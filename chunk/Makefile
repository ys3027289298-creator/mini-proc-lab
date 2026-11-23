CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra
CPPFLAGS += -I../noise/src

all: chunk

chunk: src/chunk.c ../noise/src/noise.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -o $@ -lm

test: chunk
	./chunk

clean:
	rm -f chunk
