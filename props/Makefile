CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: props

props: src/props.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: props
	./props out

clean:
	rm -f props
	rm -rf out
