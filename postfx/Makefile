CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra
ART ?= ../artifacts/postfx
PYTHON ?= python3

all: postfx

postfx: src/postfx.c
	$(CC) $(CFLAGS) $< -o $@ -lm

$(ART):
	mkdir -p $(ART)

test: postfx | $(ART)
	./postfx $(ART)/postfx.bin | tee $(ART)/report.txt
	$(PYTHON) -I tools/render.py $(ART)/postfx.bin $(ART)

clean:
	rm -f postfx
