CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: vibfield

vibfield: src/vibfield.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: vibfield
	./vibfield

runs: vibfield
	./vibfield runs FILTERED 10 > /tmp/run_filtered.csv
	./vibfield runs RAW 10 > /tmp/run_raw.csv
	python3 -I tools/repeat_report.py /tmp/run_filtered.csv /tmp/run_raw.csv \
		../artifacts/oob_repeat_report.md

clean:
	rm -f vibfield

.PHONY: test runs clean
