PY ?= python3
OUT ?= ../artifacts/bandpass

all: run

run:
	$(PY) src/bandpass.py $(OUT)

# 验收：流水线在任意一档切频超 20% 红线时以非零状态退出
test: run

clean:
	rm -rf $(OUT)

.PHONY: all run test clean
