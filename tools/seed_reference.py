#!/usr/bin/env python3
"""Independent reference implementation of the mini-proc-lab seed derivation.

Keep this file in sync with common/seed.c and docs/seed-derivation.md.
Run:  python3 -I tools/seed_reference.py [master ...]
"""
import sys

MASK = (1 << 64) - 1
GOLDEN = 0x9E3779B97F4A7C15
FNV_OFFSET = 0xCBF29CE484222325
FNV_PRIME = 0x100000001B3
PREFIX = b"mini-proc-lab/v1/"
LAYERS = ("noise", "scatter", "erosion", "biome", "hydro")


def mix64(z):
    z = (z + GOLDEN) & MASK
    z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & MASK
    z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & MASK
    return (z ^ (z >> 31)) & MASK


def fnv1a64(data):
    h = FNV_OFFSET
    for b in data:
        h ^= b
        h = (h * FNV_PRIME) & MASK
    return h


def subseed(master, layer):
    return mix64((master & MASK) ^ fnv1a64(PREFIX + layer.encode("ascii")))


def derive_all(master):
    return {name: subseed(master, name) for name in LAYERS}


def main(argv):
    masters = [int(a, 0) for a in argv[1:]] or [0, 1, 123456789, 0xDEADBEEFCAFEF00D]
    for master in masters:
        print(f"master {master} 0x{master & MASK:016x}")
        for name, value in derive_all(master).items():
            print(f"  {name:<8} {value:20d}  0x{value:016x}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
