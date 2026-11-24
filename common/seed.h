#ifndef SP_SEED_H
#define SP_SEED_H

#include <stddef.h>
#include <stdint.h>

/*
 * One master seed, five layer sub-seeds.
 *
 * The derivation is a public contract: see docs/seed-derivation.md and the
 * independent reference in tools/seed_reference.py.  Changing any constant,
 * the prefix or the layer name strings invalidates every recorded world.
 */

#define SP_LAYER_COUNT 5
#define SP_SEED_PREFIX "mini-proc-lab/v1/"

typedef enum {
    SP_LAYER_NOISE = 0,
    SP_LAYER_SCATTER = 1,
    SP_LAYER_EROSION = 2,
    SP_LAYER_BIOME = 3,
    SP_LAYER_HYDRO = 4
} SpLayer;

extern const char *const sp_layer_names[SP_LAYER_COUNT];

/* SplitMix64 finalizer, including the 0x9e3779b97f4a7c15 increment. */
uint64_t sp_mix64(uint64_t z);

/* FNV-1a, 64 bit, over raw bytes. */
uint64_t sp_fnv1a64(const void *data, size_t len);
uint64_t sp_fnv1a64_str(const char *s);

/* sp_mix64(master XOR FNV1a64(SP_SEED_PREFIX || layer_name)) */
uint64_t sp_subseed(uint64_t master, const char *layer_name);

/* Fills out[SP_LAYER_*] in enum order. */
void sp_derive_all(uint64_t master, uint64_t out[SP_LAYER_COUNT]);

#endif
