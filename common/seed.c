#include "seed.h"

#include <string.h>

const char *const sp_layer_names[SP_LAYER_COUNT] = {
    "noise", "scatter", "erosion", "biome", "hydro"
};

uint64_t sp_mix64(uint64_t z) {
    z += UINT64_C(0x9E3779B97F4A7C15);
    z = (z ^ (z >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)) * UINT64_C(0x94D049BB133111EB);
    return z ^ (z >> 31);
}

static uint64_t fnv_update(uint64_t h, const void *data, size_t len) {
    const unsigned char *p = (const unsigned char *)data;
    for (size_t i = 0; i < len; i++) {
        h ^= (uint64_t)p[i];
        h *= UINT64_C(0x100000001B3);
    }
    return h;
}

uint64_t sp_fnv1a64(const void *data, size_t len) {
    return fnv_update(UINT64_C(0xCBF29CE484222325), data, len);
}

uint64_t sp_fnv1a64_str(const char *s) {
    return sp_fnv1a64(s, strlen(s));
}

uint64_t sp_subseed(uint64_t master, const char *layer_name) {
    uint64_t h = UINT64_C(0xCBF29CE484222325);
    h = fnv_update(h, SP_SEED_PREFIX, sizeof(SP_SEED_PREFIX) - 1u);
    h = fnv_update(h, layer_name, strlen(layer_name));
    return sp_mix64(master ^ h);
}

void sp_derive_all(uint64_t master, uint64_t out[SP_LAYER_COUNT]) {
    for (int i = 0; i < SP_LAYER_COUNT; i++) {
        out[i] = sp_subseed(master, sp_layer_names[i]);
    }
}
