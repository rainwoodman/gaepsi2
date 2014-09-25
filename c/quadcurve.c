#include <stdint.h>
#include <stdio.h>

#include "gaepsi.h"

#define BITS 3 /* < 32, > 1 */

uint32_t translate_table[UINT16_MAX] = {0};
uint32_t UINTMAX = 0;

static uint32_t muddle_bits(uint16_t v) {
    unsigned int b = 0;
    uint32_t r = 0;
    for (b = 0; b < BITS; b ++) {
        r += ((v >> b) & 1u) << (2 * b);
    }
    return r;
}

static void init_table() {
    if(translate_table[1] != 0) return;
    UINTMAX = (1u << (BITS - 1));
    int i;
    for(i = 0; i < UINT16_MAX; i ++) {
        translate_table[i] = muddle_bits(i);
    }
}

/* 
 * translate x, y to quadindex 
 * 0 <= x < 2
 * 0 <= y < 2 
 * usually just do 0 <= x <= 1 is OK.
 *
 * */
quadindex_t quadindex(double x, double y) {
    init_table();
    uint32_t xi = x * UINTMAX;
    uint32_t yi = y * UINTMAX;
    uint32_t r = 0;
    r |= muddle_bits((uint16_t) xi);
    r |= muddle_bits((uint16_t) yi) << 1u;
#if BITS > 16
    r |= muddle_bits((uint16_t) (xi >> 16u)) << 32u;
    r |= muddle_bits((uint16_t) (yi >> 16u)) << 33u;
#endif
    return r;
}

#if 0
char * bin(quadindex_t qi) {
    char * rt = calloc(128, 1);
    unsigned int b;
    for(b = 0; b < BITS * 2 + 4; b ++) {
        rt[BITS * 2 + 3 - b] = ((qi >> b) & 1u) + '0';
    }
    return rt;
}
int main(int argc, char * argv) {
    double x = 0;
    double y = 0;
    for(x = 0; x <= 1.0; x += 0.25) {
        for(y = 0; y <= 1.0; y += 0.25) {
            printf("%10.5g %10.5g % 20s\n", x, y, bin(quadindex(x, y)));
        }
    }
    return 0;
}
#endif
