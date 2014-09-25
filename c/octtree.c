#include <stdio.h>
#include <stdlib.h>

typedef int64_t key_t;
typedef struct Array Array;

typedef key_t (*keyfunc_t)(Array * array, ptrdiff_t ind);
static inline key_t ikey(int ipos[3]);

struct Array {
    void * data;
    size_t shape[2];
    ptrdiff_t strides[2];
    size_t elsize;
};

typedef struct {
    union {
        int c[8]; 
        struct {
            int parent;
            int next;
            int sibling;
            int start;
            int size;
        };
    };
} Treenode;

struct Tree {
    TreeNode * nodes;
    int used;
    int size;
    ptrdiff_t * ind;
    double shift[3];
    double inv[3];
    keyfunc_t keyfunc;
};

key_t keyfunc_double(Array * array, ptrdiff_t ind) {
    int k;
    int ipos[3];
    for(k = 0; k < 3; k ++) {
        double * ptr = (double*) (((char*) array->data) + ind * array->strides[0] + k * array->strides[1]);
        double foo = *ptr;
        ipos[k] = ((1 << 21) - 1) * (foo - array->shift[k]) * inv[k];
        if(ipos[k] < 0) ipos[k] = 0;
        if(ipos[k] > (1 << 21) - 1) ipos[k] = (1 << 21) - 1;
    }
    return ikey(ipos);
}
void array_init(Array * array, void * data, char * dtype, ptrdiff_t shape[], ptrdiff_t strides[]) {
    if(!strstr(dtype, "f8")) {
        array->elsize = 8;
    } else {
        array->elsize = 4;
    }

    array->data = data;
    array->shape[0] = shape[0];
    array->shape[1] = shape[1];
    array->strides[0] = strides[0];
    array->strides[1] = strides[1];
}

void tree_init(Tree * tree, Array * array) {
    tree->array = array;
    tree->keyfunc = keyfunc_double;
    tree->ind = malloc(sizeof(ptrdiff_t) * array->shape[0]);
    struct sortitem {
        key_t key;
        ptrdiff_t ind;
    } * sortitems = malloc(sizeof(struct sortitem) * array->shape[0]);

    ptrdiff_t i;
    for(i = 0; i < array->shape[0]; i ++) {
        sortitems[i].key = array->keyfunc(array, i);
    }
    qsort(sortitems, array->shape[0], sizeof(struct sortitem), compare_sortitem_by_key);
    for(i = 0; i < array->shape[0]; i ++) {
        tree->ind[i] = sortitem[i].ind;
    }
    for(i = 0; i < array->shape[0]; i ++) {
        sortitem[k].key - sortitem[
    }
    

    free(sortitems);
}

int main(int argc, char * argv[]) {
    double test[1000][3];
    int i, j, k, u;
    u = 0;
    for(i = 0; i < 10; i ++)
    for(j = 0; j < 10; j ++)
    for(k = 0; k < 10; k ++) {
        test[i * 100 + j * 10 + k][0] = i;
        test[i * 100 + j * 10 + k][1] = j;
        test[i * 100 + j * 10 + k][2] = k;
    }
    Array array;
    ptrdiff_t strides[2] = {24, 8};
    ptrdiff_t shape[2] = {1000, 3};
    array_init(&array, &test[0][0], "f8", shape, strides);

    Tree tree;
    tree_init(&tree, array);
    return 0;
}

static ipos_t utab[256] = {
    0x000000, 0x000001, 0x000008, 0x000009,
    0x000040, 0x000041, 0x000048, 0x000049,
    0x000200, 0x000201, 0x000208, 0x000209,
    0x000240, 0x000241, 0x000248, 0x000249,
    0x001000, 0x001001, 0x001008, 0x001009,
    0x001040, 0x001041, 0x001048, 0x001049,
    0x001200, 0x001201, 0x001208, 0x001209,
    0x001240, 0x001241, 0x001248, 0x001249,
    0x008000, 0x008001, 0x008008, 0x008009,
    0x008040, 0x008041, 0x008048, 0x008049,
    0x008200, 0x008201, 0x008208, 0x008209,
    0x008240, 0x008241, 0x008248, 0x008249,
    0x009000, 0x009001, 0x009008, 0x009009,
    0x009040, 0x009041, 0x009048, 0x009049,
    0x009200, 0x009201, 0x009208, 0x009209,
    0x009240, 0x009241, 0x009248, 0x009249,
    0x040000, 0x040001, 0x040008, 0x040009,
    0x040040, 0x040041, 0x040048, 0x040049,
    0x040200, 0x040201, 0x040208, 0x040209,
    0x040240, 0x040241, 0x040248, 0x040249,
    0x041000, 0x041001, 0x041008, 0x041009,
    0x041040, 0x041041, 0x041048, 0x041049,
    0x041200, 0x041201, 0x041208, 0x041209,
    0x041240, 0x041241, 0x041248, 0x041249,
    0x048000, 0x048001, 0x048008, 0x048009,
    0x048040, 0x048041, 0x048048, 0x048049,
    0x048200, 0x048201, 0x048208, 0x048209,
    0x048240, 0x048241, 0x048248, 0x048249,
    0x049000, 0x049001, 0x049008, 0x049009,
    0x049040, 0x049041, 0x049048, 0x049049,
    0x049200, 0x049201, 0x049208, 0x049209,
    0x049240, 0x049241, 0x049248, 0x049249,
    0x200000, 0x200001, 0x200008, 0x200009,
    0x200040, 0x200041, 0x200048, 0x200049,
    0x200200, 0x200201, 0x200208, 0x200209,
    0x200240, 0x200241, 0x200248, 0x200249,
    0x201000, 0x201001, 0x201008, 0x201009,
    0x201040, 0x201041, 0x201048, 0x201049,
    0x201200, 0x201201, 0x201208, 0x201209,
    0x201240, 0x201241, 0x201248, 0x201249,
    0x208000, 0x208001, 0x208008, 0x208009,
    0x208040, 0x208041, 0x208048, 0x208049,
    0x208200, 0x208201, 0x208208, 0x208209,
    0x208240, 0x208241, 0x208248, 0x208249,
    0x209000, 0x209001, 0x209008, 0x209009,
    0x209040, 0x209041, 0x209048, 0x209049,
    0x209200, 0x209201, 0x209208, 0x209209,
    0x209240, 0x209241, 0x209248, 0x209249,
    0x240000, 0x240001, 0x240008, 0x240009,
    0x240040, 0x240041, 0x240048, 0x240049,
    0x240200, 0x240201, 0x240208, 0x240209,
    0x240240, 0x240241, 0x240248, 0x240249,
    0x241000, 0x241001, 0x241008, 0x241009,
    0x241040, 0x241041, 0x241048, 0x241049,
    0x241200, 0x241201, 0x241208, 0x241209,
    0x241240, 0x241241, 0x241248, 0x241249,
    0x248000, 0x248001, 0x248008, 0x248009,
    0x248040, 0x248041, 0x248048, 0x248049,
    0x248200, 0x248201, 0x248208, 0x248209,
    0x248240, 0x248241, 0x248248, 0x248249,
    0x249000, 0x249001, 0x249008, 0x249009,
    0x249040, 0x249041, 0x249048, 0x249049,
    0x249200, 0x249201, 0x249208, 0x249209,
    0x249240, 0x249241, 0x249248, 0x249249,
};
static inline key_t ikey(int ipos[3]) {
    int x = ipos[0];
    int y = ipos[1];
    int z = ipos[2];
    key_t ind = 0;
    x &= (1L << BITS) - 1;
    y &= (1L << BITS) - 1;
    z &= (1L << BITS) - 1;
    int i;
    for (i = 0; i < BITS; i+= 8) {
      ind |= (key_t) utab[(uint8_t) x] << (i*3+0);
      ind |= (key_t) utab[(uint8_t) y] << (i*3+1);
      ind |= (key_t) utab[(uint8_t) z] << (i*3+2);
      x >>= 8;
      y >>= 8;
      z >>= 8;
    }
    return ind;
}

