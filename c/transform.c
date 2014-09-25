#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "gaepsi.h"

void gtmat_apply(GTMatrix T, double v1[4], double v2[4]) {
    int i, j;
    for(i = 0; i < 4; i ++) {
        v2[i] = 0;
        for(j = 0; j < 4; j ++) {
            v2[i] += v1[j] * T[i][j];
        }
    }
}

void gtmat_mul(GTMatrix T, GTMatrix A, GTMatrix B) {
    int i, j, k;
    for(i = 0; i < 4; i ++) {
        for(j = 0; j < 4; j ++) {
            T[i][j] = 0;
            for(k = 0; k < 4; k ++) {
                T[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

void gtmat_copy(GTMatrix T, GTMatrix A) {
    memcpy(T, A, sizeof(double) * 16);
}
void gtmat_inverse(GTMatrix T, GTMatrix A) {
    int i, j;
    for(i = 0; i < 4; i ++) {
    for(j = 0; j < 4; j ++) {
        T[i][j] = A[j][i];
    }
    }
}

char * gtmat_format(GTMatrix T) {
    int width = 10;
    char * rt = calloc(40 * width * 4 + 1, 1);
    memset(rt, ' ', 40 * width * 4);
    int i, j;
    for(i = 0; i < 4; i ++) {
        for(j = 0; j < 4; j ++) {
            char buf[128];
            sprintf(buf, "%6.2e", T[i][j]);
            memcpy(&rt[i * width * 4 + j * width], buf, strlen(buf));
        }
        rt[i * width * 4 + width * 4 - 1] = (i == 3)?0:'\n';
    }
    return rt;
}

#if 0
int main(int argc, char * argv[]) {
    GTMatrix r = {
        1, 0, 0, 0,
        0, 2, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1,
    };
    GTMatrix s = GTMAT_SHIFT(1, 1, 1);
    GTMatrix t;

    gtmat_mul(t, r, s);
    double v1[4] = {0, 0, 0, 1};
    double v2[4] = {1, 1, 1, 0};
    gtmat_apply(t, v1, v2);
    printf("s %s\n", gtmat_format(s));
    printf("r %s\n", gtmat_format(r));
    printf("t %s\n", gtmat_format(t));
    printf("%g %g %g %g -> %g %g %g %g\n", 
            v1[0], v1[1], v1[2], v1[3],
            v2[0], v2[1], v2[2], v2[3]);
    return 0;
}
#endif
