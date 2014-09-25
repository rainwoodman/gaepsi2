#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <gaepsi.h>

#include <mpi.h>
#include <math.h>
#include <bigfile-mpi.h>

typedef struct {
    double pos[3];
    float mass;
    float value;
    float sml;
    int task;
} PStruct;

int NTask2D[2];
int ThisTask2D[2];
int NTask;
int ThisTask;

int ThisTaskOffsetPx[2];
int ThisTaskSizePx[2];
int64_t FullSizePx[2];

double BoxSize;
MPI_Comm Cart2D;
SVRemap svr = {0};
char * filename = "/physics2/yfeng1/BWSim/TEST/TEST-fof4/PART_027";
//char * filename = "/home/yfeng1/PART_032";
//char * filename = "fakedata";
char * imagefilename = "image";
int NumWriters = 16;
int remap[3][3] = {
    {1, 1, 0}, 
    {0, 1, 0},
    {0, 0, 1}};


int TilePadding = 256;
int PIXEL_WIDTH = 1600;
int DEBUG = 0;

PStruct * Pread = NULL;
PStruct * Precv = NULL;

int64_t NumPartRead = 0;
int64_t NumPartMax = 0;
int64_t NumPartTotal = 0;
int64_t NumPart = 0;
int64_t NumPartRecv;

struct {
    float * V;
    float * W;
} MyImg;

BigFile Snapshot = {0};
double (*sphkernel)(double r);

static void init() {

    sphkernel = gsph_spline_query(SPLINE_2D_PROJ_CUBIC);

    /* initialize the remap and geometry */
    svremap_init(&svr, remap);

    double aspect = svr.size[0] / svr.size[1];

    int per[2] = {0};
    
    if(big_file_mpi_open(&Snapshot, filename, MPI_COMM_WORLD)) {
        fprintf(stderr, big_file_get_error_message());
        abort();
    }

    /* initialize the Cartisan connection */

    MPI_Comm_size(MPI_COMM_WORLD, &NTask);
    MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);

    int foo = sqrt(NTask * aspect) + 1;
    while(NTask % foo) {
        foo --;
    }
    NTask2D[0] = foo;
    NTask2D[1] = NTask / foo;

    MPI_Cart_create(MPI_COMM_WORLD, 2, NTask2D, per, 0, &Cart2D);

    MPI_Cart_get(Cart2D, 2, NTask2D, per, ThisTask2D);

    FullSizePx[0] = aspect * PIXEL_WIDTH;
    FullSizePx[1] = PIXEL_WIDTH;

    /* partition the image */
    int i;
    for(i = 0; i < 2; i ++) {
        /* align to TilePadding pixel edges, but the last process will do the padding */
        ThisTaskOffsetPx[i] = 
            (ThisTask2D[i] * (int64_t) (FullSizePx[i] / TilePadding)) / NTask2D[i] * TilePadding;
        ThisTaskSizePx[i] = 
            ((ThisTask2D[i] + 1 ) * (int64_t) (FullSizePx[i] / TilePadding)) 
            / NTask2D[i] * TilePadding - ThisTaskOffsetPx[i];
        if(ThisTask2D[i] == NTask2D[i] - 1) {
            ThisTaskSizePx[i] = FullSizePx[i] - ThisTaskOffsetPx[i];
        }
    }
    int mysize = ThisTaskSizePx[0] * ThisTaskSizePx[1];
    MyImg.V = calloc(mysize, sizeof(float));
    MyImg.W = calloc(mysize, sizeof(float));

    if(ThisTask == 0) {
        printf("Box: %g %g %g\n", svr.size[0], svr.size[1], svr.size[2]);
        printf("Layout: %d x %d\n", NTask2D[0], NTask2D[1]);
        printf("Full Image : %ld x %ld \n", FullSizePx[0], FullSizePx[1]);
    }
    printf("ThisTask = %d My Image : %d x %d + %d x %d \n", 
            ThisTask,
            ThisTaskOffsetPx[0], ThisTaskOffsetPx[1], 
            ThisTaskSizePx[0], ThisTaskSizePx[1]);

}

static void read() {
    BigBlock header = {0};
    int64_t NumPartTotalAll[6];

    if(0 != big_file_mpi_open_block(&Snapshot, &header, "header", MPI_COMM_WORLD)) {
        fprintf(stderr, "%s", big_file_get_error_message());
        abort();
    }

    if(0 != big_block_get_attr(&header, "TotNumPart", NumPartTotalAll, "i8", 6) 
     && 0 != big_block_get_attr(&header, "NumPartInGroupTotal", NumPartTotalAll, "i8", 6)) {
        fprintf(stderr, "%s", big_file_get_error_message());
        abort();
    }

    if(0 != big_block_get_attr(&header, "BoxSize", &BoxSize, "f8", 1)) {
        fprintf(stderr, "%s", big_file_get_error_message());
        abort();
    }
    /* only do gas */
    NumPartTotal = NumPartTotalAll[0];
    int64_t OffPartRead = ThisTask * NumPartTotal / NTask;

    NumPartRead = (ThisTask + 1) * NumPartTotal / NTask
             - ThisTask * NumPartTotal / NTask;
    NumPartMax = NumPartTotal * 8 / NTask + 16;

    Pread = malloc(sizeof(PStruct) * NumPartMax);

    BigArray array = {0};

    BigBlock bb = {0};

    BigBlockPtr ptr = {0};

    size_t dims[2] = {NumPartRead, 3};
    ptrdiff_t strides[2] = {sizeof(PStruct), sizeof(double)};

    /* read Position */
    big_file_mpi_open_block(&Snapshot, &bb, "0/Position", MPI_COMM_WORLD);
    big_array_init(&array, &Pread[0].pos[0], "f8", 2, dims, strides);
    big_block_seek(&bb, &ptr, OffPartRead);
    big_block_read(&bb, &ptr, &array);
    big_block_mpi_close(&bb, MPI_COMM_WORLD);

    /* read Mass */
    dims[1] = 1;
    strides[1] = sizeof(float);

    big_file_mpi_open_block(&Snapshot, &bb, "0/Mass", MPI_COMM_WORLD);
    big_array_init(&array, &Pread[0].mass, "f4", 2, dims, strides);
    big_block_seek(&bb, &ptr, OffPartRead);
    big_block_read(&bb, &ptr, &array);
    big_block_mpi_close(&bb, MPI_COMM_WORLD);

    /* read Sml */
    dims[1] = 1;
    strides[1] = sizeof(float);

    big_file_mpi_open_block(&Snapshot, &bb, "0/SmoothingLength", MPI_COMM_WORLD);
    big_array_init(&array, &Pread[0].sml, "f4", 2, dims, strides);
    big_block_seek(&bb, &ptr, OffPartRead);
    big_block_read(&bb, &ptr, &array);
    big_block_mpi_close(&bb, MPI_COMM_WORLD);

    /* read Value */
    dims[1] = 1;
    strides[1] = sizeof(float);

    big_file_mpi_open_block(&Snapshot, &bb, "0/Entropy", MPI_COMM_WORLD);
    big_array_init(&array, &Pread[0].value, "f4", 2, dims, strides);
    big_block_seek(&bb, &ptr, OffPartRead);
    big_block_read(&bb, &ptr, &array);
    big_block_mpi_close(&bb, MPI_COMM_WORLD);
}
static int findtask(int x, int * offsets, int N) {
    int left = 0;
    int right = N;
    if(x >= offsets[N - 1]) return N - 1;
    if(x < offsets[0]) return -1;

    while(right > left + 1) {
        int mid = left + ((right - left) >> 1);
        if(offsets[mid] >= x) {
            right = mid;
        } else {
            left = mid;
        }
    }
    //printf("find task: offsets[%d] = %d, x = %d\n", left, offsets[left], x);
    return left;
}

static void remap_and_tag() {
    int * (TaskOffsetPx[2]);
    int k;
    for(k = 0; k < 2; k ++) {
        TaskOffsetPx[k] = malloc(sizeof(int) * NTask2D[k]);
        int remaindims[2] = {0};
        remaindims[k] = 1;
        MPI_Comm newcomm;
        MPI_Cart_sub(Cart2D, remaindims, &newcomm);
        MPI_Allgather(&ThisTaskOffsetPx[k], 1, MPI_INT, 
            TaskOffsetPx[k], 1, MPI_INT, newcomm);
        MPI_Comm_free(&newcomm);
        int i;

        /*
        printf("TaskOffsetPx ");
        for(i = 0; i < NTask2D[k]; i ++) {
            printf(" %d", TaskOffsetPx[k][i]);
        }
        printf("\n");
        fflush(stdout);
        */
    }

    /* tag particles for image tiles dup if needed */
    int64_t i;
    NumPart = NumPartRead;
    int tmpI[3] = {0};
    /* remap and transform to pixel coordinate */
    for(i = 0; i < NumPartRead; i ++) {
        int k;
        double tmp[3];
        for(k = 0; k < 3; k ++) {
            tmp[k] = Pread[i].pos[k] / BoxSize;
            while(tmp[k] < 0.0) tmp[k] += 1.0;
            while(tmp[k] >= 1.0) tmp[k] -= 1.0;
        }
        svremap_apply(&svr, tmp, &Pread[i].pos[0], tmpI);
        for(k = 0; k < 2; k ++) {
            Pread[i].pos[k] *= (FullSizePx[k] / svr.size[k]);
        }

        /* transform sml to pixel coordinate */
        Pread[i].sml /= BoxSize;
        Pread[i].sml *= (FullSizePx[0] / svr.size[0]);

        /* now find ThisTask2D*/
        int left[2];
        int right[2];

        for(k = 0; k < 2; k ++) {
            left[k] = findtask((int) floor(Pread[i].pos[k] - Pread[i].sml),
                    TaskOffsetPx[k], NTask2D[k]);
            right[k] = findtask((int) ceil(Pread[i].pos[k] + Pread[i].sml),
                    TaskOffsetPx[k], NTask2D[k]);
            if(left[k] < 0) left[k] = 0;
            if(left[k] >= NTask2D[k]) left[k] = NTask2D[k] - 1;
            if(right[k] < 0) right[k] = 0;
            if(right[k] >= NTask2D[k]) right[k] = NTask2D[k] - 1;
        }

        int task[2];
        int t = 0;
        for(task[0] = left[0]; task[0] <= right[0]; task[0] ++) {
        for(task[1] = left[1]; task[1] <= right[1]; task[1] ++) {
            int rank = 0;
            MPI_Cart_rank(Cart2D, task, &rank);
            
            if(t == 0) {
                Pread[i].task = rank;
            } else {
                Pread[NumPart] = Pread[i];
                Pread[NumPart].task = rank;
                NumPart ++;
                if(NumPart >= NumPartMax) {
                    fprintf(stderr, "out of storage Pread"); 
                    abort();
                }
            }
            t ++;
        }
        }
    }
}

static int cmptask(const void * v1, const void * v2) {
    int t1 = ((PStruct * )v1)->task;
    int t2 = ((PStruct * )v2)->task;
    return (t1 > t2) - (t1 < t2);
}

static void sort_and_exchange() {
    MPI_Datatype MPI_TYPE_PART;
    Precv = Pread + NumPart;

    /* first sort */
    qsort(Pread, NumPart, sizeof(Pread[0]), cmptask);
    /* now set up the communcation layout */
    MPI_Type_contiguous(sizeof(Pread[0]), MPI_BYTE, &MPI_TYPE_PART);
    MPI_Type_commit(&MPI_TYPE_PART);


    int SendCount[NTask];
    int RecvCount[NTask];
    int SendDispl[NTask];
    int RecvDispl[NTask];

    memset(SendCount, 0, sizeof(int) * NTask);
    int64_t i;

    for(i = 0; i < NumPart; i ++) {
        SendCount[Pread[i].task] ++;
    }

    MPI_Alltoall(SendCount, 1, MPI_INT,
            RecvCount, 1, MPI_INT, 
            MPI_COMM_WORLD);

    SendDispl[0] = 0;
    RecvDispl[0] = 0;

    NumPartRecv = RecvCount[0];
    for(i = 1; i < NTask; i ++) {
        SendDispl[i] = SendDispl[i - 1] + SendCount[i - 1];
        RecvDispl[i] = RecvDispl[i - 1] + RecvCount[i - 1];
        NumPartRecv += RecvCount[i];
    }

    if(NumPartRecv > NumPartMax - NumPart) {
        fprintf(stderr, "out of storage, need %g * NumPartTotal / NTask", 1.0 * NumPartRecv * NTask / NumPartTotal);
        abort();
    }

    MPI_Alltoallv(
            Pread, SendCount, SendDispl, MPI_TYPE_PART,
            Precv, RecvCount, RecvDispl, MPI_TYPE_PART,
            MPI_COMM_WORLD);

    printf("ThisTask = %d Recv %ld particles\n", ThisTask, NumPartRecv);

    MPI_Type_free(&MPI_TYPE_PART);
}

static void paint() {
    GSPHImage image = {0};
    GSPHPainter painter = {0};
    gsph_image_init(&image, 
        ThisTaskSizePx, "f4", NULL, MyImg.W, MyImg.V);
    
    gsph_painter_init(&painter,
            ThisTaskSizePx,
            gsph_image_write,
            sphkernel,
            &image);

    if(ThisTaskSizePx[0] == 0 || ThisTaskSizePx[1] == 0) return;
    int64_t i;

    for(i = 0; i < NumPartRecv; i ++) {
        double sml = Precv[i].sml;
        double mass = Precv[i].mass;
        double mvalue = mass * Precv[i].value;
        double * pos = Precv[i].pos;

        int k;
        for(k = 0; k < 2; k ++) {
            pos[k] -= ThisTaskOffsetPx[k];
        }

        gsph_painter_rasterize(&painter, pos, sml, mass, &mvalue, 1);
    }
    printf("ThisTask=%d, done painting\n", ThisTask);
}

static void write() {
    BigFile imagefile = {0};
    BigBlock header = {0};
    BigBlock weight = {0};
    BigBlock value = {0};

    if(0 != big_file_mpi_create(&imagefile, imagefilename, MPI_COMM_WORLD)) {
        fprintf(stderr, big_file_get_error_message());
        abort();
    }

    if(0 != big_file_mpi_create_block(&imagefile, &header, "header", NULL, 1, 
                0, 0, MPI_COMM_WORLD)) {
        fprintf(stderr, big_file_get_error_message());
        abort();
    }
    
    int NTile = 1;
    int NTile2D[2];
    int k;
    for(k = 0; k < 2; k ++) {
        NTile2D[k] = FullSizePx[k] / TilePadding + (FullSizePx[k] % TilePadding != 0);
        NTile *= NTile2D[k];
    }
    big_block_set_attr(&header, "ImageSize", FullSizePx, "i8", 2);
    big_block_set_attr(&header, "TilePadding", &TilePadding, "i4", 1);
    big_block_set_attr(&header, "NTile", NTile2D, "i4", 2);

    if(0 != big_file_mpi_create_block(&imagefile, &value, "V", "f4", 
                TilePadding * TilePadding,
                NumWriters, NTile, MPI_COMM_WORLD)) {
        fprintf(stderr, big_file_get_error_message());
        abort();
    }

    if(0 != big_file_mpi_create_block(&imagefile, &weight, "W", "f4", 
                TilePadding * TilePadding,
                NumWriters, NTile, MPI_COMM_WORLD)) {
        fprintf(stderr, big_file_get_error_message());
        abort();
    }
    if(DEBUG) {
        int i;
        for(i = 0; i < ThisTaskSizePx[0]; i ++) {
            MyImg.V[i * ThisTaskSizePx[1]] = -1.0;
            MyImg.W[i * ThisTaskSizePx[1]] = -1.0;
        }
        for(i = 0; i < ThisTaskSizePx[1]; i ++) {
            MyImg.V[i] = -1.0;
            MyImg.W[i] = -1.0;
        }
    }

    float bufferv[TilePadding * TilePadding];
    float bufferw[TilePadding * TilePadding];
    BigArray arrayv = {0};
    BigArray arrayw = {0};
    size_t dims[2] = {1, TilePadding * TilePadding};

    big_array_init(&arrayv, bufferv, "f4", 2, dims, NULL);
    big_array_init(&arrayw, bufferw, "f4", 2, dims, NULL);

    int it[2];
    int start[2], end[2];
    for(k = 0; k < 2; k ++) {
        start[k] = ThisTaskOffsetPx[k] / TilePadding;
        end[k] = start[k] + 
               ThisTaskSizePx[k] / TilePadding + (ThisTaskSizePx[k] % TilePadding != 0);

    }
    for(it[0] = start[0]; it[0] < end[0]; it[0]++) {
    for(it[1] = start[1]; it[1] < end[1]; it[1]++) {
        BigBlockPtr ptr = {0};
        int p = it[0] * NTile2D[1] + it[1];
        int x, y;
        ptrdiff_t offsetpx[2];
        for(k = 0; k < 2; k ++) {
            offsetpx[k] = (it[k] - start[k]) * TilePadding;
        }
        for(y = 0; y < TilePadding; y ++) {
            for(x = 0; x < TilePadding; x ++) {
                if( (y + offsetpx[0] >= ThisTaskSizePx[0])
                  ||(x + offsetpx[1] >= ThisTaskSizePx[1])) {
                    bufferv[y * TilePadding + x] = 0;
                    bufferw[y * TilePadding + x] = 0;
                } else{
                    ptrdiff_t pp = (offsetpx[0] + y) * ThisTaskSizePx[1]
                        + offsetpx[1] + x;
                    bufferv[y * TilePadding + x] = MyImg.V[pp];
                    bufferw[y * TilePadding + x] = MyImg.W[pp];
                }
            }
        }
        big_block_seek(&weight, &ptr, p);
        big_block_write(&weight, &ptr, &arrayw);

        big_block_seek(&value, &ptr, p);
        big_block_write(&value, &ptr, &arrayv);
    }
    }

    printf("ThisTask = %d Done writing\n", ThisTask);
    big_block_mpi_close(&value, MPI_COMM_WORLD);
    big_block_mpi_close(&weight, MPI_COMM_WORLD);
    big_block_mpi_close(&header, MPI_COMM_WORLD);
    big_file_mpi_close(&imagefile, MPI_COMM_WORLD);
}

int main(int argc, char * argv[]) {
    MPI_Init(&argc, &argv);

    init();

    read();
    remap_and_tag();
    sort_and_exchange();

    paint();

    write();

    MPI_Finalize();
    return 0;
}
