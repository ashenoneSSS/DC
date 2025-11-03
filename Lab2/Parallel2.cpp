#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <mpi.h>

// ======= общие штуки =======
int ProcNum = 0;      // сколько процессов
int ProcRank = 0;     // мой ранг

// простая инициализация: все единицы
void DummyDataInitialization(double* A, double* B, int N) {
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            A[i*N + j] = 1.0;
            B[i*N + j] = 1.0;
        }
}

// печать матрицы (если нужно для отладки)
void PrintMatrix(double* M, int R, int C) {
    for (int i = 0; i < R; ++i) {
        for (int j = 0; j < C; ++j) printf("%7.4f ", M[i*C + j]);
        printf("\n");
    }
}

// обычное тройное умножение C += A*B (полные матрицы)
void SerialMatMul(double* A, double* B, double* C, int N) {
    for (int i = 0; i < N; ++i)
        for (int k = 0; k < N; ++k) {
            double aik = A[i*N + k];
            for (int j = 0; j < N; ++j)
                C[i*N + j] += aik * B[k*N + j];
        }
}

// ======= Ветвь 1: квадратная сетка процессов (Fox как у тебя) =======
/* Используем твою логику Fox/«шашечки» почти без изменений,
   только обёрнуто в функции и аккуратно с проверками ввода. */

static int GridSize;            // sqrt(ProcNum)
static int GridCoords[2];       // координаты в сетке
static MPI_Comm GridComm;       // декартов коммуникатор
static MPI_Comm RowComm, ColComm;

void CreateGridCommunicators() {
    int dims[2] = { GridSize, GridSize };
    int periods[2] = { 0, 0 };
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &GridComm);
    MPI_Cart_coords(GridComm, ProcRank, 2, GridCoords);

    int sub[2];
    sub[0] = 0; sub[1] = 1; MPI_Cart_sub(GridComm, sub, &RowComm);
    sub[0] = 1; sub[1] = 0; MPI_Cart_sub(GridComm, sub, &ColComm);
}

void CheckerboardScatter(double* M, double* Mblk, int N, int BS) {
    double* rowbuf = new double[BS * N];
    if (GridCoords[1] == 0)
        MPI_Scatter(M, BS*N, MPI_DOUBLE, rowbuf, BS*N, MPI_DOUBLE, 0, ColComm);
    for (int i = 0; i < BS; ++i)
        MPI_Scatter(&rowbuf[i*N], BS, MPI_DOUBLE,
                    &Mblk[i*BS], BS, MPI_DOUBLE, 0, RowComm);
    delete[] rowbuf;
}

void ABlockBcast(int iter, double* Abuf, double* Ablk, int BS) {
    int pivot = (GridCoords[0] + iter) % GridSize;
    if (GridCoords[1] == pivot)
        for (int i = 0; i < BS*BS; ++i) Abuf[i] = Ablk[i];
    MPI_Bcast(Abuf, BS*BS, MPI_DOUBLE, pivot, RowComm);
}

void BBlockShift(double* Bblk, int BS) {
    MPI_Status st;
    int next = (GridCoords[0] + 1) % GridSize;
    int prev = (GridCoords[0] - 1 + GridSize) % GridSize;
    MPI_Sendrecv_replace(Bblk, BS*BS, MPI_DOUBLE, next, 0, prev, 0, ColComm, &st);
}

void BlockMulAdd(double* Ablk, double* Bblk, double* Cblk, int BS) {
    for (int i = 0; i < BS; ++i)
        for (int k = 0; k < BS; ++k) {
            double aik = Ablk[i*BS + k];
            for (int j = 0; j < BS; ++j)
                Cblk[i*BS + j] += aik * Bblk[k*BS + j];
        }
}

void GatherCheckerboard(double* C, double* Cblk, int N, int BS) {
    double* rowbuf = new double[N * BS];
    for (int i = 0; i < BS; ++i)
        MPI_Gather(&Cblk[i*BS], BS, MPI_DOUBLE, &rowbuf[i*N], BS, MPI_DOUBLE, 0, RowComm);
    if (GridCoords[1] == 0)
        MPI_Gather(rowbuf, BS*N, MPI_DOUBLE, C, BS*N, MPI_DOUBLE, 0, ColComm);
    delete[] rowbuf;
}

void RunSquareGrid(int N) {
    int BS = N / GridSize;

    double *A = nullptr, *B = nullptr, *C = nullptr;
    double *Ablk = new double[BS*BS];
    double *Bblk = new double[BS*BS];
    double *Cblk = new double[BS*BS]();
    double *Abuf = new double[BS*BS];

    if (ProcRank == 0) {
        A = new double[N*N];
        B = new double[N*N];
        C = new double[N*N]();
        DummyDataInitialization(A, B, N);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double t0 = MPI_Wtime();

    // раздать блоки
    CheckerboardScatter(A, Ablk, N, BS);
    CheckerboardScatter(B, Bblk, N, BS);

    // Fox
    for (int it = 0; it < GridSize; ++it) {
        ABlockBcast(it, Abuf, Ablk, BS);
        BlockMulAdd(Abuf, Bblk, Cblk, BS);
        BBlockShift(Bblk, BS);
    }

    // собрать результат
    GatherCheckerboard(C, Cblk, N, BS);

    MPI_Barrier(MPI_COMM_WORLD);
    double t1 = MPI_Wtime();

    if (ProcRank == 0) {
        // проверка (по желанию — можно отключить)
        int ok = 1;
        double* Cref = (double*)calloc(N*N, sizeof(double));
        SerialMatMul(A, B, Cref, N);
        for (int i = 0; i < N*N; ++i) if (fabs(C[i]-Cref[i]) > 1e-6) { ok = 0; break; }
        printf("The results of serial and parallel algorithms are %s.\n", ok ? "identical" : "NOT identical");
        printf("Time of execution = %f\n", t1 - t0);
        delete[] A; delete[] B; delete[] C; free(Cref);
    }
    delete[] Ablk; delete[] Bblk; delete[] Cblk; delete[] Abuf;
}

// ======= Ветвь 2: НЕ квадрат процессов — строчно-разделённый вариант =======
void RunRowStriped(int N) {
    // ранг 0 создаёт полные A,B; всем нужна только B (рассылаем broadcast)
    double *A = nullptr, *B = nullptr, *C = nullptr;
    if (ProcRank == 0) {
        A = new double[N*N];
        B = new double[N*N];
        C = new double[N*N]();
        DummyDataInitialization(A, B, N);
    } else {
        B = new double[N*N]; // на остальных храним копию B
    }

    // раздать строки A: балансируем counts/displs
    int base = N / ProcNum, rem = N % ProcNum;
    int my_rows = base + (ProcRank < rem ? 1 : 0);

    int* counts = nullptr; 
    int* displs = nullptr;
    if (ProcRank == 0) {
        counts = new int[ProcNum];
        displs = new int[ProcNum];
        int off = 0;
        for (int r = 0; r < ProcNum; ++r) {
            int rows = base + (r < rem ? 1 : 0);
            counts[r] = rows * N;      // элементов double
            displs[r] = off * N;       // смещение в элементах
            off += rows;
        }
    }

    double* A_local = new double[(size_t)my_rows * N];
    double* C_local = new double[(size_t)my_rows * N]();

    // рассылаем строки A
    MPI_Scatterv(A, counts, displs, MPI_DOUBLE,
                 A_local, my_rows * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // рассылаем B всем
    MPI_Bcast(B, N*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    double t0 = MPI_Wtime();

    // локальное умножение (мои строки A * полная B -> мои строки C)
    for (int i = 0; i < my_rows; ++i)
        for (int k = 0; k < N; ++k) {
            double aik = A_local[i*N + k];
            for (int j = 0; j < N; ++j)
                C_local[i*N + j] += aik * B[k*N + j];
        }

    MPI_Barrier(MPI_COMM_WORLD);
    double t1 = MPI_Wtime();

    // собираем C
    MPI_Gatherv(C_local, my_rows * N, MPI_DOUBLE,
                C, counts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (ProcRank == 0) {
        // проверка (по желанию — можно отключить)
        int ok = 1;
        double* Cref = (double*)calloc(N*N, sizeof(double));
        SerialMatMul(A, B, Cref, N);
        for (int i = 0; i < N*N; ++i) if (fabs(C[i]-Cref[i]) > 1e-6) { ok = 0; break; }
        printf("The results of serial and parallel algorithms are %s.\n", ok ? "identical" : "NOT identical");
        printf("Time of execution = %f\n", t1 - t0);
        delete[] counts; delete[] displs;
    }

    delete[] A_local; delete[] C_local;
    if (ProcRank == 0) { delete[] A; delete[] B; delete[] C; }
    else { delete[] B; }
}

int main(int argc, char** argv) {
    setvbuf(stdout, 0, _IONBF, 0);
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

    if (ProcRank == 0) printf("Parallel matrix multiplication program\n");

    // читаем N на корне с проверкой
    int N = 0;
    if (ProcRank == 0) {
        printf("\nEnter the size of matrices: ");
        if (scanf("%d", &N) != 1) { fprintf(stderr, "Input error\n"); N = 0; }
    }
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (N <= 0) { MPI_Finalize(); return 0; }

    // если квадрат процессов — запускаем Fox; иначе — строчный вариант
    int g = (int)(sqrt((double)ProcNum) + 0.5);
    if (g * g == ProcNum) {
    GridSize = g;

    // если N не делится на GridSize, просто переходим в RowStriped
    if (N % GridSize != 0) {
        if (ProcRank == 0)
            printf("Matrix size %d not divisible by sqrt(np)=%d → using row-striped version\n", N, GridSize);
        RunRowStriped(N);
    } else {
        CreateGridCommunicators();
        RunSquareGrid(N);
    }
    } else {
        RunRowStriped(N);
    }


    MPI_Finalize();
    return 0;
}
