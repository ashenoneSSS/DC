#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <mpi.h>

using namespace std;

// ===== Глобальні =====
int ProcRank = 0;   // ранг процесу
int ProcNum  = 0;   // кількість процесів
const double InfinitiesPercent = 50.0;

// ===== Прототипи (щоб порядок не мав значення) =====
int  Min(int A, int B);
void ProcessInitialization(int *&pMatrix, int *&pProcRows, int &Size, int &RowNum);
void ProcessInitializationTest(int *&pMatrix, int *&pProcRows, int &Size, int &RowNum);
void ProcessTermination(int *pMatrix, int *pProcRows);
void DummyDataInitialization(int *pMatrix, int Size);
void RandomDataInitialization(int *pMatrix, int Size);
void DataDistribution(int *pMatrix, int *pProcRows, int Size, int RowNum);
void ResultCollection(int *pMatrix, int *pProcRows, int Size, int RowNum);
void ParallelFloyd(int *pProcRows, int Size, int RowNum);
void RowDistribution(int *pProcRows, int Size, int RowNum, int k, int *pRow);
void ParallelPrintMatrix(int *pProcRows, int Size, int RowNum);
void TestDistribution(int *pMatrix, int *pProcRows, int Size, int RowNum);
void TestResult(int *pMatrix, int *pSerialMatrix, int Size);

// Допоміжні для тестів/відладки
void PrintMatrix(int *M, int R, int C);
void CopyMatrix(int *src, int N, int *dst);
bool CompareMatrices(int *A, int *B, int N);
void SerialFloyd(int *pMatrix, int Size);

// ================= main =================
int main(int argc, char* argv[]) {
    int *pMatrix = nullptr;   // повна матриця на 0-му
    int Size = 0;             // розмір графа
    int *pProcRows = nullptr; // смуга рядків на кожному процесі
    int RowNum = 0;           // кількість рядків у процесі
    double start, finish, duration;
    int *pSerialMatrix = nullptr;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

    // 7 запусків: 10, 500, 600, 700, 800, 900, 1000
    for (int t = 0; t < 7; ++t) {
        Size = (t == 0) ? 10 : (500 + (t - 1) * 100);

        ProcessInitializationTest(pMatrix, pProcRows, Size, RowNum);

        if (ProcRank == 0) {
            pSerialMatrix = new int[Size * Size];
            CopyMatrix(pMatrix, Size, pSerialMatrix);
        }

        start = MPI_Wtime();

        DataDistribution(pMatrix, pProcRows, Size, RowNum);
        //TestDistribution(pMatrix, pProcRows, Size, RowNum);

        ParallelFloyd(pProcRows, Size, RowNum);

        ResultCollection(pMatrix, pProcRows, Size, RowNum);

        finish = MPI_Wtime();
        duration = finish - start;

        if (ProcRank == 0) {
            //TestResult(pMatrix, pSerialMatrix, Size);
            printf("Time of execution: %f\n", duration);
            delete [] pSerialMatrix;
            pSerialMatrix = nullptr;
        }

        ProcessTermination(pMatrix, pProcRows);
    }

    MPI_Finalize();
    return 0;
}

// ================= Реалізації =================
int Min(int A, int B) {
    int Result = (A < B) ? A : B;
    if ((A < 0) && (B >= 0)) Result = B;
    if ((B < 0) && (A >= 0)) Result = A;
    if ((A < 0) && (B < 0)) Result = -1;
    return Result;
}

void ProcessInitialization(int *&pMatrix, int *&pProcRows, int &Size, int &RowNum) {
    setvbuf(stdout, 0, _IONBF, 0);

    if (ProcRank == 0) {
        do {
            printf("Enter the number of vertices: ");
            if (scanf("%d", &Size) != 1) {
                int c; while ((c = getchar()) != '\n' && c != EOF) {}
                Size = 0;
            }
            if (Size < ProcNum)
                printf("The number of vertices should be greater then the number of processes\n");
        } while (Size < ProcNum);
        printf("Using the graph with %d vertices\n", Size);
    }

    MPI_Bcast(&Size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int RestRows = Size;
    for (int i = 0; i < ProcRank; ++i) RestRows -= RestRows / (ProcNum - i);
    RowNum = RestRows / (ProcNum - ProcRank);

    pProcRows = new int[Size * RowNum];

    if (ProcRank == 0) {
        pMatrix = new int[Size * Size];
        DummyDataInitialization(pMatrix, Size);
        // RandomDataInitialization(pMatrix, Size);
    }
}

void ProcessInitializationTest(int *&pMatrix, int *&pProcRows, int &Size, int &RowNum) {
    setvbuf(stdout, 0, _IONBF, 0);

    if (ProcRank == 0) {
        if (Size < ProcNum)
            printf("The number of vertices should be greater then the number of processes\n");
        printf("Using the graph with %d vertices\n", Size);
    }

    MPI_Bcast(&Size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int RestRows = Size;
    for (int i = 0; i < ProcRank; ++i) RestRows -= RestRows / (ProcNum - i);
    RowNum = RestRows / (ProcNum - ProcRank);

    pProcRows = new int[Size * RowNum];

    if (ProcRank == 0) {
        pMatrix = new int[Size * Size];
        DummyDataInitialization(pMatrix, Size);
        // RandomDataInitialization(pMatrix, Size);
    }
}

void ProcessTermination(int *pMatrix, int *pProcRows) {
    if (ProcRank == 0) delete [] pMatrix;
    delete [] pProcRows;
}

void DummyDataInitialization(int *pMatrix, int Size) {
    for (int i = 0; i < Size; ++i)
        for (int j = i; j < Size; ++j) {
            if (i == j) pMatrix[i*Size + j] = 0;
            else if (i == 0) pMatrix[i*Size + j] = j;
            else pMatrix[i*Size + j] = -1;
            pMatrix[j*Size + i] = pMatrix[i*Size + j];
        }
}

void RandomDataInitialization(int *pMatrix, int Size) {
    srand((unsigned)time(nullptr));
    for (int i = 0; i < Size; ++i)
        for (int j = 0; j < Size; ++j)
            if (i != j) {
                if ((rand() % 100) < InfinitiesPercent)
                    pMatrix[i*Size + j] = -1;
                else
                    pMatrix[i*Size + j] = rand() + 1;
            } else {
                pMatrix[i*Size + j] = 0;
            }
}

void DataDistribution(int *pMatrix, int *pProcRows, int Size, int RowNum) {
    int *pSendNum = new int[ProcNum];
    int *pSendInd = new int[ProcNum];

    int RestRows = Size;
    RowNum       = Size / ProcNum;
    pSendNum[0]  = RowNum * Size;
    pSendInd[0]  = 0;

    for (int i = 1; i < ProcNum; ++i) {
        RestRows   -= RowNum;
        RowNum      = RestRows / (ProcNum - i);
        pSendNum[i] = RowNum * Size;
        pSendInd[i] = pSendInd[i-1] + pSendNum[i-1];
    }

    MPI_Scatterv(pMatrix, pSendNum, pSendInd, MPI_INT,
                 pProcRows, pSendNum[ProcRank], MPI_INT,
                 0, MPI_COMM_WORLD);

    delete [] pSendNum;
    delete [] pSendInd;
}

void ResultCollection(int *pMatrix, int *pProcRows, int Size, int RowNum) {
    int *pReceiveNum = new int[ProcNum];
    int *pReceiveInd = new int[ProcNum];

    int RestRows = Size;
    RowNum = Size / ProcNum;
    pReceiveInd[0] = 0;
    pReceiveNum[0] = RowNum * Size;

    for (int i = 1; i < ProcNum; ++i) {
        RestRows      -= RowNum;
        RowNum         = RestRows / (ProcNum - i);
        pReceiveNum[i] = RowNum * Size;
        pReceiveInd[i] = pReceiveInd[i-1] + pReceiveNum[i-1];
    }

    MPI_Gatherv(pProcRows, pReceiveNum[ProcRank], MPI_INT,
                pMatrix, pReceiveNum, pReceiveInd, MPI_INT,
                0, MPI_COMM_WORLD);

    delete [] pReceiveNum;
    delete [] pReceiveInd;
}

void RowDistribution(int *pProcRows, int Size, int RowNum, int k, int *pRow) {
    int ProcRowRank, ProcRowNum;

    int RestRows = Size;
    int Ind = 0;
    int Num = Size / ProcNum;

    for (ProcRowRank = 1; ProcRowRank < ProcNum + 1; ++ProcRowRank) {
        if (k < Ind + Num) break;
        RestRows -= Num;
        Ind      += Num;
        Num       = RestRows / (ProcNum - ProcRowRank);
    }
    ProcRowRank -= 1;
    ProcRowNum   = k - Ind;

    if (ProcRowRank == ProcRank) {
        std::copy(&pProcRows[ProcRowNum*Size],
                  &pProcRows[(ProcRowNum+1)*Size],
                  pRow);
    }

    MPI_Bcast(pRow, Size, MPI_INT, ProcRowRank, MPI_COMM_WORLD);
}

void ParallelFloyd(int *pProcRows, int Size, int RowNum) {
    int *pRow = new int[Size];
    int t1, t2;

    for (int k = 0; k < Size; ++k) {
        RowDistribution(pProcRows, Size, RowNum, k, pRow);

        for (int i = 0; i < RowNum; ++i)
            for (int j = 0; j < Size; ++j)
                if ((pProcRows[i*Size + k] != -1) && (pRow[j] != -1)) {
                    t1 = pProcRows[i*Size + j];
                    t2 = pProcRows[i*Size + k] + pRow[j];
                    pProcRows[i*Size + j] = Min(t1, t2);
                }
    }
    delete [] pRow;
}

// ===== Відладочні/тести (не обов'язково викликати) =====
void ParallelPrintMatrix(int *pProcRows, int Size, int RowNum) {
    for (int r = 0; r < ProcNum; ++r) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (ProcRank == r) {
            printf("ProcRank = %d\n", ProcRank);
            printf("Proc rows:\n");
            PrintMatrix(pProcRows, RowNum, Size);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

void TestDistribution(int *pMatrix, int *pProcRows, int Size, int RowNum) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (ProcRank == 0) {
        printf("Initial adjacency matrix:\n");
        PrintMatrix(pMatrix, Size, Size);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    ParallelPrintMatrix(pProcRows, Size, RowNum);
}

void TestResult(int *pMatrix, int *pSerialMatrix, int Size) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (ProcRank == 0) {
        SerialFloyd(pSerialMatrix, Size);
        if (!CompareMatrices(pMatrix, pSerialMatrix, Size))
            printf("Results of serial and parallel algorithms are NOT identical. Check your code\n");
        else
            printf("Results of serial and parallel algorithms are identical\n");
    }
}

void PrintMatrix(int *M, int R, int C) {
    for (int i = 0; i < R; ++i) {
        for (int j = 0; j < C; ++j) printf("%6d ", M[i*C + j]);
        printf("\n");
    }
}

void CopyMatrix(int *src, int N, int *dst) {
    std::copy(src, src + N* (size_t)N, dst);
}

bool CompareMatrices(int *A, int *B, int N) {
    for (int i = 0; i < N*N; ++i)
        if (A[i] != B[i]) return false;
    return true;
}

// Послідовний Флойд для перевірки
void SerialFloyd(int *pMatrix, int Size) {
    int t1, t2;
    for (int k = 0; k < Size; ++k)
        for (int i = 0; i < Size; ++i)
            for (int j = 0; j < Size; ++j)
                if ((pMatrix[i*Size + k] != -1) &&
                    (pMatrix[k*Size + j] != -1)) {
                    t1 = pMatrix[i*Size + j];
                    t2 = pMatrix[i*Size + k] + pMatrix[k*Size + j];
                    pMatrix[i*Size + j] = Min(t1, t2);
                }
}
