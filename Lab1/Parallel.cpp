#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

int ProcNum = 0;   // number of processes
int ProcRank = 0;  // rank of current process

// random init
void RandomDataInitialization(double* pMatrix, double* pVector, int Size) {
    srand((unsigned)time(NULL));
    for (int i = 0; i < Size; i++) {
        pVector[i] = rand() / 1000.0;
        for (int j = 0; j < Size; j++)
            pMatrix[i * Size + j] = rand() / 1000.0;
    }
}

// alloc + input with validation
void ProcessInitialization(double*& pMatrix, double*& pVector,
                           double*& pResult, double*& pProcRows,
                           double*& pProcResult, int& Size, int& RowNum) {
    int RestRows;
    setvbuf(stdout, 0, _IONBF, 0);

    if (ProcRank == 0) {
        int ok = 0;
        do {
            printf("\nEnter the size of the matrix and vector: ");
            int read = scanf("%d", &Size);

            if (read != 1) {                 // ввели не число
                int ch;
                while ((ch = getchar()) != '\n' && ch != EOF) {}
                printf("Please enter an integer value.\n");
                continue;
            }
            if (Size < ProcNum) {
                printf("Size must be >= number of processes (%d).\n", ProcNum);
                continue;
            }
            if (Size <= 0) {
                printf("Size must be positive.\n");
                continue;
            }
            ok = 1;
        } while (!ok);
    }

    // раздать Size всем
    MPI_Bcast(&Size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // сколько строк у этого процесса
    RestRows = Size;
    for (int i = 0; i < ProcRank; i++)
        RestRows = RestRows - RestRows / (ProcNum - i);
    RowNum = RestRows / (ProcNum - ProcRank);

    // память
    pVector      = new double[Size];
    pResult      = new double[Size];
    pProcRows    = new double[RowNum * Size];
    pProcResult  = new double[RowNum];

    if (ProcRank == 0) {
        pMatrix = new double[Size * Size];
        RandomDataInitialization(pMatrix, pVector, Size);
    }
}

// раздача данных
void DataDistribution(double* pMatrix, double* pProcRows, double* pVector,
                      int Size, int RowNum) {
    MPI_Bcast(pVector, Size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    int* pSendNum = new int[ProcNum];
    int* pSendInd = new int[ProcNum];

    int RestRows = Size;
    int chunk    = Size / ProcNum;
    pSendNum[0]  = chunk * Size;
    pSendInd[0]  = 0;

    for (int i = 1; i < ProcNum; i++) {
        RestRows -= chunk;
        chunk     = RestRows / (ProcNum - i);
        pSendNum[i] = chunk * Size;
        pSendInd[i] = pSendInd[i - 1] + pSendNum[i - 1];
    }

    MPI_Scatterv(pMatrix, pSendNum, pSendInd, MPI_DOUBLE,
                 pProcRows, pSendNum[ProcRank], MPI_DOUBLE,
                 0, MPI_COMM_WORLD);

    delete[] pSendNum;
    delete[] pSendInd;
}

// локное умножение
void ParallelResultCalculation(double* pProcRows, double* pVector,
                               double* pProcResult, int Size, int RowNum) {
    for (int i = 0; i < RowNum; i++) {
        double acc = 0.0;
        for (int j = 0; j < Size; j++)
            acc += pProcRows[i * Size + j] * pVector[j];
        pProcResult[i] = acc;
    }
}

// собрать результат на всех
void ResultReplication(double* pProcResult, double* pResult,
                       int Size, int RowNum) {
    int* pReceiveNum = new int[ProcNum];
    int* pReceiveInd = new int[ProcNum];

    int RestRows = Size;
    pReceiveInd[0] = 0;
    pReceiveNum[0] = Size / ProcNum;
    for (int i = 1; i < ProcNum; i++) {
        RestRows      -= pReceiveNum[i - 1];
        pReceiveNum[i] = RestRows / (ProcNum - i);
        pReceiveInd[i] = pReceiveInd[i - 1] + pReceiveNum[i - 1];
    }

    MPI_Allgatherv(pProcResult, pReceiveNum[ProcRank], MPI_DOUBLE,
                   pResult, pReceiveNum, pReceiveInd, MPI_DOUBLE,
                   MPI_COMM_WORLD);

    delete[] pReceiveNum;
    delete[] pReceiveInd;
}

// последовательная проверка (только на 0-м)
void SerialResultCalculation(double* pMatrix, double* pVector,
                             double* pResult, int Size) {
    for (int i = 0; i < Size; i++) {
        double acc = 0.0;
        for (int j = 0; j < Size; j++)
            acc += pMatrix[i * Size + j] * pVector[j];
        pResult[i] = acc;
    }
}

void TestResult(double* pMatrix, double* pVector, double* pResult, int Size) {
    if (ProcRank == 0) {
        double* pSerialResult = new double[Size];
        SerialResultCalculation(pMatrix, pVector, pSerialResult, Size);
        int ok = 1;
        for (int i = 0; i < Size; i++) {
            if (pResult[i] != pSerialResult[i]) { ok = 0; break; }
        }
        printf("Serial vs Parallel: %s\n", ok ? "OK" : "DIFF");
        delete[] pSerialResult;
    }
}

void ProcessTermination(double* pMatrix, double* pVector, double* pResult,
                        double* pProcRows, double* pProcResult) {
    if (ProcRank == 0) delete[] pMatrix;
    delete[] pVector;
    delete[] pResult;
    delete[] pProcRows;
    delete[] pProcResult;
}

int main(int argc, char* argv[]) {
    double *pMatrix = nullptr, *pVector = nullptr, *pResult = nullptr;
    double *pProcRows = nullptr, *pProcResult = nullptr;
    int Size = 0, RowNum = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

    if (ProcRank == 0)
        printf("Parallel matrix-vector multiplication program\n");

    ProcessInitialization(pMatrix, pVector, pResult, pProcRows, pProcResult,
                          Size, RowNum);

    double Start = MPI_Wtime();

    DataDistribution(pMatrix, pProcRows, pVector, Size, RowNum);
    ParallelResultCalculation(pProcRows, pVector, pProcResult, Size, RowNum);
    ResultReplication(pProcResult, pResult, Size, RowNum);

    double Finish = MPI_Wtime();
    double Duration = Finish - Start;

    TestResult(pMatrix, pVector, pResult, Size);

    if (ProcRank == 0)
        printf("\nTime of execution = %f\n", Duration);

    ProcessTermination(pMatrix, pVector, pResult, pProcRows, pProcResult);
    MPI_Finalize();
    return 0;
}
