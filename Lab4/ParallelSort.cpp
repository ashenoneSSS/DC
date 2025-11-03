#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <mpi.h>

using namespace std;

const double RandomDataMultiplier = 1000.0;

int ProcNum = 0;   // число процессов
int ProcRank = -1; // ранг процесса

// ---------- прототипы наших функций ----------
void ProcessInitialization(double *&pData, int& DataSize,
                           double *&pProcData, int& BlockSize);
void ProcessTermination(double *pData, double *pProcData);
void DummyDataInitialization(double*& pData, int& DataSize);
void RandomDataInitialization(double *&pData, int& DataSize);
void DataDistribution(double *pData, int DataSize,
                      double *pProcData, int BlockSize);
void DataCollection(double *pData, int DataSize,
                    double *pProcData, int BlockSize);
void ParallelBubble(double *pProcData, int BlockSize);
void ExchangeData(double *pProcData, int BlockSize, int DualRank,
                  double *pDualData, int DualBlockSize);
void TestDistribution(double *pData, int DataSize,
                      double *pProcData, int BlockSize);
void ParallelPrintData(double *pProcData, int BlockSize);
void TestResult(double *pData, double *pSerialData, int DataSize);

// вспомогательные функции (то, что раньше было в *Test.h/*.cpp)
void CopyData(const double *src, int n, double *dst);
void PrintData(const double *data, int n);
void SerialBubbleSort(double *data, int n);
void SerialStdSort(double *data, int n);
bool CompareData(const double *a, const double *b, int n);

// ---------- main ----------
int main(int argc, char *argv[]) {
    double *pData       = nullptr;
    double *pProcData   = nullptr;
    int     DataSize    = 0;
    int     BlockSize   = 0;
    double *pSerialData = nullptr;
    double  start, finish;
    double  duration = 0.0;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

    if (ProcRank == 0)
        printf("Parallel bubble sort program\n");

    // Инициализация
    ProcessInitialization(pData, DataSize, pProcData, BlockSize);

    if (ProcRank == 0) {
        pSerialData = new double[DataSize];
        CopyData(pData, DataSize, pSerialData);
    }

    start = MPI_Wtime();

    // Распределение данных
    DataDistribution(pData, DataSize, pProcData, BlockSize);

    // Parallel bubble sort
    ParallelBubble(pProcData, BlockSize);

    // Сбор результатов
    DataCollection(pData, DataSize, pProcData, BlockSize);
    // TestResult(pData, pSerialData, DataSize);

    finish = MPI_Wtime();
    duration = finish - start;

    if (ProcRank == 0)
        printf("Time of execution: %f\n", duration);

    if (ProcRank == 0)
        delete [] pSerialData;

    ProcessTermination(pData, pProcData);
    MPI_Finalize();

    return 0;
}

// ---------- реализация функций ----------

// Выделение памяти и задание начальных значений
void ProcessInitialization(double *&pData, int& DataSize,
                           double *&pProcData, int& BlockSize) {
    // отключаем буферизацию stdout, чтобы printf сразу показывались
    setvbuf(stdout, nullptr, _IONBF, 0);

    if (ProcRank == 0) {
        while (true) {
            printf("Enter the size of data to be sorted: ");
            fflush(stdout);

            int rc = scanf("%d", &DataSize);
            if (rc != 1) {
                printf("Input error: please enter an integer value\n");
                int ch;
                while ((ch = getchar()) != '\n' && ch != EOF) { }
                DataSize = 0;
                continue;
            }

            if (DataSize < ProcNum) {
                printf("Data size should be greater than or equal to number of processes\n");
                continue;
            }
            break;
        }

        printf("Sorting %d data items\n", DataSize);
    }

    // Рассылаем размер всем процессам
    MPI_Bcast(&DataSize, 1, MPI_INT, 0, MPI_COMM_WORLD);


    int RestData = DataSize;
    for (int i = 0; i < ProcRank; ++i)
        RestData -= RestData / (ProcNum - i);
    BlockSize = RestData / (ProcNum - ProcRank);

    pProcData = new double[BlockSize];

    if (ProcRank == 0) {
        pData = new double[DataSize];
        // RandomDataInitialization(pData, DataSize);
        DummyDataInitialization(pData, DataSize);
    }
}

// Завершение вычислений
void ProcessTermination(double *pData, double *pProcData) {
    if (ProcRank == 0) {
        delete [] pData;
    }
    delete [] pProcData;
}

// Простая инициализация данных
void DummyDataInitialization(double*& pData, int& DataSize) {
    for (int i = 0; i < DataSize; ++i)
        pData[i] = DataSize - i;
}

// Инициализация случайными данными
void RandomDataInitialization(double *&pData, int& DataSize) {
    srand((unsigned)time(nullptr));
    for (int i = 0; i < DataSize; ++i)
        pData[i] = double(rand()) / RAND_MAX * RandomDataMultiplier;
}

// Распределение данных между процессами
void DataDistribution(double *pData, int DataSize,
                      double *pProcData, int BlockSize) {
    (void)BlockSize; // на случай придирчивого компилятора

    int *pSendInd = new int[ProcNum];
    int *pSendNum = new int[ProcNum];

    int RestData    = DataSize;
    int CurrentSize = DataSize / ProcNum;

    pSendNum[0] = CurrentSize;
    pSendInd[0] = 0;

    for (int i = 1; i < ProcNum; ++i) {
        RestData   -= CurrentSize;
        CurrentSize = RestData / (ProcNum - i);
        pSendNum[i] = CurrentSize;
        pSendInd[i] = pSendInd[i - 1] + pSendNum[i - 1];
    }

    MPI_Scatterv(pData, pSendNum, pSendInd, MPI_DOUBLE,
                 pProcData, pSendNum[ProcRank], MPI_DOUBLE,
                 0, MPI_COMM_WORLD);

    delete [] pSendNum;
    delete [] pSendInd;
}

// Сбор данных на нулевом процессе
void DataCollection(double *pData, int DataSize,
                    double *pProcData, int BlockSize) {
    (void)BlockSize;

    int *pReceiveNum = new int[ProcNum];
    int *pReceiveInd = new int[ProcNum];

    int RestData = DataSize;
    pReceiveInd[0] = 0;
    pReceiveNum[0] = DataSize / ProcNum;

    for (int i = 1; i < ProcNum; ++i) {
        RestData        -= pReceiveNum[i - 1];
        pReceiveNum[i]   = RestData / (ProcNum - i);
        pReceiveInd[i]   = pReceiveInd[i - 1] + pReceiveNum[i - 1];
    }

    MPI_Gatherv(pProcData, pReceiveNum[ProcRank], MPI_DOUBLE,
                pData, pReceiveNum, pReceiveInd, MPI_DOUBLE,
                0, MPI_COMM_WORLD);

    delete [] pReceiveNum;
    delete [] pReceiveInd;
}

enum split_mode { KeepFirstHalf, KeepSecondHalf };

// Параллельный bubble sort (чередование обменов блоками)
void ParallelBubble(double *pProcData, int BlockSize) {
    // Локальная сортировка
    // SerialBubbleSort(pProcData, BlockSize);
    SerialStdSort(pProcData, BlockSize);

    int        Offset;
    split_mode SplitMode;

    for (int i = 0; i < ProcNum; ++i) {
        if ((i % 2) == 1) {
            if ((ProcRank % 2) == 1) {
                Offset    = 1;
                SplitMode = KeepFirstHalf;
            } else {
                Offset    = -1;
                SplitMode = KeepSecondHalf;
            }
        } else {
            if ((ProcRank % 2) == 1) {
                Offset    = -1;
                SplitMode = KeepSecondHalf;
            } else {
                Offset    = 1;
                SplitMode = KeepFirstHalf;
            }
        }

        // границы
        if ((ProcRank == ProcNum - 1) && (Offset == 1))  continue;
        if ((ProcRank == 0)           && (Offset == -1)) continue;

        int DualBlockSize = 0;

        MPI_Sendrecv(&BlockSize, 1, MPI_INT, ProcRank + Offset, 0,
                     &DualBlockSize, 1, MPI_INT, ProcRank + Offset, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        double *pDualData   = new double[DualBlockSize];
        double *pMergedData = new double[BlockSize + DualBlockSize];

        // обмен
        ExchangeData(pProcData, BlockSize, ProcRank + Offset,
                     pDualData, DualBlockSize);


        // слияние
        std::merge(pProcData, pProcData + BlockSize,
                   pDualData,  pDualData + DualBlockSize,
                   pMergedData);

        // разделение обратно
        if (SplitMode == KeepFirstHalf) {
            std::copy(pMergedData,
                      pMergedData + BlockSize,
                      pProcData);
        } else {
            std::copy(pMergedData + BlockSize,
                      pMergedData + BlockSize + DualBlockSize,
                      pProcData);
        }

        delete [] pDualData;
        delete [] pMergedData;
    }
}

// Обмен данными между соседними процессами
void ExchangeData(double *pProcData, int BlockSize, int DualRank,
                  double *pDualData, int DualBlockSize) {
    MPI_Sendrecv(pProcData, BlockSize, MPI_DOUBLE, DualRank, 0,
                 pDualData,  DualBlockSize, MPI_DOUBLE, DualRank, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

// Проверка распределения данных
void TestDistribution(double *pData, int DataSize,
                      double *pProcData, int BlockSize) {
    MPI_Barrier(MPI_COMM_WORLD);

    if (ProcRank == 0) {
        printf("Initial data:\n");
        PrintData(pData, DataSize);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    for (int i = 0; i < ProcNum; ++i) {
        if (ProcRank == i) {
            printf("ProcRank = %d\n", ProcRank);
            printf("Block:\n");
            PrintData(pProcData, BlockSize);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

// Параллельный вывод данных
void ParallelPrintData(double *pProcData, int BlockSize) {
    for (int i = 0; i < ProcNum; ++i) {
        if (ProcRank == i) {
            printf("ProcRank = %d\n", ProcRank);
            printf("Proc sorted data:\n");
            PrintData(pProcData, BlockSize);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

// Проверка результата параллельной сортировки
void TestResult(double *pData, double *pSerialData, int DataSize) {
    MPI_Barrier(MPI_COMM_WORLD);

    if (ProcRank == 0) {
        SerialBubbleSort(pSerialData, DataSize);
        // SerialStdSort(pSerialData, DataSize);

        if (!CompareData(pData, pSerialData, DataSize))
            printf("The results of serial and parallel algorithms are "
                   "NOT identical. Check your code\n");
        else
            printf("The results of serial and parallel algorithms are "
                   "identical\n");
    }
}

// ---------- вспомогательные функции ----------

void CopyData(const double *src, int n, double *dst) {
    for (int i = 0; i < n; ++i)
        dst[i] = src[i];
}

void PrintData(const double *data, int n) {
    for (int i = 0; i < n; ++i)
        printf("%8.3f ", data[i]);
    printf("\n");
}

void SerialBubbleSort(double *data, int n) {
    for (int i = 1; i < n; ++i)
        for (int j = 0; j < n - i; ++j)
            if (data[j] > data[j + 1]) {
                double tmp   = data[j];
                data[j]      = data[j + 1];
                data[j + 1]  = tmp;
            }
}

void SerialStdSort(double *data, int n) {
    std::sort(data, data + n);
}

bool CompareData(const double *a, const double *b, int n) {
    const double eps = 1e-9;
    for (int i = 0; i < n; ++i)
        if (std::fabs(a[i] - b[i]) > eps)
            return false;
    return true;
}
