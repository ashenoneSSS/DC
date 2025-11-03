#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <algorithm>

using namespace std;

const double RandomDataMultiplier = 1000.0;

// Прототипи можна залишити для читабельності, але вони в тому ж файлі
void ProcessInitialization(double *&pData, int& DataSize);
void ProcessTermination(double *pData);
void DummyDataInitialization(double*& pData, int& DataSize);
void RandomDataInitialization(double *&pData, int& DataSize);
void SerialBubble(double *pData, int DataSize);
void SerialStdSort(double *pData, int DataSize);

int main(int argc, char *argv[]) {
    double *pData = nullptr;
    int DataSize = 0;
    clock_t start, finish;
    double duration = 0.0;

    printf("Serial bubble sort program\n");

    // Ініціалізація
    ProcessInitialization(pData, DataSize);
    // printf("Data before sorting\n");
    // PrintData(pData, DataSize);

    start = clock();
    // Варіант 1: свій bubble sort
    // SerialBubble(pData, DataSize);

    // Варіант 2: стандартне сортування
    SerialStdSort(pData, DataSize);
    finish = clock();

    // printf("Data after sorting\n");
    // PrintData(pData, DataSize);

    duration = (finish - start) / double(CLOCKS_PER_SEC);
    printf("Time of execution: %f\n", duration);

    // Завершення
    ProcessTermination(pData);

    return 0;
}

// Функція виділення пам’яті і задання початкових значень
void ProcessInitialization(double *&pData, int& DataSize) {
    while (true) {
        printf("Enter the size of data to be sorted: ");

        int rc = scanf("%d", &DataSize);
        if (rc != 1) {
            // Нормально очищаем некорректный ввод
            printf("Input error: please enter an integer\n");
            int ch;
            while ((ch = getchar()) != '\n' && ch != EOF) { }
            continue;
        }

        if (DataSize <= 0) {
            printf("Data size should be greater than zero\n");
            continue;
        }

        break;
    }

    printf("Sorting %d data items\n", DataSize);

    pData = new double[DataSize];

    // Простое задання даних
    // DummyDataInitialization(pData, DataSize);
    // Випадкові дані
    RandomDataInitialization(pData, DataSize);
}


// Завершення обчислювального процесу
void ProcessTermination(double *pData) {
    delete [] pData;
}

// Просте задання вихідних даних (убывающая последовательность)
void DummyDataInitialization(double*& pData, int& DataSize) {
    for (int i = 0; i < DataSize; i++)
        pData[i] = DataSize - i;
}

// Ініціалізація випадковими даними
void RandomDataInitialization(double *&pData, int& DataSize) {
    srand((unsigned)time(0));

    for (int i = 0; i < DataSize; i++)
        pData[i] = double(rand()) / RAND_MAX * RandomDataMultiplier;
}

// Послідовний bubble sort
void SerialBubble(double *pData, int DataSize) {
    double Tmp;

    for (int i = 1; i < DataSize; i++)
        for (int j = 0; j < DataSize - i; j++)
            if (pData[j] > pData[j + 1]) {
                Tmp = pData[j];
                pData[j] = pData[j + 1];
                pData[j + 1] = Tmp;
            }
}

// Сортування стандартним алгоритмом
void SerialStdSort(double *pData, int DataSize) {
    sort(pData, pData + DataSize);
}
