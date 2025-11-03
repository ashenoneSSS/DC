#include <cstdio>
#include <cstdlib>
#include <ctime>

// ---- налаштування тестів ----
const double InfinitiesPercent = 50.0;

// ---- допоміжне ----
int Min(int A, int B) {
    int Result = (A < B) ? A : B;
    if ((A < 0) && (B >= 0)) Result = B;
    if ((B < 0) && (A >= 0)) Result = A;
    if ((A < 0) && (B < 0)) Result = -1;
    return Result;
}

// виділення пам'яті + ручне введення Size
void ProcessInitialization(int *&pMatrix, int &Size) {
    do {
        std::printf("Enter the number of vertices: ");
        if (std::scanf("%d", &Size) != 1) {
            // очистити сміття у stdin
            int ch; while ((ch = std::getchar()) != '\n' && ch != EOF) {}
            Size = 0;
        }
        if (Size <= 0)
            std::printf("The number of vertices should be greater then zero\n");
    } while (Size <= 0);

    std::printf("Using graph with %d vertices\n", Size);
    pMatrix = new int[Size * Size];
    // ініціалізація
    // RandomDataInitialization(pMatrix, Size);
    for (int i = 0; i < Size; i++)
        for (int j = i; j < Size; j++) {
            if (i == j) pMatrix[i * Size + j] = 0;
            else if (i == 0) pMatrix[i * Size + j] = j;
            else pMatrix[i * Size + j] = -1;
            pMatrix[j * Size + i] = pMatrix[i * Size + j];
        }
}

// виділення пам'яті + Size уже заданий (для автотестів)
void ProcessInitializationTest(int *&pMatrix, int &Size) {
    if (Size <= 0)
        std::printf("The number of vertices should be greater then zero\n");
    std::printf("Using graph with %d vertices\n", Size);
    pMatrix = new int[Size * Size];
    // DummyDataInitialization
    for (int i = 0; i < Size; i++)
        for (int j = i; j < Size; j++) {
            if (i == j) pMatrix[i * Size + j] = 0;
            else if (i == 0) pMatrix[i * Size + j] = j;
            else pMatrix[i * Size + j] = -1;
            pMatrix[j * Size + i] = pMatrix[i * Size + j];
        }
}

// випадкова ініціалізація (якщо захочеш увімкнути)
void RandomDataInitialization(int *pMatrix, int Size) {
    std::srand((unsigned)std::time(nullptr));
    for (int i = 0; i < Size; i++)
        for (int j = 0; j < Size; j++)
            if (i != j) {
                if ((std::rand() % 100) < InfinitiesPercent)
                    pMatrix[i * Size + j] = -1;
                else
                    pMatrix[i * Size + j] = std::rand() + 1;
            } else {
                pMatrix[i * Size + j] = 0;
            }
}

void ProcessTermination(int *pMatrix) {
    delete [] pMatrix;
}

// послідовний Флойд–Уоршелл для ваг з -1 як «∞»
void SerialFloyd(int *pMatrix, int Size) {
    int t1, t2;
    for (int k = 0; k < Size; k++)
        for (int i = 0; i < Size; i++)
            for (int j = 0; j < Size; j++)
                if ((pMatrix[i * Size + k] != -1) &&
                    (pMatrix[k * Size + j] != -1)) {
                    t1 = pMatrix[i * Size + j];
                    t2 = pMatrix[i * Size + k] + pMatrix[k * Size + j];
                    pMatrix[i * Size + j] = Min(t1, t2);
                }
}

int main() {
    std::printf("Serial Floyd algorithm\n");

    int *pMatrix = nullptr;
    int Size;
    std::clock_t start, finish;
    double duration = 0.0;

    
    for (int i = 0; i < 7; i++) {
        
        if (i == 0) Size = 10;
        else        Size = 500 + (i - 1) * 100;

        ProcessInitializationTest(pMatrix, Size);

        start = std::clock();
        SerialFloyd(pMatrix, Size);
        finish = std::clock();

        duration = (finish - start) / double(CLOCKS_PER_SEC);
        std::printf("Time of execution: %f\n", duration);

        ProcessTermination(pMatrix);
    }
    return 0;
}
