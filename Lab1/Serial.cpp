// Serial matrix–vector multiplication with high-precision timing
// C++17

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <chrono>

// ---------- init helpers ----------
void DummyDataInit(double* pMatrix, double* pVector, int n) {
    for (int i = 0; i < n; ++i) {
        pVector[i] = 1.0;
        for (int j = 0; j < n; ++j)
            pMatrix[i * n + j] = static_cast<double>(i);
    }
}

void RandomDataInit(double* pMatrix, double* pVector, int n) {
    std::srand(static_cast<unsigned>(std::time(nullptr)));
    for (int i = 0; i < n; ++i) {
        pVector[i] = std::rand() / 1000000.0;
        for (int j = 0; j < n; ++j)
            pMatrix[i * n + j] = std::rand() / 1000000.0;
    }
}

// ---------- io helpers ----------
void PrintMatrix(double* pMatrix, int rows, int cols) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j)
            std::printf("%9.4f ", pMatrix[i * cols + j]);
        std::printf("\n");
    }
}

void PrintVector(double* pVector, int n) {
    for (int i = 0; i < n; ++i) std::printf("%7.4f ", pVector[i]);
    std::printf("\n");
}

// ---------- compute ----------
void ResultCalculation(double* pMatrix, double* pVector, double* pResult, int n) {
    for (int i = 0; i < n; ++i) {
        double acc = 0.0;
        for (int j = 0; j < n; ++j)
            acc += pMatrix[i * n + j] * pVector[j];
        pResult[i] = acc;
    }
}

// ---------- main ----------
int main() {
    using clock = std::chrono::steady_clock;

    double* pMatrix = nullptr;
    double* pVector = nullptr;
    double* pResult = nullptr;
    int n = 0;

    std::cout << "Serial matrix-vector multiplication program\n";

    // ----- input size -----
    do {
        std::cout << "Enter the size of the initial objects: ";
        if (!(std::cin >> n)) { std::cerr << "Input error.\n"; return 1; }
        if (n <= 0) std::cout << "Size must be > 0\n";
    } while (n <= 0);
    std::cout << "Chosen object size: " << n << "\n";

    // ----- allocate -----
    pMatrix = new double[n * n];
    pVector = new double[n];
    pResult = new double[n];

    // ----- init data -----
    // поменяй на DummyDataInit(...) если хочется детерминированных значений
    RandomDataInit(pMatrix, pVector, n);

    // ----- optional prints (only small sizes) -----
    if (n <= 10) {
        std::cout << "\nInitial matrix:\n";
        PrintMatrix(pMatrix, n, n);
        std::cout << "Initial vector:\n";
        PrintVector(pVector, n);
    }

    // ----- timing -----
    // подберём число повторов, чтобы измерение было не ~0 мс
    int reps = (n <= 100 ? 2000 : (n <= 300 ? 400 : (n <= 1000 ? 100 : 10)));

    auto t0 = clock::now();
    for (int r = 0; r < reps; ++r)
        ResultCalculation(pMatrix, pVector, pResult, n);
    auto t1 = clock::now();

    // общее время и среднее за один прогон
    const double total_ms =
        std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count() / 1000.0;
    const double avg_ms = (total_ms / reps)/1000;

    // ----- output -----
    if (n <= 10) {
        std::cout << "\nResult vector:\n";
        PrintVector(pResult, n);
    }
    
    std::printf("Average time per run:   %.6f s\n", avg_ms);

    // ----- cleanup -----
    delete[] pMatrix;
    delete[] pVector;
    delete[] pResult;
    return 0;
}
