// Serial square matrix multiplication C = A * B (no matrix prints)
// - Checks available RAM and refuses too-large sizes
// - Times only the multiply

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <cstring>
#include <cmath>

// --- read MemAvailable from /proc/meminfo (Linux/WSL) ---
static unsigned long long mem_available_bytes() {
    FILE* f = std::fopen("/proc/meminfo", "r");
    if (!f) return 0ULL;
    char line[256];
    unsigned long long kb = 0;
    while (std::fgets(line, sizeof(line), f)) {
        // line format: "MemAvailable:  1234567 kB"
        if (std::strncmp(line, "MemAvailable:", 13) == 0) {
            std::sscanf(line + 13, "%llu", &kb);
            break;
        }
    }
    std::fclose(f);
    return kb * 1024ULL;
}

// --- init: A=1, B=1 ---
void DummyDataInitialization(double* A, double* B, int n) {
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            A[i * n + j] = 1.0;
            B[i * n + j] = 1.0;
        }
}

// --- C = A * B ---
void MatMul(double* A, double* B, double* C, int n) {
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            double acc = 0.0;
            // классический тройной цикл
            for (int k = 0; k < n; ++k)
                acc += A[i * n + k] * B[k * n + j];
            C[i * n + j] = acc;
        }
}

int main() {
    std::printf("Serial matrix multiplication\n");

    int n = 0;
    // --- ввод размера с проверкой ---
    do {
        std::printf("Enter the size of matrices: ");
        if (std::scanf("%d", &n) != 1) {
            std::fprintf(stderr, "Input error (expected integer)\n");
            return 1;
        }
        if (n <= 0) std::printf("Size must be > 0\n");
        if (n <= 0) continue;

        // --- проверка памяти: требуемые байты = 3 * n*n * 8 ---
        const unsigned long long need =
            3ULL * (unsigned long long)n * (unsigned long long)n * sizeof(double);

        const unsigned long long avail = mem_available_bytes(); // может быть 0, если не удалось прочитать
        // возьмём порог 60% от доступной памяти (чтобы оставить место системе)
        const double safety = 0.60;
        if (avail > 0) {
            const unsigned long long limit =
                (unsigned long long)(avail * safety);
            if (need > limit) {
                // посчитаем максимально допустимый n
                // need = 3 * n^2 * 8  => n_max = floor( sqrt(limit / 24) )
                double nmax_d = std::sqrt((double)limit / 24.0);
                int nmax = (int)nmax_d;
                if (nmax < 1) nmax = 1;
                std::printf(
                    "Not enough RAM for n=%d. Requested ~%.2f GB, allowed ~%.2f GB.\n"
                    "Try n <= %d.\n",
                    n,
                    need / (1024.0*1024.0*1024.0),
                    limit / (1024.0*1024.0*1024.0),
                    nmax
                );
                n = 0; // заставим повторить ввод
            }
        } else {
            // если не удалось определить память — просто предупредим при очень крупных n
            if (need > (unsigned long long)1.5e9) {
                std::printf(
                    "Warning: requested ~%.2f GB; may fail on your system.\n",
                    need / (1024.0*1024.0*1024.0)
                );
            }
        }
    } while (n <= 0);

    std::printf("Chosen matrices' size = %d\n", n);

    // --- аллокации с проверкой ---
    double* A = (double*)std::malloc((size_t)n * n * sizeof(double));
    double* B = (double*)std::malloc((size_t)n * n * sizeof(double));
    double* C = (double*)std::malloc((size_t)n * n * sizeof(double));
    if (!A || !B || !C) {
        std::fprintf(stderr, "Allocation failed. Reduce n.\n");
        std::free(A); std::free(B); std::free(C);
        return 1;
    }

    DummyDataInitialization(A, B, n);

    // --- измерение времени ---
    using clk = std::chrono::steady_clock;
    auto t0 = clk::now();
    MatMul(A, B, C, n);
    auto t1 = clk::now();
    double sec = std::chrono::duration<double>(t1 - t0).count(); // секунды с дробной частью
    std::printf("Time of execution: %.6f s\n", sec);


    std::free(A); std::free(B); std::free(C);
    return 0;
}
