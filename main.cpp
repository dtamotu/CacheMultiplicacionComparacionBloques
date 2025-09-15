#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>
#include <algorithm>
using namespace std;
using namespace std::chrono;

// ==============================
// Producto matriz-vector
// ==============================
pair<long long, long long> evaluarMV(int MAX) {
    vector<vector<double>> A(MAX, vector<double>(MAX, 1.0));
    vector<double> x(MAX, 1.0);
    vector<double> y(MAX, 0.0);


    fill(y.begin(), y.end(), 0.0);
    auto start1 = high_resolution_clock::now();
    for (int i = 0; i < MAX; i++)
        for (int j = 0; j < MAX; j++)
            y[i] += A[i][j] * x[j];
    auto end1 = high_resolution_clock::now();
    auto duration1 = duration_cast<milliseconds>(end1 - start1).count();

    fill(y.begin(), y.end(), 0.0);
    auto start2 = high_resolution_clock::now();
    for (int j = 0; j < MAX; j++)
        for (int i = 0; i < MAX; i++)
            y[i] += A[i][j] * x[j];
    auto end2 = high_resolution_clock::now();
    auto duration2 = duration_cast<milliseconds>(end2 - start2).count();

    return {duration1, duration2};
}

// ==============================
// Multiplicación de matrices clásica (3 bucles)
// ==============================
long long multiplicacionClasica(int N) {
    vector<vector<double>> A(N, vector<double>(N, 1.0));
    vector<vector<double>> B(N, vector<double>(N, 1.0));
    vector<vector<double>> C(N, vector<double>(N, 0.0));

    auto start = high_resolution_clock::now();
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++)
                C[i][j] += A[i][k] * B[k][j];
    auto end = high_resolution_clock::now();
    return duration_cast<milliseconds>(end - start).count();
}

// ==============================
// Multiplicación de matrices con 6 bucles (bloques / tiling)
// ==============================
long long multiplicacionBloques(int N, int B) {
    vector<vector<double>> A(N, vector<double>(N, 1.0));
    vector<vector<double>> Bm(N, vector<double>(N, 1.0));
    vector<vector<double>> C(N, vector<double>(N, 0.0));

    auto start = high_resolution_clock::now();
    for (int ii = 0; ii < N; ii += B)
        for (int jj = 0; jj < N; jj += B)
            for (int kk = 0; kk < N; kk += B)
                for (int i = ii; i < min(ii + B, N); i++)
                    for (int j = jj; j < min(jj + B, N); j++)
                        for (int k = kk; k < min(kk + B, N); k++)
                            C[i][j] += A[i][k] * Bm[k][j];
    auto end = high_resolution_clock::now();
    return duration_cast<milliseconds>(end - start).count();
}

// ==============================
// MAIN
// ==============================
int main() {
    vector<int> tamanos = {100, 500, 1000}; // puedes agregar más tamaños
    int bloque = 64; // tamaño de bloque para 6 bucles

    ofstream archivo("resultados_comparacion.csv");
    archivo << "Tamano,"
            << "RowMajor_ms,ColMajor_ms,"
            << "Mat3Bucles_ms,Mat6Bucles_ms\n";

   for (int N : tamanos) {
        cout << "\nMatriz " << N << "x" << N << ":\n";


        auto [rowMajor, colMajor] = evaluarMV(N);
        cout << "  Producto MV -> Row-major: " << rowMajor
             << " ms | Col-major: " << colMajor << " ms\n";


        auto tiempo3 = multiplicacionClasica(N);
        cout << "  Multiplicación 3 bucles: " << tiempo3 << " ms\n";


        auto tiempo6 = multiplicacionBloques(N, bloque);
        cout << "  Multiplicación 6 bucles: " << tiempo6 << " ms\n";

        archivo << N << ","
                << rowMajor << "," << colMajor << ","
                << tiempo3 << "," << tiempo6 << "\n";
    }

    archivo.close();
    cout << "\nResultados guardados en 'resultados_comparacion.csv'\n";
    return 0;
}