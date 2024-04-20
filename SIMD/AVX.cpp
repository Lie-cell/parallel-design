#include <iostream>
#include <cstdlib>
#include <ctime>
#include <Windows.h>
#include <immintrin.h>
using namespace std;

constexpr int N = 10;

float m[N][N];

void m_reset() {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            m[i][j] = 0;
        }
        m[i][i] = 1.0;
        for (int j = i + 1; j < N; j++) {
            m[i][j] = static_cast<float>(rand()) / RAND_MAX;
        }
    }
    for (int k = 0; k < N; k++) {
        for (int i = k + 1; i < N; i++) {
            for (int j = 0; j < N; j++) {
                m[i][j] += m[k][j];
            }
        }
    }
}
//AVX指令集
/*void gaussian_elimination() {
    for (int k = 0; k < N; ++k) {
        __m256 vt = _mm256_set1_ps(m[k][k]); // duplicate m[k,k] into a AVX register

        // Normalize the current row
        for (int j = k + 1; j < N; j += 8) { // Use 8-way vectorization
            __m256 va = _mm256_loadu_ps(&m[k][j]); // load 8 floats from m[k,j]
            va = _mm256_div_ps(va, vt);    // va = va / vt
            _mm256_storeu_ps(&m[k][j], va); // store va back to m[k,j]
        }

        m[k][k] = 1.0f; // Set m[k,k] to 1.0

        // Reduce the other rows
        for (int i = k + 1; i < N; ++i) {
            __m256 vaik = _mm256_set1_ps(m[i][k]); // duplicate m[i,k] into a AVX register
            for (int j = k + 1; j < N; j += 8) { // Use 8-way vectorization
                __m256 vakj = _mm256_loadu_ps(&m[k][j]); // load 8 floats from m[k,j]
                __m256 vaij = _mm256_loadu_ps(&m[i][j]); // load 8 floats from m[i,j]
                __m256 vx = _mm256_mul_ps(vakj, vaik);  // vx = vakj * vaik
                vaij = _mm256_sub_ps(vaij, vx);  // vaij = vaij - vx
                _mm256_storeu_ps(&m[i][j], vaij); // store vaij back to m[i,j]
            }
            m[i][k] = 0.0f; // Set m[i,k] to 0.0
        }
    }
}*/
//AVX-512
void gaussian_elimination() {
    for (int k = 0; k < N; ++k) {
        __m512 vt = _mm512_set1_ps(m[k][k]); // duplicate m[k,k] into a AVX-512 register

        // Normalize the current row
        for (int j = k + 1; j < N; j += 16) { // Use 16-way vectorization
            __m512 va = _mm512_loadu_ps(&m[k][j]); // load 16 floats from m[k,j]
            va = _mm512_div_ps(va, vt);    // va = va / vt
            _mm512_storeu_ps(&m[k][j], va); // store va back to m[k,j]
        }

        m[k][k] = 1.0f; // Set m[k,k] to 1.0

        // Reduce the other rows
        for (int i = k + 1; i < N; ++i) {
            __m512 vaik = _mm512_set1_ps(m[i][k]); // duplicate m[i,k] into a AVX-512 register
            for (int j = k + 1; j < N; j += 16) { // Use 16-way vectorization
                __m512 vakj = _mm512_loadu_ps(&m[k][j]); // load 16 floats from m[k,j]
                __m512 vaij = _mm512_loadu_ps(&m[i][j]); // load 16 floats from m[i,j]
                __m512 vx = _mm512_mul_ps(vakj, vaik);  // vx = vakj * vaik
                vaij = _mm512_sub_ps(vaij, vx);  // vaij = vaij - vx
                _mm512_storeu_ps(&m[i][j], vaij); // store vaij back to m[i,j]
            }
            m[i][k] = 0.0f; // Set m[i,k] to 0.0
        }
    }
}


int main() {
    long long head, tail, freq;
    double total_time=0;
    for (int i = 0; i < 100; i++) {
        m_reset();
        QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
        QueryPerformanceCounter((LARGE_INTEGER*)&head);
        gaussian_elimination();
        QueryPerformanceCounter((LARGE_INTEGER*)&tail);
        total_time+=(tail - head) * 1000.0 / freq;
        /*for(int s=0;s<10;s++)
        {
            for(int j=0;j<10;j++){
                cout<<m[s][j]<<" ";
            }cout<<endl;
        }*/
    }
    total_time/=100;
    std::cout << "Average time: " <<total_time << " ms" << std::endl;
    return 0;
}
