#include <iostream>
#include <cstdlib>
#include <ctime>
#include <Windows.h>
#include <immintrin.h> // For SSE intrinsics

const int N = 1000; // 矩阵大小

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
//优化算法
/*void gaussian_elimination() {
    for (int k = 0; k < N; ++k) {
        __m128 vt = _mm_set_ps1(m[k][k]); // duplicate m[k,k] into a SSE register

        // Normalize the current row
        for (int j = k + 1; j < N; j += 4) { // Use 4-way vectorization
            __m128 va = _mm_loadu_ps(&m[k][j]); // load 4 floats from m[k,j]
            va = _mm_div_ps(va, vt);    // va = va / vt
            _mm_storeu_ps(&m[k][j], va); // store va back to m[k,j]
        }

        m[k][k] = 1.0f; // Set m[k,k] to 1.0

        // Reduce the other rows
        for (int i = k + 1; i < N; ++i) {
            __m128 vaik = _mm_set_ps1(m[i][k]); // duplicate m[i,k] into a SSE register
            for (int j = k + 1; j < N; j += 4) { // Use 4-way vectorization
                __m128 vakj = _mm_loadu_ps(&m[k][j]); // load 4 floats from m[k,j]
                __m128 vaij = _mm_loadu_ps(&m[i][j]); // load 4 floats from m[i,j]
                __m128 vx = _mm_mul_ps(vakj, vaik);  // vx = vakj * vaik
                vaij = _mm_sub_ps(vaij, vx);  // vaij = vaij - vx
                _mm_storeu_ps(&m[i][j], vaij); // store vaij back to m[i,j]
            }
            m[i][k] = 0.0f; // Set m[i,k] to 0.0
        }
    }
}*/
//优化对齐
void gaussian_elimination() {
    for (int k = 0; k < N - 1; k++) {
        __m128 vt = _mm_set1_ps(m[k][k]);
        int j_aligned = (k + 1 + 3) & ~3; // Align to a multiple of 4

        // Process aligned part, iterating over 4 elements at a time
        for (int j = j_aligned; j < N; j += 4) {
            __m128 va = _mm_load_ps(&m[k][j]);
            va = _mm_div_ps(va, vt);
            _mm_store_ps(&m[k][j], va);
        }

        // Process unaligned part
        for (int j = k + 1; j < j_aligned; j++) {
            m[k][j] /= m[k][k];
        }

        m[k][k] = 1.0f;

        // Process remaining rows
        for (int i = k + 1; i < N; i++) {
            __m128 vaik = _mm_set1_ps(m[i][k]);

            // Process aligned part, iterating over 4 elements at a time
            for (int j = j_aligned; j < N; j += 4) {
                __m128 vakj = _mm_set_ps(m[k][j + 3], m[k][j + 2], m[k][j + 1], m[k][j]);
                __m128 vaij = _mm_load_ps(&m[i][j]);
                vaij = _mm_sub_ps(vaij, _mm_mul_ps(vaik, vakj));
                _mm_store_ps(&m[i][j], vaij);
            }

            // Process unaligned part
            for (int j = k + 1; j < j_aligned; j++) {
                m[i][j] -= m[i][k] * m[k][j];
            }

            m[i][k] = 0.0f;
        }
    }
}
/*void gaussian_elimination() {
    for (int k = 0; k < N; ++k) {
        __m128 pivot = _mm_set_ps1(m[k][k]);

        // Normalize the current row
        for (int j = k + 1; j < N; ++j) {
            m[k][j] /= m[k][k];
        }
        m[k][k] = 1.0f; // Set m[k,k] to 1.0

        // Reduce the other rows
        for (int i = k + 1; i < N; ++i) {
            float factor = m[i][k];
            for (int j = k + 1; j < N; ++j) {
                m[i][j] -= factor * m[k][j];
            }
            m[i][k] = 0.0f; // Set m[i,k] to 0.0
        }
    }
}*/

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
        /*for(int s=0;s<10;s++){
        	for (int j=0;j<10;j++)
        	std::cout<<m[s][j]<<" ";
        	std::cout<<std::endl;
		}*/
    }
    total_time/=100;
    std::cout << "Average time: " <<total_time << " ms" << std::endl;
    return 0;
}



