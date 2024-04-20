#include <iostream>
#include <cstdlib>
#include <ctime>
#include <Windows.h>
#include <immintrin.h>
using namespace std;

constexpr int N = 100;

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
void gaussian_elimination() {
    for (int k = 0; k < N; ++k) {
        float vt = m[k][k]; e
        for (int j = k + 1; j < N; ++j) {
            m[k][j] /= vt;
        }

        m[k][k] = 1.0f;

        for (int i = k + 1; i < N; ++i) {
            float vaik = m[i][k];
            for (int j = k + 1; j < N; ++j) {
                m[i][j] -= m[k][j] * vaik;
            }
            m[i][k] = 0.0f;
        }
    }
}
int main() {
    long long head, tail, freq;
    double total_time=0;
    for (int i = 0; i < 100; i++) {

        QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
        QueryPerformanceCounter((LARGE_INTEGER*)&head);
        m_reset();
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
