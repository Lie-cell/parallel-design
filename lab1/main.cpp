#include <iostream>
#include <windows.h>
#include <stdlib.h>

using namespace std;
const int N= 10240;
double b[N][N], col_sum[N],a[N];
void init(int n)
{
for (int i = 0; i <N; i++){
for (int j = 0; j <N; j++)
b[i][j] = i + j;}
for (int i = 0; i <N; i++)
a[i]=i;
}
int main()
{/*double sumtime;
for(int i=0;i<100;i++){
long long head, tail, freq;*/
init(N);
//QueryPerformanceFrequency((LARGE_INTEGER *)&freq );
      //QueryPerformanceCounter((LARGE_INTEGER *)&head);
for (int i = 0; i < 10000; i++){
col_sum[i] = 0.0;
for (int j =0; j<10000; j++)
col_sum[i] +=b[j][i]*a[j];
}}//QueryPerformanceCounter ((LARGE_INTEGER *)& tail) ;
//sumtime+=(tail-head)*1000.0 / freq;}
 //cout << "\nordCol:" <<sumtime/100<< "ms" << endl;}

