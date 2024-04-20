#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<map>
#include<windows.h>
#include<tmmintrin.h>
#include<xmmintrin.h>
#include<emmintrin.h>
#include<pmmintrin.h>
#include<smmintrin.h>
#include<nmmintrin.h>
#include<immintrin.h>
using namespace std;

const int maxsize = 3000;
const int maxrow = 60000;
const int numBasis = 40000;
long long head, tail, freq;

map<int, int*>iToBasis;
map<int, int>iToFirst;
map<int, int*>ans;

fstream RowFile("1.txt", ios::in | ios::out);
fstream BasisFile("2.txt", ios::in | ios::out);


int gRows[maxrow][maxsize];
int gBasis[numBasis][maxsize];

void reset() {
	memset(gRows, 0, sizeof(gRows));
	memset(gBasis, 0, sizeof(gBasis));
	RowFile.close();
	BasisFile.close();
	RowFile.open("1.txt", ios::in | ios::out);
	BasisFile.open("2.txt", ios::in | ios::out);
	iToBasis.clear();
	iToFirst.clear();
	ans.clear();

}

void readBasis() {
    for (int i = 0; i < numBasis; ++i) {
        string line;
        if (getline(BasisFile, line)) {
            stringstream ss(line);
            int currentRow = 0;
            bool isNewRow = false;
            int pos;
            while (ss >> pos) {
                if (!isNewRow) {
                    currentRow = pos;
                    isNewRow = true;
                    iToBasis.insert(pair<int, int*>(currentRow, gBasis[currentRow]));
                }
                int index = pos / 32;
                int offset = pos % 32;
                gBasis[currentRow][index] |= (1 << offset);
            }
        } else {
            cout << "Failed to read basis line " << i << endl;
            return;
        }
    }
}


int readRowsFrom(int pos) {
    iToFirst.clear();
    if (RowFile.is_open())
        RowFile.close();
    RowFile.open("1.txt", ios::in | ios::out);
    memset(gRows, 0, sizeof(gRows));
    string line;
    for (int i = 0; i < pos; i++) {
        getline(RowFile, line);
    }
    for (int i = pos; i < pos + maxrow; i++) {
        int tmp;
        getline(RowFile, line);
        if (line.empty()) {
            cout << "Read elimination row " << i << " line" << endl;
            return i;
        }
        bool isFirstElementFound = false;
        stringstream s(line);
        int rowNumber = i - pos;
        while (s >> tmp) {
            if (!isFirstElementFound) {
                iToFirst[rowNumber] = tmp;
                isFirstElementFound = true;
            }
            int index = tmp / 32;
            int offset = tmp % 32;
            gRows[rowNumber][index] |= (1 << offset);
        }
    }
    cout << "Read max rows" << endl;
    return -1;
}


void update(int row) {
    bool isFirstElementFound = false;
    for (int i = maxsize - 1; i >= 0; i--) {
        if (gRows[row][i] == 0) {
            continue;
        } else {
            isFirstElementFound = true;
            int position = i * 32;
            int offset = 0;
            for (int k = 31; k >= 0; k--) {
                if (gRows[row][i] & (1 << k)) {
                    offset = k;
                    break;
                }
            }
            int newFirstElement = position + offset;
            iToFirst.erase(row);
            iToFirst.insert(pair<int, int>(row, newFirstElement));
            break;
        }
    }
    if (!isFirstElementFound) {
        iToFirst.erase(row);
    }
}

void writeResult(ofstream& out) {
    for (auto it = ans.rbegin(); it != ans.rend(); it++) {
        int* result = it->second;
        int max = it->first / 32 + 1;
        for (int i = max; i >= 0; i--) {
            if (result[i] == 0) {
                continue;
            }
            int position = i * 32;
            for (int k = 31; k >= 0; k--) {
                if (result[i] & (1 << k)) {
                    out << k + position << " ";
                }
            }
        }
        out << endl;
    }
}

void grobner() {
	long long t1, t2;
	int b = 0;
	int f;
	f = readRowsFrom(b);
	int n = (f == -1) ? maxrow : f;
	for (int i = 0; i < n; i++) {
		while (iToFirst.find(i) != iToFirst.end()) {
			int fi = iToFirst.find(i)->second;
			if (iToBasis.find(fi) != iToBasis.end()) {
				int* bs = iToBasis.find(fi)->second;
				for (int j = 0; j < maxsize; j++) {
					gRows[i][j] = gRows[i][j] ^ bs[j];
				}
			}
			else {
				for (int j = 0; j < maxsize; j++) {
					gBasis[fi][j] = gRows[i][j];
				}
				iToBasis.insert(pair<int, int*>(fi, gBasis[fi]));
				ans.insert(pair<int, int*>(fi, gBasis[fi]));
				iToFirst.erase(i);
				continue;
			}
			update(i);
		}
	}
}

void AVXgrobner() {
	long long rb, re;
	int b = 0;
	int f = readRowsFrom(b);
	int n = (f == -1) ? maxrow : f;
	for (int i = 0; i < n; i++) {
		while (iToFirst.find(i) != iToFirst.end()) {
			int fst = iToFirst.find(i)->second;
			if (iToBasis.find(fst) != iToBasis.end()) {
				int* bs = iToBasis.find(fst)->second;
				int j = 0;
				for (; j + 8 < maxsize; j += 8) {
					__m256i vij = _mm256_loadu_si256((__m256i*) & gRows[i][j]);
					__m256i vj = _mm256_loadu_si256((__m256i*) & bs[j]);
					__m256i vx = _mm256_xor_si256(vij, vj);
					_mm256_storeu_si256((__m256i*) & gRows[i][j], vx);
				}
				for (; j < maxsize; j++) {
					gRows[i][j] = gRows[i][j] ^ bs[j];
				}
				update(i);
			}
			else {
				int j = 0;
				for (; j + 8 < maxsize; j += 8) {
					__m256i vij = _mm256_loadu_si256((__m256i*) & gRows[i][j]);
					_mm256_storeu_si256((__m256i*) & gBasis[fst][j], vij);
				}
				for (; j < maxsize; j++) {
					gBasis[fst][j] = gRows[i][j];
				}
				iToBasis.insert(pair<int, int*>(fst, gBasis[fst]));
				ans.insert(pair<int, int*>(fst, gBasis[fst]));
				iToFirst.erase(i);
			}
		}
	}
}

int main() {
	double time1 = 0;
	double time2 = 0;
	for (int i = 0; i < 1; i++) {
		ofstream out("123.txt");
		ofstream out1("123(AVX).txt");
		QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
		readBasis();
		QueryPerformanceCounter((LARGE_INTEGER*)&head);
		grobner();
		QueryPerformanceCounter((LARGE_INTEGER*)&tail);
		time1 += (tail - head) * 1000 / freq;
		writeResult(out);
		reset();
		readBasis();
		QueryPerformanceCounter((LARGE_INTEGER*)&head);
		AVXgrobner();
		QueryPerformanceCounter((LARGE_INTEGER*)&tail);
		time2 += (tail - head) * 1000 / freq;
		writeResult(out1);
		reset();
		out.close();
		out1.close();
	}
	cout << "time1:" << time1 / 1 << endl << "time2:" << time2 / 1;
}

