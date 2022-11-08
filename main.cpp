#include <iostream>
#include <cinttypes>
#include <cmath>
#include <vector>
#include <cstring>
#include <algorithm>

#define baisset(ba, i) ((ba)[(i)>>6] & (1ull << ((i) & 63))) != 0
#define baisclear(ba, i) ((ba)[(i)>>6] & (1ull << ((i) & 63))) == 0
#define baset(ba, i) (ba)[(i)>>6] |= (1ull << ((i) & 63))

using namespace std;

int maxp;
#define MAXSIEVE 8000000

#define NUMPRIMES 2262 //Computed, can't know beforehand

int primesList[NUMPRIMES];
int numprimes;
int primes1indexes[NUMPRIMES];
int nprimes1;
int primes2indexes[NUMPRIMES];
int nprimes2;


void fillprimes(int upperBound);
int is_prime(uint32_t x);

int mult10[10];
inline int whichmult(int v, int i) {
	while (mult10[i] <= v)
		i++;
	return i;
}

vector<int> concatenated[NUMPRIMES];

void init() {
	int *primesIndexLists[2] = { primes1indexes, primes2indexes };
	int *nprimes[2] = { &nprimes1, &nprimes2 };

	mult10[0] = 1;
	for (int i = 1; i < 10; i++)
		mult10[i] = mult10[i - 1] * 10;

	fillprimes(MAXSIEVE*1.01);

	primes1indexes[nprimes1++] = 1;
	primes2indexes[nprimes2++] = 1;
	for (int i = 3; i < numprimes; i++) {
		int p = primesList[i];
		if (p % 3 == 1) {
			primes1indexes[nprimes1++] = i;
		}
		else {
			primes2indexes[nprimes2++] = i;
		}
	}

	int p1, imultp1, multp1, imultp2, p2;
	for (int z = 0; z < 2; z++) {
		int *primeIndexList = primesIndexLists[z];
		int np = nprimes[z][0];
		imultp1 = 0;
		for (int i = 0; i < np; i++) {
			p1 = primesList[primeIndexList[i]];
			imultp1 = whichmult(p1, imultp1);
			imultp2 = imultp1;
			for (int j = i + 1; j < np; j++) {
				p2 = primesList[primeIndexList[j]];
				imultp2 = whichmult(p2, imultp2);

				if (is_prime(mult10[imultp2]*p1 + p2) && is_prime(p2*mult10[imultp1] + p1)) {
					concatenated[primeIndexList[i]].push_back(primeIndexList[j]);
				}
			}
		}
	}
}


int solutions[1000000];
int nsol;

void solve(int k) {
	sort(concatenated[1].begin(), concatenated[1].end());
	for (int i1 = 0; i1 < numprimes; i1++) {
		if (!concatenated[i1].size()) continue;
		for (int i2 : concatenated[i1]) {
			if (!concatenated[i2].size()) continue;
			for (int i3 : concatenated[i2]) {
				if (!binary_search(concatenated[i1].begin(), concatenated[i1].end(), i3)) continue;

				if (k == 3) {
					solutions[nsol++] = primesList[i1] + primesList[i2] + primesList[i3];
					continue;
				}


				if (!concatenated[i3].size()) continue;

				for (int i4 : concatenated[i3]) {
					if (!binary_search(concatenated[i1].begin(), concatenated[i1].end(), i4)) continue;
					if (!binary_search(concatenated[i2].begin(), concatenated[i2].end(), i4)) continue;

					if (k == 4) {
						solutions[nsol++] = primesList[i1] + primesList[i2] + primesList[i3] + primesList[i4];
						continue;
					}

					if (!concatenated[i4].size()) continue;

					for (int i5 : concatenated[i4]) {
						if (!binary_search(concatenated[i1].begin(), concatenated[i1].end(), i5)) continue;
						if (!binary_search(concatenated[i2].begin(), concatenated[i2].end(), i5)) continue;
						if (!binary_search(concatenated[i3].begin(), concatenated[i3].end(), i5)) continue;
						solutions[nsol++] = primesList[i1] + primesList[i2] + primesList[i3] + primesList[i4] + primesList[i5];
					}

				}
			}
		}


	}

}


int main()
{
	int n = 20000; 
    int k = 5;
	maxp = n;
	init();
	solve(k);

	sort(solutions, solutions + nsol);
	for (int i = 0; i < nsol; i++)
		cout << solutions[i] << endl;

	return 0;
}



/*****************************************************************************/

uint64_t iscompound[130000];
void fillprimes(int upperBound) {
	int i, j;

	baset(iscompound, 0);
	baset(iscompound, 1);
	int sqrtUb = sqrt(upperBound);
	for (i = 2; i <= sqrtUb; i++) {
		if (baisclear(iscompound, i)) {
			for (j = i*i; j <= upperBound; j += i) {
				baset(iscompound, j);
			}
		}
	}

	for (i = 2; i <= maxp; i++) {
		if (baisclear(iscompound, i)) {
			primesList[numprimes++] = i;
		}
	}

}

//Miller-Rabin, from https://people.ksp.sk/~misof/primes/FJ32_1024.cc
uint8_t bases_2[] = { 17, 11, 6, 60, 7, 13, 11, 34, 13, 2, 3, 37, 13, 11,
38, 2, 7, 105, 2, 7, 42, 11, 7, 3, 6, 15, 53, 44, 6, 6, 5, 15, 54, 7, 35,
10, 10, 15, 10, 10, 17, 17, 11, 10, 15, 43, 7, 5, 5, 3, 7, 43, 34, 2, 34,
2, 68, 53, 39, 10, 7, 6, 11, 2, 5, 2, 7, 2, 6, 5, 15, 40, 3, 5, 5, 2, 2,
10, 47, 13, 7, 43, 6, 7, 5, 6, 6, 13, 6, 35, 6, 15, 6, 13, 40, 10, 11, 2,
7, 2, 2, 3, 13, 3, 11, 15, 10, 5, 11, 14, 7, 11, 47, 5, 2, 2, 6, 2, 5, 55,
6, 5, 7, 2, 6, 58, 35, 11, 5, 12, 17, 6, 10, 12, 6, 6, 2, 53, 2, 2, 13, 5,
14, 7, 15, 6, 13, 62, 10, 6, 3, 7, 7, 3, 14, 5, 14, 73, 15, 11, 11, 6, 5,
17, 10, 5, 3, 37, 51, 10, 7, 5, 38, 12, 5, 11, 5, 7, 6, 5, 6, 40, 43, 57,
10, 13, 7, 15, 2, 10, 34, 7, 39, 10, 5, 3, 6, 13, 11, 5, 10, 43, 10, 5, 3,
14, 5, 2, 5, 41, 5, 39, 46, 2, 10, 2, 5, 12, 3, 2, 2, 5, 15, 43, 17, 41, 2,
13, 15, 38, 11, 11, 3, 34, 5, 6, 3, 7, 2, 37, 5, 6, 10, 17, 35, 2, 15, 6,
7, 5, 3, 13, 13, 12, 34, 2, 12, 10, 15, 13, 2, 2, 34, 6, 6, 5, 2, 7, 13, 3,
6, 11, 39, 42, 7, 2, 6, 39, 47, 3, 17, 5, 13, 7, 2, 47, 3, 7, 6, 11, 17,
37, 48, 7, 37, 11, 7, 10, 3, 14, 39, 14, 15, 43, 17, 2, 12, 7, 13, 5, 3, 6,
34, 37, 3, 17, 13, 2, 5, 10, 10, 44, 37, 2, 2, 10, 10, 7, 3, 7, 2, 7, 5,
43, 43, 11, 15, 51, 13, 17, 10, 11, 2, 5, 34, 17, 2, 2, 42, 6, 6, 5, 47,
15, 2, 12, 7, 3, 10, 15, 3, 7, 12, 12, 15, 43, 14, 7, 58, 13, 10, 6, 6, 38,
34, 5, 5, 13, 38, 6, 11, 10, 6, 7, 2, 55, 2, 13, 5, 11, 44, 15, 17, 2, 40,
2, 15, 13, 6, 2, 3, 3, 3, 3, 6, 39, 5, 11, 17, 37, 5, 7, 6, 10, 6, 12, 7,
5, 14, 10, 12, 71, 10, 35, 6, 11, 3, 2, 38, 3, 2, 34, 10, 17, 42, 2, 12, 6,
6, 11, 40, 12, 10, 6, 10, 2, 3, 3, 56, 11, 7, 42, 2, 38, 12, 2, 2, 13, 40,
12, 6, 5, 5, 59, 15, 38, 5, 5, 5, 7, 2, 10, 7, 2, 17, 10, 11, 6, 6, 6, 2,
10, 6, 54, 2, 82, 3, 34, 14, 15, 44, 5, 46, 2, 13, 5, 12, 13, 11, 10, 39,
5, 40, 3, 60, 3, 42, 11, 3, 46, 17, 3, 2, 37, 6, 42, 12, 14, 3, 12, 66, 13,
34, 7, 3, 13, 3, 11, 2, 13, 12, 38, 34, 5, 40, 10, 14, 6, 14, 11, 38, 58,
2, 48, 5, 15, 5, 73, 3, 37, 5, 11, 10, 5, 5, 13, 2, 10, 13, 34, 17, 3, 7,
47, 2, 2, 10, 15, 3, 3, 13, 6, 34, 13, 10, 13, 3, 6, 41, 10, 6, 2, 6, 2, 6,
2, 6, 6, 37, 10, 44, 35, 13, 51, 2, 7, 53, 5, 40, 5, 2, 37, 11, 15, 11, 13,
2, 5, 2, 6, 10, 17, 15, 43, 39, 17, 2, 12, 10, 15, 17, 7, 13, 3, 7, 15, 37,
5, 15, 7, 6, 10, 51, 2, 2, 40, 61, 2, 13, 13, 11, 2, 5, 34, 5, 5, 7, 2, 2,
2, 11, 3, 6, 13, 6, 17, 11, 10, 7, 46, 15, 7, 14, 35, 11, 7, 10, 6, 11, 40,
11, 2, 39, 7, 6, 66, 5, 3, 6, 5, 11, 10, 2, 10, 7, 13, 2, 45, 34, 6, 35, 2,
11, 5, 59, 75, 10, 17, 14, 17, 17, 17, 2, 11, 7, 10, 6, 11, 6, 56, 34, 35,
11, 14, 12, 41, 40, 17, 40, 3, 11, 7, 37, 14, 7, 13, 7, 5, 2, 10, 6, 39, 2,
7, 37, 35, 10, 5, 15, 2, 7, 38, 34, 11, 17, 5, 6, 10, 3, 6, 7, 7, 43, 14,
2, 43, 3, 2, 47, 7, 35, 7, 3, 53, 2, 10, 10, 10, 60, 10, 6, 2, 6, 10, 5, 7,
57, 53, 13, 3, 35, 38, 15, 42, 3, 3, 12, 2, 10, 3, 38, 54, 13, 10, 11, 7,
13, 7, 2, 12, 39, 10, 54, 2, 12, 38, 10, 12, 12, 5, 15, 6, 10, 13, 5, 15,
10, 13, 6, 41, 40, 14, 12, 10, 11, 40, 5, 11, 10, 2, 5, 2, 13, 6, 2, 13, 5,
2, 10, 15, 5, 5, 10, 34, 13, 2, 5, 14, 5, 6, 5, 13, 3, 43, 6, 13, 11, 50,
3, 6, 6, 12, 15, 11, 37, 7, 69, 11, 14, 14, 7, 43, 5, 35, 11, 35, 11, 11,
34, 34, 39, 14, 11, 2, 10, 53, 6, 11, 2, 11, 60, 39, 11, 6, 15, 40, 17, 47,
34, 50, 7, 59, 47, 5, 13, 39, 5, 6, 53, 10, 14, 5, 51, 5, 7, 5, 6, 77, 7,
12, 7, 42, 2, 5, 2, 6, 60, 10, 13, 10, 6, 47, 6, 15, 17, 10, 11, 10, 12, 7,
7, 10, 17, 34, 5, 10, 7, 7, 2, 6, 10, 38, 2, 15, 6, 13, 7, 13, 2, 3, 13, 5,
3, 17, 2, 5, 15, 11, 39, 7, 39, 10, 10, 2, 6, 13, 3, 5, 17, 6, 14, 10, 37,
44, 3, 34, 5, 11, 7, 12, 2, 5, 3, 12, 3, 2, 3, 133, 12, 2, 2, 2, 3, 34, 14,
41, 2, 37, 11, 2, 6, 11, 6, 7, 15, 11, 35, 13, 6, 5, 2, 14, 7, 2 };


uint64_t modpow(uint64_t a, uint64_t e, uint64_t n)
{
	uint64_t ans = 1;
	while (e>0) {
		if (e & 1) ans = (ans*a) % n;
		a = (a*a) % n;
		e >>= 1;
	}
	return ans;
}


int is_SPRP(uint64_t a, uint64_t n) {
	int d = 0;
	uint64_t t = n - 1;
	while (t % 2 == 0) {
		++d;
		t >>= 1;
	}
	uint64_t x = modpow(a, t, n);
	if (x == 1) return 1;
	for (int i = 0; i<d; ++i) {
		if (x == n - 1) return 1;
		x = (x*x) % n;
		if (x == 1) return 0;
	}
	return 0;
}


int is_prime(uint32_t x) {
	if (x < MAXSIEVE)
		return baisclear(iscompound, x);

	if (x % 2 == 0 || x % 3 == 0 || x % 5 == 0 || x % 7 == 0) return 0;

	return is_SPRP(bases_2[(((long long)0xAFF7B4 * x) >> 7) & 1023], (uint64_t)x);
}