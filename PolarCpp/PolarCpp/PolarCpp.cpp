// PolarCpp.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include "SC.h"

#include <vector>
#include <string>
#include <random>
#include <algorithm>
#include <cmath>

using namespace std;

void polar_encode(const bit* u, bit* x, int N, const bit* frozen_bits)
{
	ASSERT(!(N&(N - 1)));
	int n = log2(N);
	
	int k = 0;
	for (int i = 0; i < N; i++)
	{
		if (!frozen_bits[i])
		{
			x[i] = u[k++];
		}
		else
		{
			x[i] = false;
		}
	}

	// x=uF^n. Calculate on-site.
	for (int layer = 1; layer <= n; layer++)
	{
		int step = (0x1) << layer;	// start from 2.
		int inner = step >> 1;
		int blocks = N >> layer;

		for (int b = 0; b < blocks; b++)
		{
			int block_up = b * step;
			int block_down = block_up + inner;
			for (int i = 0; i < inner; i++)
			{
				x[block_up + i] = x[block_up + i] ^ x[block_down + i];
				x[block_down + i] = x[block_down + i];
			}
		}
	}
}

int main()
{
	std::default_random_engine e;
	std::normal_distribution<double> n(0, 1);

	// Study polar codes using normal distribution generator.
	bit frozen_bits[] = { 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,1,1,1,1,0,0,0,1,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,0,1,0,0,0,1,1,1,0,1,0,0,0,1,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,1,1,1,1,1,1,1,0,1,1,1,0,1,0,0,0,1,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,0,1,1,1,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, };

	int N = 512, M = 256;
	int N_sim = 5000;

	double Ebn0 = 2.5;
	double sigma = 1 / sqrt(2 * (M / ((double)(N))))*pow(10, -Ebn0 / 20);
	cout << "Evaluate SC @ Eb/n0 = " << Ebn0 << endl;

	bit* bits_to_encode = new bit[M];
	bit* bits_encoded = new bit[N];
	bit* bits_decoded = new bit[M];
	LLR* y = new LLR[N];

	// Construct SC decoder.
	SC_Decoder scd(N, frozen_bits);

	int block_error_cnt = 0;

	/* START simulation */
	for (int sim_iter = 0; sim_iter < N_sim; sim_iter++)
	{
		// Step1: Generate random bits to be encoded.
		for (int i = 0; i < M; i++)
			bits_to_encode[i] = (rand() & 0x1);

		// Step2: Polar-encode.
		polar_encode(bits_to_encode, bits_encoded, N, frozen_bits);

		// Step3: Add noise.
		for (int i = 0; i < N; i++)
		{
			y[i] = 1 - 2 * bits_encoded[i];		// BPSK
			y[i] += sigma * n(e);					// Add noise
			y[i] = 2 * y[i] / (sigma*sigma);		// Calculate LLR.
		}

		// Step4: Perform SC-decoding.
		scd.sc_decode(y, bits_decoded);

		// Step5: Display number of error bits.
		
		for (int i = 0; i < M; i++)
		{
			if (bits_decoded[i] != bits_to_encode[i])
			{
				block_error_cnt++;
				break;
			}
		}
	}

	cout << "BLER = " << ((double)block_error_cnt) / N_sim << endl;

	delete[] bits_to_encode;
	delete[] bits_encoded;
	delete[] y;
	delete[] bits_decoded;
	cout << "Hello SC!" << endl;
	return 0;
}

