// PolarCpp.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include "SC.h"

#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <time.h>

using namespace std;

void polar_encode(const bit* d, bit* x, int N, const bit* frozen_bits)
{
	ASSERT(!(N&(N - 1)));
	int n = log2(N);
	
	int k = 0;
	for (int i = 0; i < N; i++)
	{
		if (!frozen_bits[i])
		{
			x[i] = d[k++];
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
				//x[block_down + i] = x[block_down + i];
			}
		}
	}
}

void qary_polar_encode(const bit* d, GF* x, int N, int m, const bit* frozen_bits, const GF& alpha)
{
	// input: bitstream, output: GF symbol stream. Fully frozen.
	// m bits once for encoding.
	ASSERT(!(N&(N - 1)));
	int n = log2(N);

	int k = 0;
	int N_info_syms = 0;
	for (int i = 0; i < N; i++)
	{
		if (!frozen_bits[i])
		{
			x[i] = GF(m, 0);
			for (int j = 0; j < m; j++)
			{
				if (d[k++])
				{
					x[i].x |= (0x1 << j);
				}
			}
			N_info_syms++;
		}
		else
		{
			x[i] = GF(m, 0);			// frozen.
		}
	}

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
				x[block_up + i] = x[block_up + i] + alpha * x[block_down + i];
				x[block_down + i] = x[block_down + i];				// may be unuseful.
			}
		}
	}
}

void qary_polar_encode(const bit* d, GF* x, int N, int m, const GF* frozen_syms, const GF& alpha)
{
	// input: bitstream, output: GF symbol stream. Partially frozen.
	ASSERT(!(N&(N - 1)));
	int n = log2(N);

	int k = 0;
	int N_info_syms = 0;

	// Initialize.
	for (int i = 0; i < N; i++)
	{
		const GF& fs = frozen_syms[i];	// frozen sym indicator.
		GF& xi = x[i];
		xi = GF(m, 0);

		for (int j = 0; j < m; j++)
		{
			if (!(fs.x & ((0x1) << j)))			// if not frozen:
			{
				if (d[k++])
				{
					xi.x |= ((0x1) << j);
				}
			}
		}
	}

	for (int layer = 1; layer <= n; layer++)
	{
		int step = (0x1) << layer;			// start from 2.
		int inner = step >> 1;
		int blocks = N >> layer;

		for (int b = 0; b < blocks; b++)
		{
			int block_up = b * step;
			int block_down = block_up + inner;
			for (int i = 0; i < inner; i++)
			{
				x[block_up + i] = x[block_up + i] + alpha * x[block_down + i];
			}
		}
	}

}

void qary_modem_bpsk(const GF* x, qary_distribution* y, int N, double sigma_channel_bpsk, default_random_engine& e, normal_distribution<double>& n)
{
	// using independent BPSK.
	// N: number of GF symbols.
	int m = x[0].m;
	GF_ASSERT(m == y[0].m);
	double* bpsk_posteriori_0 = new double[m];
	double* bpsk_posteriori_1 = new double[m];
	double a = 2.0 / (sigma_channel_bpsk*sigma_channel_bpsk);

	for (int i = 0; i < N; i++)
	{
		const GF& xi = x[i];
		double s;
		for (int j = 0; j < m; j++)
		{
			s = (xi.x & ((0x1) << j)) ? -1.0 : 1.0;
			s += sigma_channel_bpsk * n(e);

			double lr = exp(a*s);
			bpsk_posteriori_1[j] = 1 / (1 + lr);
			bpsk_posteriori_0[j] = 1 - bpsk_posteriori_1[j];
		}

		qary_distribution& qd = y[i];
		for (int k = 0; k < (0x1 << m); k++)qd.dist[k] = 1.0;	// initialize as 1.0 to enable multiplications.
		for (int j = 0; j < m; j++)
		{
			for (int k = 0; k < (0x1 << m); k++)
			{
				if (k & (0x1 << j))
				{
					qd.dist[k] *= bpsk_posteriori_1[j];
				}
				else
				{
					qd.dist[k] *= bpsk_posteriori_0[j];
				}
			}
		}
	}

	delete[] bpsk_posteriori_0;
	delete[] bpsk_posteriori_1;
}

double test_SC(int min_errors = 800)
{
	std::default_random_engine e;
	std::normal_distribution<double> n(0, 1);

	// Study polar codes using normal distribution generator.
	// bit frozen_bits[] = { 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,1,1,1,1,0,0,0,1,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,0,1,0,0,0,1,1,1,0,1,0,0,0,1,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,1,1,1,1,1,1,1,0,1,1,1,0,1,0,0,0,1,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,0,1,1,1,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, };
	
	bit frozen_bits[] = {1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0};
	int N = 16, M = 8;

	ASSERT(sizeof(frozen_bits) == N);

	double Ebn0 = 2.5;
	double sigma = 1 / sqrt(2 * (M / ((double)(N))))*pow(10, -Ebn0 / 20);
	cout << "Evaluate SC @ Eb/n0 = " << Ebn0 << endl;
	cout << "Code configuration: (" << N << ", " << M << ") binary Polar code." << endl;

	bit* bits_to_encode = new bit[M];
	bit* bits_encoded = new bit[N];
	bit* bits_decoded = new bit[M];
	LLR* y = new LLR[N];

	// Construct SC decoder.
	SC_Decoder scd(N, frozen_bits);

	int block_error_cnt = 0;
	int N_runs = 0;

	/* START simulation */
	while(block_error_cnt < min_errors)
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

		// Step5: Count number of error bits.
		for (int i = 0; i < M; i++)
		{
			if (bits_decoded[i] != bits_to_encode[i])
			{
				block_error_cnt++;
				break;
			}
		}
		N_runs++;
	}

	double r;
	cout << "BLER = " << (r = ((double)block_error_cnt) / N_runs) << endl;

	delete[] bits_to_encode;
	delete[] bits_encoded;
	delete[] y;
	delete[] bits_decoded;

	cout << "Test SC complete!" << endl;
	return r;
}

double test_qary_SC(int min_errors = 800)
{
	std::default_random_engine e;
	std::normal_distribution<double> n(0, 1);

	// Study polar codes using normal distribution generator.
	// Use (8,4) 4-ary code with partially frozen.
	GF frozen_bits[] = { GF(2,3), GF(2,3), GF(2,3), GF(2,0), GF(2,3), GF(2,0), GF(2,0), GF(2,0) };

	int N = 8, M = 4;
	int N_sim = 10000;
	int m = 2;						// 4-ary code.

	double R = 0.5;			// Code rate = 0.5
	double Ebn0 = 2.5;
	double sigma = 1 / sqrt(2 * R)*pow(10, -Ebn0 / 20);
	cout << "Evaluate 4-ary SC @ Eb/n0 = " << Ebn0 << endl;
	cout << "Code configuration: (" << N << ", " << M << ") quaternary Polar code" << endl;

	int N_bits = N * m;
	int M_bits = 8;

	bit* bits_to_encode = new bit[M_bits];
	GF* syms_encoded = new GF[N];
	bit* bits_decoded = new bit[M_bits];

	// generate qary_distribution objects.
	qary_distribution* y = qary_distribution::newqd(m, N);

	// Find primitive alpha.
	GF t(m, 0);
	GF alpha = t.get_prim();

	// Construct q-ary SC decoder.
	SC_Decoder_qary scd(N, m, frozen_bits, alpha);

	int block_error_cnt = 0;
	int N_runs = 0;

	/* START simulation */
	while(block_error_cnt < min_errors)
	{
		// Step1: Generate random bits to be encoded.
		for (int i = 0; i < M_bits; i++)
			bits_to_encode[i] = (rand() & 0x1);

		// Step2: Polar-encode.
		//polar_encode(bits_to_encode, bits_encoded, N, frozen_bits);
		qary_polar_encode(bits_to_encode, syms_encoded, N, m, frozen_bits, alpha);

		// Step3: Add noise.
		qary_modem_bpsk(syms_encoded, y, N, sigma, e, n);

		// Step4: Perform SC-decoding.
		//scd.sc_decode(y, bits_decoded);
		scd.sc_decode_qary(y, bits_decoded);

		// Step5: Count number of error blocks.
		for (int i = 0; i < M_bits; i++)
		{
			if (bits_decoded[i] != bits_to_encode[i])
			{
				block_error_cnt++;
				break;
			}
		}
		N_runs++;
	}

	double r;
	cout << "BLER = " << (r = ((double)block_error_cnt) / N_runs) << endl;

	delete[] bits_to_encode;
	delete[] syms_encoded;
	//delete[] y;
	qary_distribution::destroyqd(y);
	delete[] bits_decoded;

	cout << "Test q-ary SC complete!" << endl;
	return r;
}

int main()
{
	test_SC(1600);
	cout << endl;
	test_qary_SC(1600);
	return 0;
}

