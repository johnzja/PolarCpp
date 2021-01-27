#include "SC.h"

using namespace std;

#define min(x,y) ((x<y)?(x):(y))

SC_Decoder::SC_Decoder(int N, const bit* frozen_bits):N(N)
{
	ASSERT(!(N&(N - 1)));	// Ensure N is power of 2.
	n = log2(N);
	K = 0;
	// Step1: Count K.
	_frozen_bits = new bit[N];
	for (int i = 0; i < N; i++)
	{
		if (!frozen_bits[i]) K++;
		_frozen_bits[i] = frozen_bits[i];
	}

	P = new LLR[N - 1];
	CL = new bit[N - 1];
	CR = new bit[N - 1];


	// Step2: Generate bit_layer_vec and llr_layer_vec vectors.
	llr_layer_vec = new int[N];
	llr_layer_vec[0] = 0;
	for (int i = 1; i < N; i++)
	{
		// Count the consecutive 0's from the LSB of the binary expansion representation of i.
		int layer = 0;
		int num = i;
		while ((num & 0x1) == 0)
		{
			layer++;
			num >>= 1;
		}
		llr_layer_vec[i] = layer;
	}
	// llr_layer_vec[i] means the layer index from which to read LLR when decoding u[i].

	bit_layer_vec = new int[N];
	for (int i = 0; i < N; i++)
	{
		int num = i >> 1;
		int layer = 0;
		while ((num & 0x1) == 1)
		{
			layer++;
			num >>= 1;
		}
		bit_layer_vec[i] = layer;
	}
	// bit_layer_vec[i] means the number of consecutive 1's  - 1of an odd number.
	// Odd numbers represent leaf nodes that need partial-sum return.

}

SC_Decoder::~SC_Decoder()
{
	delete[] P;
	delete[] CL;
	delete[] CR;
	delete[] _frozen_bits;
	delete[] llr_layer_vec;
	delete[] bit_layer_vec;
}

void SC_Decoder::up_calculate(const LLR* llr_x1, const LLR* llr_x2, LLR* result, int len)
{
	for (int i = 0; i < len; i++)
	{
		double L1 = llr_x1[i], L2 = llr_x2[i];
		bool s1 = (L1 >= 0), s2 = (L2 >= 0);
		double Lans = min(abs(L1), abs(L2));
		result[i] = (s1^s2) ? -Lans : Lans;
	}
}

void SC_Decoder::down_calculate(const LLR* llr_x1, const LLR* llr_x2, const bit* u2, LLR* result, int len)
{
	for (int i = 0; i < len; i++)
	{
		double L1 = llr_x1[i], L2 = llr_x2[i];
		result[i] = (u2[i]) ? (L2 - L1) : (L2 + L1);
	}
}

// {0, 1,2, 3,4,5,6}.
void SC_Decoder::sc_decode(const LLR* llr, bit* estimated_info_bits)
{
	// Assume the array "estimated_info_bits" has length K.
	int index_1, index_2, op_len;
	int k = 0;

	for (int phi = 0; phi < N; phi++)
	{
		if (phi == 0)
		{
			// Use channel LLRs.
			index_1 = (0x1) << (n - 1);
			up_calculate(llr, llr + index_1, P + index_1 - 1, index_1);

			for (int i = n - 2; i >= 0; i--)
			{
				index_1 = (0x1) << i;				// 2
				index_2 = index_1<<1;			// 4
				up_calculate(P + index_2 - 1, P + index_2 - 1 + index_1, P + index_1 - 1, index_1);
			}
		}
		else if (phi == (N / 2))
		{
			index_1 = (0x1) << (n - 1);
			down_calculate(llr, llr + index_1, CL + index_1 - 1, P + index_1 - 1, index_1);

			for (int i = n - 2; i >= 0; i--)
			{
				index_1 = (0x1) << i;				
				index_2 = index_1 << 1;
				up_calculate(P + index_2 - 1, P + index_2 - 1 + index_1, P + index_1 - 1, index_1);
			}
		}
		else
		{
			int llr_layer = llr_layer_vec[phi];
			index_1 = (0x1) << (llr_layer + 1);				// 4
			index_2 = index_1 + ((0x1) << llr_layer);		// 6
			op_len = (0x1) << llr_layer;							// 2

			// Perform g function once.
			down_calculate(P + index_1 - 1, P + index_2 - 1, CL + index_1 / 2 - 1, P + index_1 / 2 - 1, op_len);

			// Peform llr_layer f functions.
			for (int i = llr_layer - 1; i >= 0; i--)
			{
				index_1 = (0x1) << i;			// 1
				index_2 = index_1 << 1;		// 2
				up_calculate(P + index_2 - 1, P + index_2 - 1 + index_1, P + index_1 - 1, index_1);
			}
		}

		// down to the leaf node HERE.
		int phi_mod_2 = phi & 0x1;
		if (_frozen_bits[phi])
		{
			if (phi_mod_2 == 0)CL[0] = 0;
			else CR[0] = 0;
		}
		else
		{
			// non-frozen, info bit.
			estimated_info_bits[k] = (P[0] < 0);
			if (phi_mod_2 == 0)
				CL[0] = estimated_info_bits[k];
			else
				CR[0] = estimated_info_bits[k];

			k++;
		}

		// partial-sum return.
		if (phi_mod_2 == 1 && phi != (N - 1))
		{
			int bit_layer = bit_layer_vec[phi];
			for (int i = 0; i < bit_layer; i++)
			{
				index_1 = (0x1) << i;		// 1
				index_2 = index_1 << 1;	// 2
				
				// do not call memcpy here.
				for (int start = index_1 - 1; start < index_2 - 1; start++)
				{
					CR[index_1 + start] = CL[start] ^ CR[start];
					CR[index_2 + start] = CR[start];
				}
			}

			index_1 = (0x1) << bit_layer;	// 2
			index_2 = index_1 << 1;			// 4

			for (int start = index_1 - 1; start < index_2 - 1; start++)
			{
				CL[index_1 + start] = CL[start] ^ CR[start];
				CL[index_2 + start] = CR[start];
			}
		}
	}
}

// TODO: Adapt "bit" in the code into any q-ary symbol.

/*           q-ary probabilities                 */
qary_distribution::qary_distribution(int m):m(m)
{
	GF_ASSERT(m >= GF_M_MIN && m <= GF_M_MAX);
	L = (0x1) << m;
	dist = new double[L];
}

qary_distribution::~qary_distribution()
{
	delete[] dist;
}

void SC_Decoder_qary::up_calculate(const qary_distribution* llr_x1, const qary_distribution* llr_x2, qary_distribution* result, GF alpha, int len)
{
	ASSERT(llr_x1->m == llr_x2->m);
	ASSERT(llr_x1->m == alpha.m);

	int _m = llr_x1->m;
	int q = (0x1) << _m;

	for (int i = 0; i < len; i++)
	{
		double* result_dist = result[i].dist;
		double* x1_dist = llr_x1[i].dist;
		double* x2_dist = llr_x2[i].dist;

		for (short u1 = 0; u1 < q; u1++)
		{
			result_dist[u1] = 0.0;
			for (short u2 = 0; u2 < q; u2++)
			{
				GF r = GF(_m, u1) + alpha * GF(_m, u2);
				result_dist[u1] += x1_dist[r.x] * x2_dist[u2];
			}
			
			// normalization may not be useful.
			result_dist[u1] /= q;
		}
	}
}

void SC_Decoder_qary::down_calculate(const qary_distribution* llr_x1, const qary_distribution* llr_x2, const GF* u1, qary_distribution* result, GF alpha, int len)
{
	ASSERT(llr_x1->m == llr_x2->m);
	ASSERT(llr_x1->m == alpha.m);

	int _m = llr_x1->m;
	int q = (0x1) << _m;

	for (int i = 0; i < len; i++)
	{
		double* result_dist = result[i].dist;
		double* x1_dist = llr_x1[i].dist;
		double* x2_dist = llr_x2[i].dist;
		const GF& u1_i = u1[i];

		for (short u2 = 0; u2 < q; u2++)
		{
			GF r = u1_i + alpha * GF(_m, u2);
			result_dist[u2] = x1_dist[r.x] * x2_dist[u2];
			result_dist[u2] /= q;
		}

	}
}


