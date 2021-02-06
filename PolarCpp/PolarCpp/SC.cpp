#include <intrin.h>
#include "SC.h"

using namespace std;

#define min(x,y) ((x<y)?(x):(y))


/* SC-Frame functions */
SCFrame::SCFrame(int N) : copy(false), N(N)
{
	ASSERT(!(N & (N - 1)));	// Ensure N is power of 2.
	n = (int)(log2(N));
	phi = 0;

	P = new LLR[N - 1];
	CL = new bit[N - 1];
	CR = new bit[N - 1];

	P_src = new LLR * [n];
	CL_src = new bit * [n];
	for (int i = 0; i < n; i++)
	{
		P_src[i] = P;
		CL_src[i] = CL;	// initialize sources.
	}

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

// copy constructor, used in initializing.
SCFrame::SCFrame(const SCFrame& src) :copy(true)
{
	bit_layer_vec = src.bit_layer_vec;
	llr_layer_vec = src.llr_layer_vec;
	N = src.N;
	n = src.n;
	phi = 0;		
	_channel_recv = src._channel_recv;

	P = new LLR[N - 1];
	CL = new bit[N - 1];
	CR = new bit[N - 1];

	P_src = new LLR * [n];
	CL_src = new bit * [n];
	for (int i = 0; i < n; i++)
	{
		P_src[i] = P;
		CL_src[i] = CL;	// initialize sources.
	}

}

SCFrame::~SCFrame()
{
	if (!copy)
	{
		delete[] llr_layer_vec;
		delete[] bit_layer_vec;
	}
	delete[] P;
	delete[] CL;
	delete[] CR;
	delete[] P_src;
	delete[] CL_src;
}


void SCFrame::up_calculate(const LLR* llr_x1, const LLR* llr_x2, LLR* result, int len)
{
	for (int i = 0; i < len; i++)
	{
		double L1 = llr_x1[i], L2 = llr_x2[i];
		bool s1 = (L1 >= 0), s2 = (L2 >= 0);
		double Lans = min(abs(L1), abs(L2));
		result[i] = (s1 ^ s2) ? -Lans : Lans;
	}
}

void SCFrame::down_calculate(const LLR* llr_x1, const LLR* llr_x2, const bit* u2, LLR* result, int len)
{
	for (int i = 0; i < len; i++)
	{
		double L1 = llr_x1[i], L2 = llr_x2[i];
		result[i] = (u2[i]) ? (L2 - L1) : (L2 + L1);
	}
}

int SCFrame::get_current_index() const
{
	return phi;
}

void SCFrame::reset_index()
{
	phi = 0;
}

LLR SCFrame::left_propagate()
{
	// Perform g functions and consecutive f functions.
	int index_1, index_2, op_len;
	if (phi == 0)
	{
		// Use channel LLRs. HERE lazy_copy is not needed, because all the source data come from the physical channel.
		// The same reasoning applies to the second case where phi == N/2.
		index_1 = (0x1) << (n - 1);
		up_calculate(_channel_recv, _channel_recv + index_1, P + index_1 - 1, index_1);

		for (int i = n - 2; i >= 0; i--)
		{
			index_1 = (0x1) << i;			// 2
			index_2 = index_1 << 1;			// 4
			up_calculate(P + index_2 - 1, P + index_2 - 1 + index_1, P + index_1 - 1, index_1);
		}
	}
	else if (phi == (N / 2))
	{
		index_1 = (0x1) << (n - 1);
		down_calculate(_channel_recv, _channel_recv + index_1, CL + index_1 - 1, P + index_1 - 1, index_1);

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
		index_1 = (0x1) << (llr_layer + 1);					// 4
		index_2 = index_1 + ((0x1) << llr_layer);			// 6
		op_len = (0x1) << llr_layer;						// 2

		// Perform g function once. HERE: lazy_copy must be used.
		LLR* P_data_src = P_src[llr_layer + 1];
		down_calculate(P_data_src + index_1 - 1, P_data_src + index_2 - 1, CL + index_1 / 2 - 1, P + index_1 / 2 - 1, op_len);

		// Peform llr_layer f functions.
		for (int i = llr_layer - 1; i >= 0; i--)
		{
			index_1 = (0x1) << i;			// 1
			index_2 = index_1 << 1;			// 2
			up_calculate(P + index_2 - 1, P + index_2 - 1 + index_1, P + index_1 - 1, index_1);
		}
	}

	return P[0];
}

void SCFrame::right_propagate(bit bit_decision)
{
	int index_1, index_2;
	int phi_mod_2 = phi & 0x1;

	if (phi_mod_2) CL[0] = bit_decision; else CR[0] = bit_decision;		

	if (phi_mod_2 == 1 && phi != (N - 1))			// partial-sum returning.
	{
		int bit_layer = bit_layer_vec[phi];
		for (int i = 0; i < bit_layer; i++)
		{
			index_1 = (0x1) << i;		// 1
			index_2 = index_1 << 1;		// 2

			bit* CL_data_src = CL_src[i];
			// do not call memcpy here.
			for (int start = index_1 - 1; start < index_2 - 1; start++)
			{
				CR[index_1 + start] = CL_data_src[start] ^ CR[start];
				CR[index_2 + start] = CR[start];
			}
		}

		index_1 = (0x1) << bit_layer;	// 2
		index_2 = index_1 << 1;			// 4

		bit* CL_data_src = CL_src[bit_layer];
		for (int start = index_1 - 1; start < index_2 - 1; start++)
		{
			CL[index_1 + start] = CL_data_src[start] ^ CR[start];
			CL[index_2 + start] = CR[start];
		}
	}

	if (phi < (N - 1))
	{
		for (int i_layer = 0; i_layer <= llr_layer_vec[phi + 1]; i_layer++)
		{
			P_src[i_layer] = P;
			CL_src[i_layer] = CL;		// Prepare for the next P-array-related g function call.
		}
	}

	phi++;
}

void SCFrame::copy_from(const SCFrame& src)
{
	phi = src.phi;
	for (int i = 0; i < n; i++)
	{
		P_src[i] = src.P_src[i];
		CL_src[i] = src.CL_src[i];
	}
}

void SCFrame::setup_channel_recv(const LLR* channel_recv)
{
	_channel_recv = channel_recv;
}

/* Decoders */
SC_Decoder::SC_Decoder(int N, const bit* frozen_bits, bool binary):N(N)
{
	ASSERT(!(N&(N - 1)));	// Ensure N is power of 2.
	n = log2(N);
	K = 0;

	// Step1: Count K.
	if (frozen_bits)
	{
		_frozen_bits = new bit[N];
		for (int i = 0; i < N; i++)
		{
			if (!frozen_bits[i]) K++;
			_frozen_bits[i] = frozen_bits[i];
		}
	}
	else
		_frozen_bits = NULL;

	if (binary)
	{
		P = new LLR[N - 1];
		CL = new bit[N - 1];
		CR = new bit[N - 1];
	}	// Notice: if not binary, then P, CL and CR should be re-initialized.
	else
	{
		P = NULL; CL = NULL; CR = NULL; // delete a nullptr is allowed.
	}

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

int SC_Decoder::get_K()
{
	return K;
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
// Posteriori probability distribution on GF(q).
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

qary_distribution* qary_distribution::newqd(int m, int N)
{
	// generate qary_distribution objects.
	qary_distribution* y = (qary_distribution*) operator new(N * sizeof(qary_distribution));
	//qary_distribution* y = new qary_distribution[N];
	for (int i = 0; i < N; i++)
	{
		new (y + i)qary_distribution(m);	// initialize, using "placement new".
	}
	return y;
}

void qary_distribution::destroyqd(qary_distribution* pqd, int N)
{
	for (int i = 0; i < N; i++)
	{
		(pqd + i)->~qary_distribution();
	}
	operator delete ((void*)pqd);		// operator new.
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
			// result_dist[u1] /= q;
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

		double sum = 0;
		for (short u2 = 0; u2 < q; u2++)
		{
			GF r = u1_i + alpha * GF(_m, u2);
			sum += (result_dist[u2] = x1_dist[r.x] * x2_dist[u2]);
		}

		for (short u2 = 0; u2 < q; u2++)
			result_dist[u2] /= sum;
	}
}


SC_Decoder_qary::SC_Decoder_qary(int N, int m, const bit* frozen_bits, const GF& alpha):m(m), SC_Decoder(N, frozen_bits, false), alpha(alpha)
{
	GF_ASSERT(m >= GF_M_MIN && m <= GF_M_MAX);
	partially_frozen = false;

	// _frozen_bits and the number K is available in SC decoder.
	_qary_frozen_syms = NULL;
	K *= m;									// Return bits.

	// construct P, CL and CR.
	P = qary_distribution::newqd(m, N-1);
	CL = new GF[N - 1];
	CR = new GF[N - 1];
}

SC_Decoder_qary::SC_Decoder_qary(int N, int m, const GF* frozen_syms, const GF& alpha) :m(m), SC_Decoder(N, NULL, false), alpha(alpha)
{
	GF_ASSERT(m >= GF_M_MIN && m <= GF_M_MAX);
	partially_frozen = true;

	// count K.
	_qary_frozen_syms = new GF[N];
	K = 0;
	for (int i = 0; i < N; i++)
	{
		_qary_frozen_syms[i] = frozen_syms[i];
		K += (m - __popcnt16(frozen_syms[i].x));		// count the number of 1's.
	}

	P = qary_distribution::newqd(m, N - 1);
	CL = new GF[N - 1];
	CR = new GF[N - 1];
}

SC_Decoder_qary:: ~SC_Decoder_qary()
{
	delete[] _qary_frozen_syms;
	qary_distribution::destroyqd(P, N - 1);
	delete[] CL;
	delete[] CR;

	// Automatically call ~SC_Decoder.
}

int SC_Decoder_qary::find_max(const qary_distribution& dist)
{
	double* prob_dist = dist.dist;
	int index = 0;
	double p = prob_dist[0];

	for (int i = 1; i < ((0x1) << m); i++)
	{
		if (prob_dist[i] > p)
		{
			p = prob_dist[i];
			index = i;
		}
	}
	return index;
}

void SC_Decoder_qary::sc_decode_qary(const qary_distribution* probs, bit* estimated_info_bits)
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
			up_calculate(probs, probs + index_1, P + index_1 - 1, alpha, index_1);

			for (int i = n - 2; i >= 0; i--)
			{
				index_1 = (0x1) << i;				// 2
				index_2 = index_1 << 1;			// 4
				up_calculate(P + index_2 - 1, P + index_2 - 1 + index_1, P + index_1 - 1, alpha, index_1);
			}
		}
		else if (phi == (N / 2))
		{
			index_1 = (0x1) << (n - 1);
			down_calculate(probs, probs + index_1, CL + index_1 - 1, P + index_1 - 1, alpha, index_1);

			for (int i = n - 2; i >= 0; i--)
			{
				index_1 = (0x1) << i;
				index_2 = index_1 << 1;
				up_calculate(P + index_2 - 1, P + index_2 - 1 + index_1, P + index_1 - 1, alpha, index_1);
			}
		}
		else
		{
			int llr_layer = llr_layer_vec[phi];
			index_1 = (0x1) << (llr_layer + 1);				// 4
			index_2 = index_1 + ((0x1) << llr_layer);		// 6
			op_len = (0x1) << llr_layer;							// 2

			// Perform g function once.
			down_calculate(P + index_1 - 1, P + index_2 - 1, CL + index_1 / 2 - 1, P + index_1 / 2 - 1, alpha, op_len);

			// Peform llr_layer f functions.
			for (int i = llr_layer - 1; i >= 0; i--)
			{
				index_1 = (0x1) << i;			// 1
				index_2 = index_1 << 1;		// 2
				up_calculate(P + index_2 - 1, P + index_2 - 1 + index_1, P + index_1 - 1, alpha, index_1);
			}
		}

		// down to the leaf node HERE.
		// HERE it is a q-ary leaf node.

		int phi_mod_2 = phi & 0x1;

		if (!partially_frozen)
		{
			if (_frozen_bits[phi])
			{
				if (phi_mod_2 == 0)CL[0] = GF(m, 0);
				else CR[0] = GF(m, 0);
			}
			else
			{
				// non-frozen, info symbol.
				// find argmax.
				GF hard_decision = GF(m, find_max(P[0]));
				short temp_x = hard_decision.x;
				for (int j = k; j < k + m; j++)
				{
					estimated_info_bits[j] = bool(temp_x & 0x1);
					temp_x >>= 1;
				}
				if (phi_mod_2 == 0)
					CL[0] = hard_decision;
				else
					CR[0] = hard_decision;

				k += m;
			}
		}
		else
		{
			// partially frozen.
			short partially_frozen_pattern = _qary_frozen_syms[phi].x;
			int one_cnt = __popcnt16(partially_frozen_pattern);
			if (m == one_cnt)
			{
				// completely frozen.
				if (phi_mod_2 == 0)CL[0] = GF(m, 0);
				else CR[0] = GF(m, 0);
			}
			else
			{
				double* P0_dist = P[0].dist;
				for (short i = 0; i < (0x1 << m); i++)
				{
					if (partially_frozen_pattern & i)
					{
						P0_dist[i] = -0.0;	// impossible.
					}
				}

				GF hard_decision = GF(m, find_max(P[0]));
				short temp_x = hard_decision.x;
				for (int i = 0; i < m; i++)
				{
					if ((~partially_frozen_pattern) & (0x1 << i))
					{
						estimated_info_bits[k++] = bool(temp_x & 0x1);
					}
					temp_x >>= 1;
				}

				if (phi_mod_2 == 0)
					CL[0] = hard_decision;
				else
					CR[0] = hard_decision;
			}
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
					CR[index_1 + start] = CL[start] + alpha * CR[start];
					CR[index_2 + start] = CR[start];
				}
			}

			index_1 = (0x1) << bit_layer;	// 2
			index_2 = index_1 << 1;			// 4

			for (int start = index_1 - 1; start < index_2 - 1; start++)
			{
				CL[index_1 + start] = CL[start] + alpha * CR[start];
				CL[index_2 + start] = CR[start];
			}
		}
	}
}

// toolkits.
qary_distribution* SC_Decoder_qary::convert_llr_into_qdist(int N_qary, int m, double* llr_arr)
{
	// llr_arr has length N_qary * m.
	qary_distribution* result = qary_distribution::newqd(m, N_qary);

	double* bpsk_posteriori_0 = new double[m];
	double* bpsk_posteriori_1 = new double[m];

	for (int i = 0; i < N_qary; i++)
	{
		for (int j = 0; j < m; j++)
		{
			double lr = exp(llr_arr[i*m + j]);
			bpsk_posteriori_1[j] = 1 / (1 + lr);
			bpsk_posteriori_0[j] = 1 - bpsk_posteriori_1[j];
		}

		qary_distribution& qd = result[i];
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
	return result;
}


/* SCL decoder, can be independent of the above SC and q-ary SC decoders. */
SCL_decoder::SCL_decoder(int N, const bit* frozen_bits, int list_size):L(list_size), N(N)
{
	ASSERT((N & (N - 1)) == 0);
	ASSERT((L & (L - 1)) == 0);	// Ensure that L is a power of two.

	if (frozen_bits)
	{
		_frozen_bits = new bit[N];
		for (int i = 0; i < N; i++)
		{
			if (!frozen_bits[i]) K++;
			_frozen_bits[i] = frozen_bits[i];
		}
	}
	else
		_frozen_bits = NULL;

	// Construct SC list.
	// using operator new.
	SCList = (SCFrame*) operator new(L * sizeof(SCFrame));
	new (SCList + 0) SCFrame(N);
	for (int i = 1; i < L; i++)
	{
		new (SCList + i) SCFrame(*(SCList));
	}

}

SCL_decoder::~SCL_decoder()
{
	for (int i = 0; i < L; i++)
	{
		(SCList + i)->~SCFrame();		// Call the destructor explicitly.
	}
	operator delete((void*)SCList);
}

void SCL_decoder::scl_decode(const LLR* llr, bit* estimated_info_bits)
{
	// Step1: Load channel LLRs.
	double* PM = new double[L];
	bool* is_active = new bool[L];
	bit** u = new bit * [L];

	for (int l_index = 0; l_index < L; l_index++)
	{
		SCList[l_index].setup_channel_recv(llr);
		SCList[l_index].reset_index();
		PM[l_index] = 0.0;
		is_active[l_index] = false;
		u[l_index] = new bit[K];
	}
	is_active[0] = true;

	int k = 0;
	// Step2: iteratively decode each bits.
	for (int phi = 0; phi < N; phi++)
	{
		for (int l_index = 0; l_index < L; l_index++)
		{
			if (!is_active[l_index]) continue;
			double ui_llr = SCList[l_index].left_propagate();

			bit bit_decision;
			if (_frozen_bits[phi])
				bit_decision = false;
			else
				bit_decision = (ui_llr < 0);
		}
	}
	// TODO: unfinished.


	delete[] PM;
	delete[] is_active;
	for (int l_index = 0; l_index < L; l_index++) delete[] u[l_index];
	delete[] u;
}



