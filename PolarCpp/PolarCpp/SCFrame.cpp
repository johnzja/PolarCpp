#include "SCFrame.h"


// MACRO definitions.
#define min(x,y) (((x)<(y))?(x):(y))
#define ASSERT assert

// Const Switches.
const bool use_lazy_copy = true;


/* Binary SC-Frame functions */
SCFrame::SCFrame(int N, bool is_binary) : copy(false), N(N), _channel_recv(NULL), is_binary(is_binary)
{
	ASSERT(!(N & (N - 1)));	// Ensure N is power of 2.
	n = (int)(log2(N));
	phi = 0;

	if (is_binary)
	{
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
	else
	{
		P = NULL; CL = NULL; CR = NULL;
		P_src = NULL; CL_src = NULL;
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
SCFrame::SCFrame(const SCFrame& src) :copy(true), is_binary(src.is_binary)
{
	bit_layer_vec = src.bit_layer_vec;
	llr_layer_vec = src.llr_layer_vec;
	N = src.N;
	n = src.n;
	_channel_recv = src._channel_recv;
	phi = 0;

	if (is_binary)
	{
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
	else
	{
		P = NULL; CL = NULL; CR = NULL; P_src = NULL; CL_src = NULL;
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

	if (phi_mod_2) CR[0] = bit_decision; else CL[0] = bit_decision;

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
	if (use_lazy_copy)
	{
		for (int i = 0; i < n; i++)
		{
			P_src[i] = src.P_src[i];
			CL_src[i] = src.CL_src[i];
		}
	}
	else
	{
		for (int j = 0; j < N - 1; j++)
		{
			P[j] = src.P[j];
			CL[j] = src.CL[j];
		}
	}
}

void SCFrame::setup_channel_recv(const LLR* channel_recv)
{
	_channel_recv = channel_recv;
}

void SCFrame::reset_all()
{
	// Reset operations must be performed when the SCL decoder starts another decoding cycle.
	// Resetting lazy-copy pointers.
	for (int i = 0; i < n; i++)
	{
		P_src[i] = P;
		CL_src[i] = CL;
	}
	phi = 0;
}



/* Q-ary SC-Frame functions */
Qary_SCFrame::Qary_SCFrame(int N, int m, const GF& alpha):SCFrame(N, false), m(m), alpha(alpha)
{
	// Initialize Base Class as is_binary = false, disabling initialization of binary P, CL and CR arrays.
	GF_ASSERT(m >= GF_M_MIN && m <= GF_M_MAX);
	P = qary_distribution::newqd(m, N - 1);
	CL = new GF[N - 1]; CR = new GF[N - 1];

	P_src = new qary_distribution * [n];
	CL_src = new GF * [n];

	for (int i = 0; i < n; i++)
	{
		P_src[i] = P;
		CL_src[i] = CL;		// initialize sources.
	}
}

Qary_SCFrame::Qary_SCFrame(const Qary_SCFrame& src):SCFrame(src), m(src.m), _channel_recv(src._channel_recv), alpha(src.alpha)
{
	P = qary_distribution::newqd(src.m, N - 1);
	CL = new GF[N - 1]; CR = new GF[N - 1];

	P_src = new qary_distribution * [n];
	CL_src = new GF * [n];

	for (int i = 0; i < n; i++)
	{
		P_src[i] = P;
		CL_src[i] = CL;		// initialize sources.
	}
}

Qary_SCFrame::~Qary_SCFrame()
{
	qary_distribution::destroyqd(P, N - 1);
	delete[] CL;
	delete[] CR;

	delete[] P_src;
	delete[] CL_src;
}

void Qary_SCFrame::up_calculate(const qary_distribution* qdist_x1, const qary_distribution* qdist_x2, qary_distribution* result, const GF& alpha, int len)
{
	ASSERT(qdist_x1->m == qdist_x2->m);
	ASSERT(qdist_x1->m == alpha.m);

	int _m = qdist_x1->m;
	int q = (0x1) << _m;

	for (int i = 0; i < len; i++)
	{
		double* result_dist = result[i].dist;
		double* x1_dist = qdist_x1[i].dist;
		double* x2_dist = qdist_x2[i].dist;

		for (short u1 = 0; u1 < q; u1++)
		{
			result_dist[u1] = 0.0;
			for (short u2 = 0; u2 < q; u2++)
			{
				GF r = GF(_m, u1) + alpha * GF(_m, u2);
				result_dist[u1] += x1_dist[r.x] * x2_dist[u2];
			}
		}
	}
}

void Qary_SCFrame::down_calculate(const qary_distribution* qdist_x1, const qary_distribution* qdist_x2, const GF* u1, qary_distribution* result, const GF& alpha, int len)
{
	ASSERT(qdist_x1->m == qdist_x2->m);
	ASSERT(qdist_x1->m == alpha.m);

	int _m = qdist_x1->m;
	int q = (0x1) << _m;

	for (int i = 0; i < len; i++)
	{
		double* result_dist = result[i].dist;
		double* x1_dist = qdist_x1[i].dist;
		double* x2_dist = qdist_x2[i].dist;
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

const qary_distribution& Qary_SCFrame::left_propagate()
{
	// Perform g functions and consecutive f functions.
	int index_1, index_2, op_len;
	if (phi == 0)
	{
		// Use channel LLRs. HERE lazy_copy is not needed, because all the source data come from the physical channel.
		// The same reasoning applies to the second case where phi == N/2.
		index_1 = (0x1) << (n - 1);
		up_calculate(_channel_recv, _channel_recv + index_1, P + index_1 - 1, alpha, index_1);

		for (int i = n - 2; i >= 0; i--)
		{
			index_1 = (0x1) << i;			// 2
			index_2 = index_1 << 1;			// 4
			up_calculate(P + index_2 - 1, P + index_2 - 1 + index_1, P + index_1 - 1, alpha, index_1);
		}
	}
	else if (phi == (N / 2))
	{
		index_1 = (0x1) << (n - 1);
		down_calculate(_channel_recv, _channel_recv + index_1, CL + index_1 - 1, P + index_1 - 1, alpha, index_1);

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
		index_1 = (0x1) << (llr_layer + 1);					// 4
		index_2 = index_1 + ((0x1) << llr_layer);			// 6
		op_len = (0x1) << llr_layer;						// 2

		// Perform g function once. HERE: lazy_copy must be used.
		qary_distribution* P_data_src = P_src[llr_layer + 1];
		down_calculate(P_data_src + index_1 - 1, P_data_src + index_2 - 1, CL + index_1 / 2 - 1, P + index_1 / 2 - 1, alpha, op_len);

		// Peform llr_layer f functions.
		for (int i = llr_layer - 1; i >= 0; i--)
		{
			index_1 = (0x1) << i;			// 1
			index_2 = index_1 << 1;			// 2
			up_calculate(P + index_2 - 1, P + index_2 - 1 + index_1, P + index_1 - 1, alpha, index_1);
		}
	}

	return P[0];
}

void Qary_SCFrame::right_propagate(const GF& gf_decision)
{
	int index_1, index_2;
	int phi_mod_2 = phi & 0x1;

	if (phi_mod_2) CR[0] = gf_decision; else CL[0] = gf_decision;

	if (phi_mod_2 == 1 && phi != (N - 1))			// partial-sum returning.
	{
		int bit_layer = bit_layer_vec[phi];
		for (int i = 0; i < bit_layer; i++)
		{
			index_1 = (0x1) << i;		// 1
			index_2 = index_1 << 1;		// 2

			GF* CL_data_src = CL_src[i];
			// do not call memcpy here.
			for (int start = index_1 - 1; start < index_2 - 1; start++)
			{
				CR[index_1 + start] = CL_data_src[start] + alpha * CR[start];
				CR[index_2 + start] = CR[start];
			}
		}

		index_1 = (0x1) << bit_layer;	// 2
		index_2 = index_1 << 1;			// 4

		GF* CL_data_src = CL_src[bit_layer];
		for (int start = index_1 - 1; start < index_2 - 1; start++)
		{
			CL[index_1 + start] = CL_data_src[start] + alpha * CR[start];
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

void Qary_SCFrame::copy_from(const Qary_SCFrame& src)
{
	phi = src.phi;
	if (use_lazy_copy)
	{
		for (int i = 0; i < n; i++)
		{
			P_src[i] = src.P_src[i];
			CL_src[i] = src.CL_src[i];
		}
	}
	else
	{
		// copy distribution data.
		for (int j = 0; j < N - 1; j++)
		{
			P[j] = src.P[j];
			CL[j] = src.CL[j];
		}
	}
}

void Qary_SCFrame::setup_channel_recv(const qary_distribution* channel_recv)
{
	_channel_recv = channel_recv;
}

void Qary_SCFrame::reset_all()
{
	for (int i = 0; i < n; i++)
	{
		P_src[i] = P;
		CL_src[i] = CL;
	}
	phi = 0;
}
