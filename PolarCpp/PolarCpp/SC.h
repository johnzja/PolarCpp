#pragma once

#include <cassert>
#include <cmath>
#include "GF.h"

#include <stack>
using namespace std;

typedef unsigned int bits32;
typedef unsigned short bits16;
//typedef vector<bool> BitVec;
//typedef vector<double> LLRVec;
#define ASSERT assert
#define REALMAX 1e20

typedef bool bit;
typedef double LLR;

/* SC-Frame, providing P, C vectors, up/down calculations and left/right propagation functions. 
* Properly using template and class inheritance may simplify coding.
*/
class SCFrame
{
public:
	SCFrame(int N);

	SCFrame(const SCFrame& src);

	LLR left_propagate();

	void right_propagate(bit bit_decision);		

	inline int get_current_index() const;										// return index{u_i}, 0 <= i <= N-1.

	inline void set_index(int input_phi);

	inline void copy_from(const SCFrame& src);

	inline void setup_channel_recv(const LLR* channel_recv);

	virtual ~SCFrame();

protected:
	static void up_calculate(const LLR* llr_x1, const LLR* llr_x2, LLR* result, int len);
	static void down_calculate(const LLR* llr_x1, const LLR* llr_x2, const bit* u2, LLR* result, int len);

	int N, n;
	int phi;
	int* llr_layer_vec, * bit_layer_vec;
	bool copy;			// HERE copy means to share llr_layer_vec and bit_layer_vec.
	LLR* P;
	bit* CL, * CR;
	const LLR* _channel_recv;
	LLR** P_src;		// perform lazy copying.
	bit** CL_src;		// perform lazy copying.
};

/* Successive Cancellation Decoder for binary polar codes. */
class SC_Decoder
{
public:
	SC_Decoder(int N, const bit* frozen_bits, bool binary = true);		// if binary == false, then this object is in fact a q-ary decoder.

	void sc_decode(const LLR* llr, bit* estimated_info_bits);

	int get_K();

	virtual ~SC_Decoder();

protected:
	// Static functions.
	static void up_calculate(const LLR* llr_x1, const LLR* llr_x2, LLR* result, int len);
	static void down_calculate(const LLR* llr_x1, const LLR* llr_x2, const bit* u2, LLR* result, int len);

	int N;
	int n;
	int K;	// K info bits.
	LLR* P;
	bit *CL, *CR;
	int* llr_layer_vec, *bit_layer_vec;
	bit* _frozen_bits;
};

// SC-List implementation: using SC-Frame class.

// q-ary SC decoder.
class SC_Decoder_qary;

class qary_distribution
{
	friend class SC_Decoder_qary;
public:
	qary_distribution(int m);
	~qary_distribution();

	static qary_distribution* newqd(int m, int N);
	static void destroyqd(qary_distribution* pqd, int N);
public:
	double* dist;
public:
	int m;
	int L;
};

class SC_Decoder_qary : public SC_Decoder
{
public:
	SC_Decoder_qary(int N, int m, const bit* frozen_bits, const GF& alpha);					// constructor of fully-frozen decoder.
	SC_Decoder_qary(int N, int m, const GF* frozen_syms, const GF& alpha);				// constructor of partially-frozen decoder.
	void sc_decode_qary(const qary_distribution* probs, bit* estimated_info_bits);

	static qary_distribution* convert_llr_into_qdist(int N_qary, int m, double* llr_arr);
	virtual ~SC_Decoder_qary();

protected:
	static void up_calculate(const qary_distribution* llr_x1, const qary_distribution* llr_x2, qary_distribution* result, GF alpha, int len);
	static void down_calculate(const qary_distribution* llr_x1, const qary_distribution* llr_x2, const GF* u1, qary_distribution* result, GF alpha, int len);
	
	int find_max(const qary_distribution& dist);

protected:
	GF* _qary_frozen_syms;
	qary_distribution* P;
	GF *CL, *CR;
	bool partially_frozen;
	GF alpha;
	int m;
};

/* SCL decoders. */
typedef struct
{
	double x;
	int index;
}PM_with_index;

class SCL_decoder
{
public:
	SCL_decoder(int N, const bit* frozen_bits, int list_size);

	virtual ~SCL_decoder();

	void scl_decode(const LLR* llr, bit* estimated_info_bits);

protected:
	int L;	// list size.
	int N;	// code length.
	int K;	// number of information bits.
	bit* _frozen_bits;
	SCFrame* SCList;

	double* PM, *PM_0, *PM_1;
	bool* is_active, *active_0, *active_1;
	bit** u;
	LLR* ui_llr;
	PM_with_index* pwi;
	stack<int> stk_killed;
};


