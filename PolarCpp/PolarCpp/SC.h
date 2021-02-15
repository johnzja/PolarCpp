#pragma once

#include "GF.h"
#include "Qary_dist.h"

#include <cassert>
#include <cmath>


using namespace std;

typedef unsigned int bits32;
typedef unsigned short bits16;
//typedef vector<bool> BitVec;
//typedef vector<double> LLRVec;
#define ASSERT assert
#define REALMAX 1e20

typedef bool bit;
typedef double LLR;



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



class SC_Decoder_qary : public SC_Decoder
{
public:
	SC_Decoder_qary(int N, int m, const bit* frozen_bits, const GF& alpha);					// constructor of fully-frozen decoder.
	SC_Decoder_qary(int N, int m, const GF* frozen_syms, const GF& alpha);				// constructor of partially-frozen decoder.
	void sc_decode_qary(const qary_distribution* probs, bool is_Genie, const GF* true_u, bit* estimated_info_bits);

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


