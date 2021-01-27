#pragma once

#include <cassert>
#include <cmath>
#include "GF.h"
using namespace std;

typedef unsigned int bits32;
typedef unsigned short bits16;
//typedef vector<bool> BitVec;
//typedef vector<double> LLRVec;
#define ASSERT assert

typedef bool bit;
typedef double LLR;

/* Successive Cancellation Decoder for binary polar codes. */
class SC_Decoder
{
public:
	SC_Decoder(int N, const bit* frozen_bits);

	void sc_decode(const LLR* llr, bit* estimated_info_bits);

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



// q-ary SC decoder.
class SC_Decoder_qary;

class qary_distribution
{
	friend class SC_Decoder_qary;
public:
	qary_distribution(int m);
	~qary_distribution();

protected:
	double* dist;
	int m;
	int L;
};

class SC_Decoder_qary : protected SC_Decoder
{
public:
	SC_Decoder_qary(int N, int m, const bit* frozen_bits, const GF& alpha);
	void sc_decode_qary(const qary_distribution* probs, GF* estimated_info_syms);

protected:
	static void up_calculate(const qary_distribution* llr_x1, const qary_distribution* llr_x2, qary_distribution* result, GF alpha, int len);
	static void down_calculate(const qary_distribution* llr_x1, const qary_distribution* llr_x2, const GF* u1, qary_distribution* result, GF alpha, int len);

protected:
	GF alpha;
	int m;
};