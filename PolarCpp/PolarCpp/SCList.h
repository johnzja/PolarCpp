#pragma once

#include "SCFrame.h"

#include <algorithm>
#include <cassert>
#include <intrin.h>
#include <stack>

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

	double* PM, * PM_0, * PM_1;
	bool* is_active, * active_0, * active_1;
	bit** u;
	LLR* ui_llr;
	PM_with_index* pwi;
	std::stack<int> stk_killed;
};

/* Not implemented yet because of difficulty in path duplication when GF(q) is considered. */
class Qary_SCL_decoder
{
public:
	Qary_SCL_decoder(int N, int m, const bit* frozen_bits, const GF& alpha, int list_size);

	Qary_SCL_decoder(int N, int m, const GF* frozen_syms, const GF& alpha, int list_size);

	virtual ~Qary_SCL_decoder();

	virtual void scl_decode(const qary_distribution* probs, bit* estimated_info_bits) = 0;

protected:
	int L, N, K;	
	bool partially_frozen;
	GF alpha;

	bit* _frozen_bits;
	GF* _frozen_syms;

	Qary_SCFrame* Q_SCList;

	double* PM, * PM_0, * PM_1;
	bool* is_active, * active_0, * active_1;
	GF** u;
	qary_distribution* ui_qdist;
	PM_with_index* pwi;
	std::stack<int> stk_killed;
};
