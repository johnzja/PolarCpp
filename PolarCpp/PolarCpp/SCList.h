#pragma once

#include "SCFrame.h"

#include <algorithm>
#include <cassert>
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
