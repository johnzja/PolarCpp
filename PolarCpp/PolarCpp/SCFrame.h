#pragma once

#include <cmath>
#include <cassert>

#include "GF.h"
#include "Qary_dist.h"

typedef bool bit;
typedef double LLR;


/* SC-Frame, providing P, C vectors, up/down calculations and left/right propagation functions.
* Properly using template and class inheritance may simplify coding.
*/
class SCFrame
{
public:
	SCFrame(int N, bool is_binary = true);

	SCFrame(const SCFrame& src);

	LLR left_propagate();

	void right_propagate(bit bit_decision);

	inline int get_current_index() const;										// return index{u_i}, 0 <= i <= N-1.

	inline void copy_from(const SCFrame& src);

	inline void setup_channel_recv(const LLR* channel_recv);

	inline void reset_all();

	virtual ~SCFrame();

protected:
	static void up_calculate(const LLR* llr_x1, const LLR* llr_x2, LLR* result, int len);
	static void down_calculate(const LLR* llr_x1, const LLR* llr_x2, const bit* u2, LLR* result, int len);

	int N, n;
	int phi;
	int* llr_layer_vec, * bit_layer_vec;
	bool copy, is_binary;				// HERE copy means to share llr_layer_vec and bit_layer_vec.
	LLR* P;
	bit* CL, * CR;
	const LLR* _channel_recv;
	LLR** P_src;		// perform lazy copying.
	bit** CL_src;		// perform lazy copying.
};

class Qary_SCFrame :public SCFrame
{
public:
	Qary_SCFrame(int N, int m, const GF& alpha);

	Qary_SCFrame(const Qary_SCFrame& src);

	const qary_distribution& left_propagate();

	void right_propagate(const GF& gf_decision);

	inline void copy_from(const Qary_SCFrame& src);

	inline void setup_channel_recv(const qary_distribution* channel_recv);

	inline void reset_all();

	virtual ~Qary_SCFrame();

protected:
	static void up_calculate(const qary_distribution* qdist_x1, const qary_distribution* qdist_x2, qary_distribution* result, const GF& alpha, int len);
	static void down_calculate(const qary_distribution* qdist_x1, const qary_distribution* qdist_x2, const GF* u1, qary_distribution* result, const GF& alpha, int len);

	int m;
	qary_distribution* P;
	GF* CL, * CR;
	const qary_distribution* _channel_recv;
	qary_distribution** P_src;
	GF** CL_src;
	GF alpha;
};
