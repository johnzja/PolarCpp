#include "SCList.h"


#define ASSERT assert
#define REALMAX 1e20
#define min(x,y) (((x)<(y))?(x):(y))

/* Binary SCL decoder. */
SCL_decoder::SCL_decoder(int N, const bit* frozen_bits, int list_size) :L(list_size), N(N)
{
	ASSERT((N & (N - 1)) == 0);
	ASSERT((L & (L - 1)) == 0);	// Ensure that L is a power of two.

	K = 0;
	_frozen_bits = new bit[N];
	for (int i = 0; i < N; i++)
	{
		if (!frozen_bits[i]) K++;
		_frozen_bits[i] = frozen_bits[i];
	}


	// Construct SC list.
	// using operator new.
	SCList = (SCFrame*) operator new(L * sizeof(SCFrame));
	new (SCList + 0) SCFrame(N);
	for (int i = 1; i < L; i++)
	{
		new (SCList + i) SCFrame(*(SCList));
	}

	PM = new double[L];
	is_active = new bool[L];
	u = new bit * [L];

	for (int i = 0; i < L; i++)
		u[i] = new bit[K];				// Create new arrays to store L decoding sequences.

	PM_0 = new double[L];
	PM_1 = new double[L];				// path splitting.

	ui_llr = new LLR[L];

	pwi = new PM_with_index[2 * L];
	active_0 = new bool[L];
	active_1 = new bool[L];

}

SCL_decoder::~SCL_decoder()
{
	delete[] _frozen_bits;
	for (int i = 0; i < L; i++)
	{
		(SCList + i)->~SCFrame();		// Call the destructor explicitly.
	}
	operator delete((void*)SCList);

	delete[] PM;
	delete[] is_active;
	for (int i = 0; i < L; i++) delete[] u[i];

	delete[] u;
	delete[] PM_0; delete[] PM_1;
	delete[] ui_llr;
	delete[] pwi;

	delete[] active_0;	delete[] active_1;
}


void SCL_decoder::scl_decode(const LLR* llr, bit* estimated_info_bits)
{
	// Step1: Load channel LLRs.
	for (int l_index = 0; l_index < L; l_index++)
	{
		SCList[l_index].reset_all();
		SCList[l_index].setup_channel_recv(llr);
		PM[l_index] = REALMAX;						// initialize path loss as MAX.
		is_active[l_index] = false;
	}
	is_active[0] = true;
	PM[0] = 0.0;

	int k = 0;

	// Step2: iteratively decode each bits.
	for (int phi = 0; phi < N; phi++)
	{
		while (!stk_killed.empty())stk_killed.pop();
		for (int i = 0; i < L; i++)
		{
			if (!is_active[i]) continue;
			ui_llr[i] = SCList[i].left_propagate();
		}

		if (_frozen_bits[phi])
		{
			for (int l_index = 0; l_index < L; l_index++)
			{
				if (!is_active[l_index]) continue;
				if (ui_llr[l_index] < 0) PM[l_index] -= ui_llr[l_index];
			}
		}
		else
		{
			// This is an info bit. Step1: Duplicate all paths.
			int num_active = 0;
			for (int l_index = 0; l_index < L; l_index++)
			{
				if (!is_active[l_index])
				{
					PM_0[l_index] = REALMAX;
					PM_1[l_index] = REALMAX;
				}
				else
				{
					LLR temp;
					if ((temp = ui_llr[l_index]) < 0)
					{
						PM_0[l_index] = PM[l_index] - temp;	// Add loss.
						PM_1[l_index] = PM[l_index];
					}
					else
					{
						PM_1[l_index] = PM[l_index] + temp;	// Add loss.
						PM_0[l_index] = PM[l_index];
					}
					num_active++;
				}
			}

			// Step2: Find at most L paths with smallest path measure.
			int num_paths_selected = min(L, 2 * num_active);
			for (int i = 0; i < L; i++)
			{
				pwi[i].x = PM_0[i];
				pwi[i].index = i;

				pwi[i + L].x = PM_1[i];
				pwi[i + L].index = i + L;
			}

			std::sort(pwi, pwi + 2 * L, [](const PM_with_index& p1, const PM_with_index& p2)->bool {return p1.x < p2.x; });		// using lambda expression to simplify the comparing code.

			// identify the killed paths.
			for (int i = 0; i < L; i++)
			{
				active_0[i] = false;
				active_1[i] = false;
			}

			int t = 0;
			while (num_paths_selected--)		// label the survived paths with 'true'.
			{
				PM_with_index& p = pwi[t++];
				if (p.index < L)
				{
					active_0[p.index] = true;
				}
				else
				{
					active_1[p.index - L] = true;
				}
			}

			//std::stack<int> stk_killed;

			for (int i = 0; i < L; i++)
			{
				if (!active_0[i] && !active_1[i])
				{
					stk_killed.push(i);
					is_active[i] = false;
				}
			}

			for (int i = 0; i < L; i++)
			{
				if (!is_active[i]) continue;

				if (active_0[i] && active_1[i])
				{
					// Perform path duplication when both 0 and 1 are promising paths.
					int new_path = stk_killed.top();	stk_killed.pop();

					is_active[new_path] = true;
					SCList[new_path].copy_from(SCList[i]);

					for (int j = 0; j < k; j++)
					{
						u[new_path][j] = u[i][j];
					}

					u[i][k] = false;			PM[i] = PM_0[i];
					u[new_path][k] = true;				PM[new_path] = PM_1[i];
				}
				else if (!active_0[i] && active_1[i])
				{
					u[i][k] = true; PM[i] = PM_1[i];
				}
				else if (active_0[i] && !active_1[i])
				{
					u[i][k] = false; PM[i] = PM_0[i];
				}
			}

			k++;
		}

		// Partial-sum return.
		for (int i = 0; i < L; i++)
		{
			if (!is_active[i]) continue;
			if (_frozen_bits[phi])
				SCList[i].right_propagate(0);
			else
				SCList[i].right_propagate(u[i][k - 1]);		// k++ executed when it is non-frozen.
		}

		// lazy_copy automatially configured.
	}

	// Step3: Find the path with smallest path measure.
	for (int i = 0; i < L; i++)
	{
		pwi[i].index = i;
		pwi[i].x = PM[i];
	}
	std::sort(pwi, pwi + L, [](const PM_with_index& p1, const PM_with_index& p2)->bool {return p1.x < p2.x; });		// using lambda expression to simplify the comparing code.

	bit* best_sequence = u[pwi[0].index];

	for (int j = 0; j < K; j++)
	{
		estimated_info_bits[j] = best_sequence[j];
	}
	return;
}



