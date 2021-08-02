// mexmain.cpp : mexFunction is the entry point of a MEX.
#include "pch.h"
#include <assert.h>
#include <cmath>

const int q = 4;

void convert_dist_into_index(int idx[4], double probs[4], int N_bins);

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
	if (nlhs != 1)
	{
		mexErrMsgTxt("LHS must be a 3-D measure!");
	}
	if (nrhs != 4)
	{
		mexErrMsgTxt("RHS must have 4 parameters!");
	}

	const mxArray* mx_dist_1			= prhs[0];
	const mxArray* mx_dist_2			= prhs[1];
	const mxArray* mx_bin_centers		= prhs[2];
	const mxArray* mx_GF_info			= prhs[3];

	size_t N_bm_dims = mxGetNumberOfDimensions(mx_dist_1);
	if (N_bm_dims != 3) mexErrMsgTxt("Input blackwell measure dimensions must be 3.");

	// Get data size.
	assert(mxGetM(mx_bin_centers) == 1);
	int N_bins = mxGetN(mx_bin_centers);

	assert(mxIsDouble(mx_bin_centers));
	double* p_bin_centers = mxGetDoubles(mx_bin_centers);
	
	const mwSize* p_dims = mxGetDimensions(mx_dist_1);
	for (int i = 0; i < 3; i++)
	{
		if (p_dims[i] != N_bins) mexErrMsgTxt("Input blackwell measure size error!");
	}

	// Get the data.
	double* dist_1 = mxGetDoubles(mx_dist_1);
	double* dist_2 = mxGetDoubles(mx_dist_2);

	mxArray* mxAddTable, * mxMultTable, * mxAlpha;
	double* add_table = NULL, * mult_table = NULL, * p_alpha = NULL;
	int* i_add_table = NULL, * i_mult_table = NULL, alpha = 0;
	mxAddTable		= mxGetField(mx_GF_info, 0, "add_table");	if (mxAddTable && mxIsDouble(mxAddTable)) add_table = mxGetDoubles(mxAddTable);
	mxMultTable		= mxGetField(mx_GF_info, 0, "mult_table");	if (mxMultTable && mxIsDouble(mxMultTable)) mult_table = mxGetDoubles(mxMultTable);
	mxAlpha			= mxGetField(mx_GF_info, 0, "alpha");		if (mxAlpha && mxIsDouble(mxAlpha)) p_alpha = mxGetDoubles(mxAlpha);

	// GF add/mult tables.
	
	if (mxAddTable && mxMultTable && mxAlpha)
	{
		i_add_table = new int[q * q];
		i_mult_table = new int[q * q];
		for (int i = 0; i < q * q; i++)
		{
			i_add_table[i] = int(add_table[i]);
			i_mult_table[i] = int(mult_table[i]);
		}
		alpha = int(*p_alpha);
	}
	else mexErrMsgTxt("GF_info error.");

	int* proba_fetching_table_s = new int[q * q];
	int* proba_fetching_table_t = new int[q * q];
	for (int idx = 0; idx < 4; idx++)
	{
		for (int u2 = 0; u2 < 4; u2++)
		{
			int x1 = i_add_table[idx + q * (i_mult_table[alpha + q * u2])];
			int x2 = u2;
			proba_fetching_table_s[u2 + q * idx] = x1;
			proba_fetching_table_t[u2 + q * idx] = x2;
		}
	}

	// Construct the output BM.
	size_t size_vec[3];
	for (int i = 0; i < 3; i++) size_vec[i] = N_bins;
	plhs[0] = mxCreateNumericArray(3, size_vec, mxClassID::mxDOUBLE_CLASS, mxComplexity::mxREAL);
	double* dist_ret = mxGetDoubles(plhs[0]);

	// Fetch data: x[a,b,c] -> x[a+b*N_bins+c*N_bins^2].

	double s[4], t[4];
	int idx[4];

	for (int i0 = 0; i0 < N_bins; i0++)
	{
		s[0] = p_bin_centers[i0];
		for (int i1 = 0; i1 < (N_bins - i0); i1++)
		{
			s[1] = p_bin_centers[i1];
			for (int i2 = 0; i2 < (N_bins - i0 - i1); i2++)
			{
				double p1 = dist_1[i0 + N_bins * i1 + N_bins * N_bins * i2];
				s[2] = p_bin_centers[i2];
				s[3] = 1.0 - s[0] - s[1] - s[2];

				for (int j0 = 0; j0 < N_bins; j0++)
				{
					t[0] = p_bin_centers[j0];
					for (int j1 = 0; j1 < (N_bins - j0); j1++)
					{
						t[1] = p_bin_centers[j1];
						for (int j2 = 0; j2 < (N_bins - j0 - j1); j2++)
						{
							t[2] = p_bin_centers[j2];
							t[3] = 1.0 - t[0] - t[1] - t[2];
							double p2 = dist_2[j0 + N_bins * j1 + N_bins * N_bins * j2];

							double r[4] = { 0.0, 0.0, 0.0, 0.0 };
							for (int i = 0; i < 4; i++)
							{
								for (int u2 = 0; u2 < q; u2++)
								{
									r[i] += s[proba_fetching_table_s[u2 + q * i]] * t[proba_fetching_table_t[u2 + q * i]];
								}
							}
							//mexPrintf("%f\n", r[0] + r[1] + r[2] + r[3]);
							convert_dist_into_index(idx, r, N_bins);
							dist_ret[idx[0] + N_bins * idx[1] + N_bins * N_bins * idx[2]] += p1 * p2;
						}
					}
				}
				//mexPrintf("%f\n", dist_1[i2 + i1 * N_bins + i0 * N_bins * N_bins]);
			}
		}
	}
	
	delete[] proba_fetching_table_s;
	delete[] proba_fetching_table_t;
	delete[] i_add_table;
	delete[] i_mult_table;

	return;
}


void convert_dist_into_index(int idx[4], double probs[4], int N_bins)
{
	int s = 0;
	for (int i = 0; i < 4; i++)
	{
		s += (idx[i] = int(round((N_bins - 1) * probs[i])));
	}

	int d = s + 1 - N_bins;
	if (d > 0)
	{
		// Find the max index.
		int max_index = 0;
		int midx = idx[0];
		for (int i = 1; i < 4; i++)
		{
			if (midx < idx[i])
			{
				midx = idx[i];
				max_index = i;
			}
		}
		idx[max_index] -= d;
	}
	
	return;
}