// mexmain.cpp : mexFunction is the entry point of a MEX.
#include "pch.h"
#include <assert.h>
#include <cmath>

#define fetch_element(arr, N, i0, i1, i2) (arr[(i0) + (N)*(i1) + (N)*(N)*(i2)])
const int q = 4;

void convert_dist_into_index(int idx[4], double probs[4], int N_bins);

/*
* bm_polar_transform(dist_1, dist_2, bin_centers, GF_info);
* 
*/
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
	if (nlhs != 1)
	{
		mexErrMsgTxt("LHS must be a 3-D measure!");
	}
	if (nrhs != 4 && nrhs != 5)
	{
		mexErrMsgTxt("RHS must have 4 or 5 parameters!");
	}

	const mxArray* mx_dist_1			= prhs[0];
	const mxArray* mx_dist_2			= prhs[1];
	const mxArray* mx_bin_centers		= prhs[2];
	const mxArray* mx_GF_info			= prhs[3];
	const mxArray* mx_is_down_transform = NULL;
	if (nrhs == 5) mx_is_down_transform = prhs[4];

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
	bool is_down_transform = false;

	mxAddTable		= mxGetField(mx_GF_info, 0, "add_table");	if (mxAddTable && mxIsDouble(mxAddTable)) add_table = mxGetDoubles(mxAddTable);
	mxMultTable		= mxGetField(mx_GF_info, 0, "mult_table");	if (mxMultTable && mxIsDouble(mxMultTable)) mult_table = mxGetDoubles(mxMultTable);
	mxAlpha			= mxGetField(mx_GF_info, 0, "alpha");		if (mxAlpha && mxIsDouble(mxAlpha)) p_alpha = mxGetDoubles(mxAlpha);

	if (mx_is_down_transform && mxIsLogical(mx_is_down_transform)) is_down_transform = (bool)(*(mxGetLogicals(mx_is_down_transform)));

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

	// Construct the output BM.
	size_t size_vec[3];
	for (int i = 0; i < 3; i++) size_vec[i] = N_bins;
	
	if (plhs[0]) mxFree(plhs[0]);
	plhs[0] = mxCreateNumericArray(3, size_vec, mxClassID::mxDOUBLE_CLASS, mxComplexity::mxREAL);
	double* dist_ret = mxGetDoubles(plhs[0]);

	// Fetch data: x[a,b,c] -> x[a+b*N_bins+c*N_bins^2].
	// initialize this dist.
	for (int i = 0; i < N_bins * N_bins * N_bins; i++)dist_ret[i] = 0.0;

	bool* dist_1_visited = new bool[N_bins * N_bins * N_bins];
	bool* dist_2_visited = new bool[N_bins * N_bins * N_bins];

	memset(dist_1_visited, 0, sizeof(bool) * N_bins * N_bins * N_bins);

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
				int i3 = N_bins - i0 - i1 - i2 - 1;

				bool b0 = fetch_element(dist_1_visited, N_bins, i0, i1, i2);
				bool b1 = fetch_element(dist_1_visited, N_bins, i1, i0, i3);
				bool b2 = fetch_element(dist_1_visited, N_bins, i2, i3, i0);
				bool b3 = fetch_element(dist_1_visited, N_bins, i3, i2, i1);

				if (b0 || b1 || b2 || b3) continue;

				double t0 = fetch_element(dist_1, N_bins, i0, i1, i2);
				double t1 = fetch_element(dist_1, N_bins, i1, i0, i3);
				double t2 = fetch_element(dist_1, N_bins, i2, i3, i0);
				double t3 = fetch_element(dist_1, N_bins, i3, i2, i1);

				fetch_element(dist_1_visited, N_bins, i0, i1, i2) = true;
				fetch_element(dist_1_visited, N_bins, i1, i0, i3) = true;
				fetch_element(dist_1_visited, N_bins, i2, i3, i0) = true;
				fetch_element(dist_1_visited, N_bins, i3, i2, i1) = true;

				double p1 = t0 + t1 + t2 + t3;
				if (p1 == 0.0) continue;

				s[2] = p_bin_centers[i2];
				s[3] = p_bin_centers[i3];

				memset(dist_2_visited, 0, sizeof(bool) * N_bins * N_bins * N_bins);

				for (int j0 = 0; j0 < N_bins; j0++)
				{
					t[0] = p_bin_centers[j0];
					for (int j1 = 0; j1 < (N_bins - j0); j1++)
					{
						t[1] = p_bin_centers[j1];
						for (int j2 = 0; j2 < (N_bins - j0 - j1); j2++)
						{
							int j3 = N_bins - j0 - j1 - j2 - 1;

							bool b0 = fetch_element(dist_2_visited, N_bins, j0, j1, j2);
							bool b1 = fetch_element(dist_2_visited, N_bins, j1, j0, j3);
							bool b2 = fetch_element(dist_2_visited, N_bins, j2, j3, j0);
							bool b3 = fetch_element(dist_2_visited, N_bins, j3, j2, j1);

							if (b0 || b1 || b2 || b3) continue;

							double t0 = fetch_element(dist_2, N_bins, j0, j1, j2);
							double t1 = fetch_element(dist_2, N_bins, j1, j0, j3);
							double t2 = fetch_element(dist_2, N_bins, j2, j3, j0);
							double t3 = fetch_element(dist_2, N_bins, j3, j2, j1);

							fetch_element(dist_2_visited, N_bins, j0, j1, j2) = true;
							fetch_element(dist_2_visited, N_bins, j1, j0, j3) = true;
							fetch_element(dist_2_visited, N_bins, j2, j3, j0) = true;
							fetch_element(dist_2_visited, N_bins, j3, j2, j1) = true;

							double p2 = t0 + t1 + t2 + t3;
							if (p2 == 0.0) continue;

							t[2] = p_bin_centers[j2];
							t[3] = p_bin_centers[j3];

							double r[4] = { 0.0, 0.0, 0.0, 0.0 };

							if (!is_down_transform)
							{
								// up-transform.
								r[0] = s[0] * t[0] + s[1] * t[3] + s[2] * t[1] + s[3] * t[2];
								r[1] = s[0] * t[3] + s[1] * t[0] + s[2] * t[2] + s[3] * t[1];
								r[2] = s[0] * t[1] + s[1] * t[2] + s[2] * t[0] + s[3] * t[3];
								r[3] = s[0] * t[2] + s[1] * t[1] + s[2] * t[3] + s[3] * t[0];

								//mexPrintf("%f\n", r[0] + r[1] + r[2] + r[3]);
								convert_dist_into_index(idx, r, N_bins);
								double ch_prob = p1 * p2 / 4;
								fetch_element(dist_ret, N_bins, idx[0], idx[1], idx[2]) += ch_prob;
								fetch_element(dist_ret, N_bins, idx[1], idx[0], idx[3]) += ch_prob;
								fetch_element(dist_ret, N_bins, idx[2], idx[3], idx[0]) += ch_prob;
								fetch_element(dist_ret, N_bins, idx[3], idx[2], idx[1]) += ch_prob;
							}
							else
							{
								// down-transform.
								double pc2[4] = { s[0], s[2], s[3], s[1] };
								double v[4][4];
								double q[4][4] = {
									{t[0], t[1], t[2], t[3]},
									{t[1], t[0], t[3], t[2]},
									{t[2], t[3], t[0], t[1]},
									{t[3], t[2], t[1], t[0]}
								};
								double lambdas[4] = { 0,0,0,0 };
								for (int i = 0; i < 4; i++)
									for (int j = 0; j < 4; j++)
										lambdas[i] += (v[i][j] = pc2[j] * q[i][j]);

								double ch_prob = p1 * p2 / 4;
								for (int i = 0; i < 4; i++)
								{
									if (lambdas[i] > 0)
									{
										for (int j = 0; j < 4; j++)r[j] = v[i][j] / lambdas[i];
										convert_dist_into_index(idx, r, N_bins);
										fetch_element(dist_ret, N_bins, idx[0], idx[1], idx[2]) += ch_prob * lambdas[i];
										fetch_element(dist_ret, N_bins, idx[1], idx[0], idx[3]) += ch_prob * lambdas[i];
										fetch_element(dist_ret, N_bins, idx[2], idx[3], idx[0]) += ch_prob * lambdas[i];
										fetch_element(dist_ret, N_bins, idx[3], idx[2], idx[1]) += ch_prob * lambdas[i];
									}
								}
							}
						}
					}
				}
			}
		}
	}
	
	delete[] i_add_table;
	delete[] i_mult_table;

	delete[] dist_1_visited;
	delete[] dist_2_visited;
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