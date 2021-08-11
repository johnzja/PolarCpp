// mexmain.cpp : mexFunction is the entry point of a MEX.
#include "pch.h"
#include <assert.h>
#include <cmath>

#include <semaphore.h>
#include <pthread.h>
#include <queue>
#include <mutex>

#pragma comment(lib, "pthreadVC2.lib")

#define fetch_element(arr, N, i0, i1, i2) (arr[(i0) + (N)*(i1) + (N)*(N)*(i2)])
const int q = 4;
std::queue<int> idle_threads;
sem_t sem;
std::mutex mtx;

void convert_dist_into_index(int idx[4], double probs[4], int N_bins);
void* worker_func(void* pthread_info);

typedef struct
{
	int thread_idx;
	volatile bool job_finished;
	bool need_to_join;
	bool is_down_transform;

	double* private_dist_ret;
	const double* dist1;
	const double* dist2;
	const double* p_bin_centers;
	const bool* KernelIndexVec0;
	const int* KernelIndexVec1;
	const int* KernelIndexMat0;
	const int* KernelIndexMat1;

	int N_jobs;
	int* pJobIndex;
	int N_bins;
}thread_info;

// Use multi-threading.
// 
void aggregate_thread_result(int N_bins, double* dist_ret, pthread_t* thread_pool, volatile thread_info* thread_infos, int index_thread)
{
	pthread_join(thread_pool[index_thread], NULL);	// This thread has exited!! Wait until it exits.
	// Aggregate the data.
	double* thread_ret = thread_infos[index_thread].private_dist_ret;
	for (int idx0 = 0; idx0 < N_bins; idx0++)
		for (int idx1 = 0; idx1 < N_bins - idx0; idx1++)
			for (int idx2 = 0; idx2 < N_bins - idx0 - idx1; idx2++)
			{
				fetch_element(dist_ret, N_bins, idx0, idx1, idx2) += fetch_element(thread_ret, N_bins, idx0, idx1, idx2);
			}
}

void dispatch_job(int N_bins, int N_threads, int N_jobsEachThread, int N_jobsEachThreadMax, double* dist_ret, pthread_t* thread_pool, pthread_attr_t* thread_pool_attr, volatile thread_info* thread_infos, int* pJobIndex)
{
	// Wait until one worker thread is ready.
	sem_wait(&sem);

	// Find an idle thread.
	mtx.lock();
	volatile int i = idle_threads.front();
	idle_threads.pop();
	mtx.unlock();

	// If index i is in the queue "idle_threads" and this thread has been dispatched, then job_finished must be true.
	if (thread_infos[i].job_finished)
	{
		aggregate_thread_result(N_bins, dist_ret, thread_pool, thread_infos, i);
	}

	if (thread_infos[i].private_dist_ret == NULL)
		thread_infos[i].private_dist_ret = new double[N_bins * N_bins * N_bins];
	for (int k = 0; k < N_bins * N_bins * N_bins; k++) thread_infos[i].private_dist_ret[k] = 0.0;

	thread_infos[i].N_jobs = N_jobsEachThread;
	thread_infos[i].job_finished = false;

	if (thread_infos[i].pJobIndex == NULL)
		thread_infos[i].pJobIndex = new int[N_jobsEachThreadMax];
	memcpy_s(thread_infos[i].pJobIndex, N_jobsEachThread * sizeof(int),
		pJobIndex, N_jobsEachThread * sizeof(int));

	thread_infos[i].need_to_join = true;
	pthread_create(&thread_pool[i], &thread_pool_attr[i], worker_func, (void*)(&thread_infos[i]));
}

void* worker_func(void* pthread_info)
{
	thread_info* ti = (thread_info*)pthread_info;
	// do the job.

	int N_bins = ti->N_bins;
	const bool* KernelIndexVec0 = ti->KernelIndexVec0;
	const int* KernelIndexVec1 = ti->KernelIndexVec1;
	const int* KernelIndexMat0 = ti->KernelIndexMat0;
	const int* KernelIndexMat1 = ti->KernelIndexMat1;
	const double* p_bin_centers = ti->p_bin_centers;
	const double* dist_1 = ti->dist1;
	const double* dist_2 = ti->dist2;
	double* dist_ret = ti->private_dist_ret;

	bool is_down_transform = ti->is_down_transform;
	double s[4], t[4];
	int idx[4];

	// Process range.
	for (int iter = 0; iter < ti->N_jobs; iter++)
	{
		int i0 = (ti->pJobIndex)[iter];
		int lb_i1 = KernelIndexVec1[i0];
		int ub_i1 = KernelIndexVec1[i0 + N_bins];

		s[0] = p_bin_centers[i0];
		for (int i1 = lb_i1 - 1; i1 < ub_i1; i1++)
		{
			s[1] = p_bin_centers[i1];
			int lb_i2 = KernelIndexMat0[i0 + N_bins * i1];
			int ub_i2 = KernelIndexMat1[i0 + N_bins * i1];

			for (int i2 = lb_i2 - 1; i2 < ub_i2; i2++)
			{
				int i3 = N_bins - i0 - i1 - i2 - 1;

				double t0 = fetch_element(dist_1, N_bins, i0, i1, i2);
				double t1 = fetch_element(dist_1, N_bins, i1, i0, i3);
				double t2 = fetch_element(dist_1, N_bins, i2, i3, i0);
				double t3 = fetch_element(dist_1, N_bins, i3, i2, i1);

				double p1 = t0 + t1 + t2 + t3;
				if (p1 == 0.0) continue;

				s[2] = p_bin_centers[i2];
				s[3] = p_bin_centers[i3];

				for (int j0 = 0; j0 < N_bins; j0++)
				{
					if (!KernelIndexVec0[j0]) continue;
					int lb_j1 = KernelIndexVec1[j0];
					int ub_j1 = KernelIndexVec1[j0 + N_bins];

					t[0] = p_bin_centers[j0];
					for (int j1 = lb_j1 - 1; j1 < ub_j1; j1++)
					{
						t[1] = p_bin_centers[j1];
						int lb_j2 = KernelIndexMat0[j0 + N_bins * j1];
						int ub_j2 = KernelIndexMat1[j0 + N_bins * j1];

						for (int j2 = lb_j2 - 1; j2 < ub_j2; j2++)
						{
							int j3 = N_bins - j0 - j1 - j2 - 1;

							double t0 = fetch_element(dist_2, N_bins, j0, j1, j2);
							double t1 = fetch_element(dist_2, N_bins, j1, j0, j3);
							double t2 = fetch_element(dist_2, N_bins, j2, j3, j0);
							double t3 = fetch_element(dist_2, N_bins, j3, j2, j1);

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

	// The job is done HERE.
	ti->job_finished = true;

	mtx.lock();
	idle_threads.push(ti->thread_idx);
	mtx.unlock();
	
	sem_post(&sem);
	return 0;
}

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

	mxArray* mxAddTable, * mxMultTable, * mxAlpha,
		* mxKernelIndexVec0, * mxKernelIndexVec1,
		* mxKernelIndexMat0, * mxKernelIndexMat1,
		* mxNumThreads, * mxJobsEachThread;
	double* add_table = NULL, * mult_table = NULL, * p_alpha = NULL;
	int* i_add_table = NULL, * i_mult_table = NULL, alpha = 0;
	bool is_down_transform = false;
	bool* KernelIndexVec0 = NULL;
	int* KernelIndexVec1 = NULL;
	int* KernelIndexMat0 = NULL;
	int* KernelIndexMat1 = NULL;
	int N_threads = 6;
	int N_jobsEachThread = 4;

	mxAddTable			= mxGetField(mx_GF_info, 0, "add_table");			if (mxAddTable && mxIsDouble(mxAddTable)) add_table = mxGetDoubles(mxAddTable);
	mxMultTable			= mxGetField(mx_GF_info, 0, "mult_table");			if (mxMultTable && mxIsDouble(mxMultTable)) mult_table = mxGetDoubles(mxMultTable);
	mxAlpha				= mxGetField(mx_GF_info, 0, "alpha");				if (mxAlpha && mxIsDouble(mxAlpha)) p_alpha = mxGetDoubles(mxAlpha);
	mxKernelIndexVec0	= mxGetField(mx_GF_info, 0, "kernel_index_vec0");	if (mxKernelIndexVec0 && mxIsLogical(mxKernelIndexVec0)) KernelIndexVec0 = mxGetLogicals(mxKernelIndexVec0);
	mxKernelIndexVec1	= mxGetField(mx_GF_info, 0, "kernel_index_vec1");	if (mxKernelIndexVec1 && mxIsInt32(mxKernelIndexVec1)) KernelIndexVec1 = mxGetInt32s(mxKernelIndexVec1);
	mxKernelIndexMat0	= mxGetField(mx_GF_info, 0, "kernel_index_mat0");	if (mxKernelIndexMat0 && mxIsInt32(mxKernelIndexMat0)) KernelIndexMat0 = mxGetInt32s(mxKernelIndexMat0);
	mxKernelIndexMat1	= mxGetField(mx_GF_info, 0, "kernel_index_mat1");	if (mxKernelIndexMat1 && mxIsInt32(mxKernelIndexMat1)) KernelIndexMat1 = mxGetInt32s(mxKernelIndexMat1);

	mxNumThreads = mxGetField(mx_GF_info, 0, "num_threads");
	if (mxNumThreads)
	{
		if (mxIsDouble(mxNumThreads))
		{
			N_threads = int(*(mxGetDoubles(mxNumThreads)));
		}
		else if (mxIsInt32(mxNumThreads))
		{
			N_threads = *(mxGetInt32s(mxNumThreads));
		}
	}
	else mexWarnMsgTxt("GF_info.num_threads set to default value 6.");

	mxJobsEachThread = mxGetField(mx_GF_info, 0, "jobs_each_thread");
	if (mxJobsEachThread)
	{
		if (mxIsDouble(mxJobsEachThread))
		{
			N_jobsEachThread = int(*(mxGetDoubles(mxJobsEachThread)));
		}
		else if (mxIsInt32(mxNumThreads))
		{
			N_jobsEachThread = *(mxGetInt32s(mxJobsEachThread));
		}
	}
	else mexWarnMsgTxt("GF_info.jobs_each_thread set to default value 4.");

	if (mx_is_down_transform && mxIsLogical(mx_is_down_transform)) is_down_transform = (bool)(*(mxGetLogicals(mx_is_down_transform)));
	if (!KernelIndexVec0 || !KernelIndexVec1 || !KernelIndexMat0 || !KernelIndexMat1)
		mexErrMsgTxt("GF_info must contain kernel indices!");

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


	int* pJobIndex						= new int[N_jobsEachThread];	
	pthread_t* thread_pool				= new pthread_t [N_threads];
	pthread_attr_t* thread_pool_attr	= new pthread_attr_t[N_threads];
	thread_info* thread_infos			= new thread_info[N_threads];


	sem_init(&sem, 0, N_threads);	// initial value is N_threads.
	for (int i = 0; i < N_threads; i++)
	{
		thread_infos[i].thread_idx = i;
		thread_infos[i].job_finished = false;
		thread_infos[i].need_to_join = false;
		thread_infos[i].is_down_transform = is_down_transform;

		thread_infos[i].N_bins = N_bins;
		thread_infos[i].dist1 = dist_1;
		thread_infos[i].dist2 = dist_2;
		thread_infos[i].p_bin_centers = p_bin_centers;
		thread_infos[i].KernelIndexMat0 = KernelIndexMat0;
		thread_infos[i].KernelIndexMat1 = KernelIndexMat1;
		thread_infos[i].KernelIndexVec0 = KernelIndexVec0;
		thread_infos[i].KernelIndexVec1 = KernelIndexVec1;
		thread_infos[i].private_dist_ret = NULL;
		thread_infos[i].pJobIndex = NULL;

		pthread_attr_init(thread_pool_attr + i);
		pthread_attr_setscope(thread_pool_attr + i, PTHREAD_SCOPE_PROCESS);
		pthread_attr_setdetachstate(thread_pool_attr + i, PTHREAD_CREATE_JOINABLE);

		idle_threads.push(i);	// Threads are not yet dispatched, so no mutex is needed.
	}
	
	int idxJob = 0;
	for (int i0 = 0; i0 < N_bins; i0++)
	{
		// dispatch at most N_threads.
		if (KernelIndexVec0[i0])
		{
			pJobIndex[idxJob++] = i0;
		}
		if (idxJob == N_jobsEachThread)
		{
			idxJob = 0;
			dispatch_job(N_bins, N_threads, N_jobsEachThread, N_jobsEachThread, dist_ret, thread_pool, thread_pool_attr, thread_infos, pJobIndex);
		}
	}
	if (idxJob > 0)
	{
		dispatch_job(N_bins, N_threads, idxJob, N_jobsEachThread, dist_ret, thread_pool, thread_pool_attr, thread_infos, pJobIndex);
	}
	
	// Wait for all the threads.
	for (int i = 0; i < N_threads; i++)
	{
		if (thread_infos[i].need_to_join)
		{
			aggregate_thread_result(N_bins, dist_ret, thread_pool, thread_infos, i);
		}
		delete[] thread_infos[i].private_dist_ret;
		delete[] thread_infos[i].pJobIndex;
	}

	delete[] pJobIndex;
	delete[] thread_pool;
	delete[] thread_pool_attr;
	delete[] thread_infos;
	sem_destroy(&sem);

	delete[] i_add_table;
	delete[] i_mult_table;
	while (!idle_threads.empty())idle_threads.pop();
	return;
}


void convert_dist_into_index(int idx[4], double probs[4], int N_bins)
{
	int s = 0;
	for (int i = 0; i < 4; i++)
	{
		s += (idx[i] = int(round((N_bins - 1) * probs[i])));
	}

	int d = s - (N_bins - 1);
	if (d == 0)
	{
		return;
	}
	else if (d == 1)
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
		idx[max_index]--;
	}
	else if (d == -1)
	{
		int min_index = 0;
		int midx = idx[0];
		for (int i = 1; i < 4; i++)
		{
			if (midx > idx[i])
			{
				midx = idx[i];
				min_index = i;
			}
		}
		idx[min_index]++;
	}
	else mexErrMsgTxt("convert_dist_into_index: unknown error.");

	return;
}
