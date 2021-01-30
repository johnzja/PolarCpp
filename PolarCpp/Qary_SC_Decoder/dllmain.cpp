// dllmain.cpp : 定义 DLL 应用程序的入口点。
#include "pch.h"

int get_int(const mxArray* pm)
{
	if (mxIsDouble(pm))
		return (int)((mxGetDoubles(pm))[0]);
	else if (mxIsInt32(pm))
		return mxGetInt32s(pm)[0];
	else
		mexErrMsgTxt("Error! INT type expected.");
}

// MATLAB entry mexFunction.
/*
	MATLAB call for partially-frozen q-ary SC decoder:
	polar_info_esti = Qary_SC_Decoder(channel_recv, N, m, frozen_syms, alpha, decoder_config);
		channel_recv			LLR or probabilities for each component bit. LSB first.
		N						Code length in q-ary symbols. Must be power of two.
		m						q = 2^m.
		frozen_syms				q-ary array of length N, specifying which channels are frozen.
		alpha					Primitive element used in encoding.

		--------------------------------------------------------------------

	MATLAB call for completely-frozen q-ary SC decoder:
	polar_info_esti = Qary_SC_Decoder(channel_recv, N, m, frozen_bits, alpha, decoder_config);
		channel_recv			LLR or probabilities for each component symbol. LSB first.
		N						Code length in q-ary symbols. Must be power of two.
		m						q = 2^m.
		frozen_bits				bit array of length N, specifying which channels are frozen.
		alpha					Primitive element used in encoding.
		
		--------------------------------------------------------------------

		decoder_config.partially_frozen = { true, false }			default: false
		decoder_config.is_qary = { true, false }					default: true
		decoder_config.is_LLR = { true, false}						default: true. If it is false, then posteriori probabilities are input.
*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/*  Step1: Parameter format check. */
	if (nlhs != 1)
	{
		mexErrMsgTxt("Only 1 output vector expected.");
	}
	else if (nrhs != 6)
	{
		mexErrMsgTxt("Six input parameters expected");
	}

	// Fetch 6 parameters.
	const mxArray* mx_channel_recv		= prhs[0];
	const mxArray* mx_N					= prhs[1];
	const mxArray* mx_m					= prhs[2];
	const mxArray* mx_frozens			= prhs[3];
	const mxArray* mx_alpha				= prhs[4];
	const mxArray* mx_decoder_config	= prhs[5];

	// Validity check.
	/* Step2: Fetch the configuration boolean switches. */
	bool partially_frozen	= false;
	bool is_qary			= true;
	bool is_LLR				= true;
	
	
	if (!mxIsStruct(mx_decoder_config)) mexErrMsgTxt("decoder_config must be a struct.");
	mxArray *mx_partially_frozen, *mx_is_qary, *mx_is_LLR;
	mx_partially_frozen = mxGetField(mx_decoder_config, 0, "partially_frozen");		if (mx_partially_frozen && mxIsLogical(mx_partially_frozen)) partially_frozen = mxGetLogicals(mx_partially_frozen)[0];
	mx_is_qary			= mxGetField(mx_decoder_config, 0, "is_qary");				if (mx_is_qary && mxIsLogical(mx_is_qary)) is_qary = mxGetLogicals(mx_is_qary)[0];
	mx_is_LLR			= mxGetField(mx_decoder_config, 0, "is_LLR");				if (mx_is_LLR && mxIsLogical(mx_is_LLR)) is_LLR = mxGetLogicals(mx_is_LLR)[0];


	/* Step3: Fetch size parameters. */
	if (mxGetM(mx_N) != 1 || mxGetN(mx_N) != 1) mexErrMsgTxt("Input N must be an integer.");
	int N = get_int(mx_N);

	if (mxGetM(mx_m) != 1 || mxGetN(mx_m) != 1) mexErrMsgTxt("Input m must be an integer.");
	int m = get_int(mx_m);

	if (mxGetM(mx_alpha) != 1 || mxGetN(mx_alpha) != 1) mexErrMsgTxt("Input alpha must be an integer within range [0, (2^m)-1].");
	int alpha = get_int(mx_alpha);
	GF prim_element(m, alpha);


	/*  Step4: Construct the decoder according to configurations. */
	GF* frozen_syms		= NULL;
	bit* frozen_bits	= NULL;

	size_t N_rows = mxGetM(mx_frozens);
	if (N_rows != 1)mexErrMsgTxt("Input frozens must be a row vector.");
	size_t N_cols = mxGetN(mx_frozens);
	if (N_cols != N) mexErrMsgTxt("Input frozens must be of length N.");

	if (partially_frozen)
	{
		// must be double type.
		if (!mxIsDouble(mx_frozens)) mexErrMsgTxt("Input frozens must be of type double if partially-frozen.");
		double* frozen_arr = mxGetDoubles(mx_frozens);

		frozen_syms = new GF[N];
		for (int i = 0; i < N; i++)
		{
			frozen_syms[i] = GF(m, (short)(frozen_arr[i]));
		}
	}
	else
	{
		
		// must be logical or double.
		if (!mxIsLogical(mx_frozens))
		{
			if (!mxIsDouble(mx_frozens)) mexErrMsgTxt("Input frozens must be of type logical or double if completely-frozen.");
			// HERE: it is double.
			double* d_frozen_arr = mxGetDoubles(mx_frozens);
			frozen_bits = new bit[N];
			for (int i = 0; i < N; i++)
			{
				if (d_frozen_arr[i] == 0.0) frozen_bits[i] = false;
				else if (d_frozen_arr[i] == 1.0) frozen_bits[i] = true;
				else mexErrMsgTxt("Elements in frozen_bits must be binary.");
			}
		}
		else
		{
			// HERE: it is logical.
			bool* frozen_arr = mxGetLogicals(mx_frozens);

			frozen_bits = new bit[N];
			for (int i = 0; i < N; i++)
			{
				frozen_bits[i] = frozen_arr[i];
			}
		}
	}

	size_t row_vector[2]; row_vector[0] = 1;

	if (is_LLR)
	{
		size_t N_rows = mxGetM(mx_channel_recv);
		if (N_rows != 1)mexErrMsgTxt("Input LLR must be a row vector.");
		if (!mxIsDouble(mx_channel_recv)) mexErrMsgTxt("Input LLR must be doubles.");
		size_t N_cols = mxGetN(mx_channel_recv);
		if (is_qary)
		{
			// Q-ary SC. 
			ASSERT(N_cols == (N*m));
			SC_Decoder_qary* p_SC_Decoder_qary = NULL;

			// Construct q-ary distributions from input LLRs.
			qary_distribution* qdist = SC_Decoder_qary::convert_llr_into_qdist(N, m, mxGetDoubles(mx_channel_recv));

			if (partially_frozen)
			{
				p_SC_Decoder_qary = new SC_Decoder_qary(N, m, frozen_syms, prim_element);
				
			}
			else
			{
				p_SC_Decoder_qary = new SC_Decoder_qary(N, m, frozen_bits, prim_element);
			}

			row_vector[1] = p_SC_Decoder_qary->get_K();
			plhs[0] = mxCreateLogicalArray(2, row_vector);
			bit* result_ptr = mxGetLogicals(plhs[0]);
			p_SC_Decoder_qary->sc_decode_qary(qdist, result_ptr);		// Perform q-ary SC.

			delete p_SC_Decoder_qary;
			qary_distribution::destroyqd(qdist, N);						// ALERT: function convert_llr_into_qdist will use qary_distribution::newqd.
		}
		else
		{
			// Ordinary SC. Then LLR is expected to be of length N.
			ASSERT(N_cols == N);
			if (partially_frozen) mexErrMsgTxt("Binary code cannot be partially frozen.");

			SC_Decoder* p_SC_Decoder = new SC_Decoder(N, frozen_bits);		// default: Set binary = true.
			row_vector[1] = p_SC_Decoder->get_K();				// Get number of info bits.
			plhs[0] = mxCreateLogicalArray(2, row_vector);		// Create row vector of size [1, K].
			bit* result_ptr = mxGetLogicals(plhs[0]);

			double* LLRs = mxGetDoubles(mx_channel_recv);
			p_SC_Decoder->sc_decode(LLRs, result_ptr);			// Perform LLR-based binary SC.

			delete p_SC_Decoder;
		}
	}
	else
	{
		if (!is_qary) mexErrMsgTxt("Binary channel input must be in LLR form.");

		// Read out the distribution matrix directly from 
		size_t N_rows = mxGetM(mx_channel_recv);
		if (N_rows != (0x1 << m))mexErrMsgTxt("Input probabilities must be a matrix of 2^m rows.");
		if (!mxIsDouble(mx_channel_recv)) mexErrMsgTxt("Input probabilities must be doubles.");
		size_t N_cols = mxGetN(mx_channel_recv);
		if (N_cols != N) mexErrMsgTxt("Input probabilities must be a matrix of N columns.");
		double* channel_recv_probs = mxGetDoubles(mx_channel_recv);

		qary_distribution* qdist = qary_distribution::newqd(m, N);
		int q = (0x1) << m;
		for (int i = 0; i < N; i++)
		{
			qary_distribution& qd = qdist[i];
			for (int j = 0; j < q; j++)
			{
				qd.dist[j] = channel_recv_probs[i*q + j];
			}
		}

		// Perform q-ary SC.
		SC_Decoder_qary* p_SC_Decoder_qary = NULL;
		if (partially_frozen)
		{
			p_SC_Decoder_qary = new SC_Decoder_qary(N, m, frozen_syms, prim_element);

		}
		else
		{
			p_SC_Decoder_qary = new SC_Decoder_qary(N, m, frozen_bits, prim_element);
		}

		row_vector[1] = p_SC_Decoder_qary->get_K();
		plhs[0] = mxCreateLogicalArray(2, row_vector);
		bit* result_ptr = mxGetLogicals(plhs[0]);
		p_SC_Decoder_qary->sc_decode_qary(qdist, result_ptr);	// Perform q-ary SC.

		qary_distribution::destroyqd(qdist, N);
		delete p_SC_Decoder_qary;
	}


	delete[] frozen_syms;
	delete[] frozen_bits;
	//GF::destroy_GFTable();
	return;
}


BOOL APIENTRY DllMain( HMODULE hModule,
                       DWORD  ul_reason_for_call,
                       LPVOID lpReserved
                     )
{
    switch (ul_reason_for_call)
    {
    case DLL_PROCESS_ATTACH:
    case DLL_THREAD_ATTACH:
    case DLL_THREAD_DETACH:
    case DLL_PROCESS_DETACH:
        break;
    }
    return TRUE;
}
