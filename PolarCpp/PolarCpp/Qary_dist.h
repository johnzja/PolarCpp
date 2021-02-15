#pragma once

#include "GF.h"

class qary_distribution
{
public:
	qary_distribution(int m);

	~qary_distribution();

	qary_distribution& operator=(const qary_distribution& src);

	static qary_distribution* newqd(int m, int N);

	static void destroyqd(qary_distribution* pqd, int N);

public:
	double* dist;
	int m;
	int L;
};
