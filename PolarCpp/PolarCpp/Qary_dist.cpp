#include "Qary_dist.h"

/*           q-ary probabilities                 
*	Posteriori probability distribution on GF(q).
*/
qary_distribution::qary_distribution(int m) :m(m)
{
	GF_ASSERT(m >= GF_M_MIN && m <= GF_M_MAX);
	L = (0x1) << m;
	dist = new double[L];
}

qary_distribution::~qary_distribution()
{
	delete[] dist;
}

qary_distribution& qary_distribution::operator=(const qary_distribution& src)
{
	if (this != &src)
	{
		GF_ASSERT(m = src.m);
		for (int i = 0; i < L; i++)
		{
			dist[i] = src.dist[i];
		}
	}

	return *this;
}

qary_distribution* qary_distribution::newqd(int m, int N)
{
	// generate qary_distribution objects.
	qary_distribution* y = (qary_distribution*) operator new(N * sizeof(qary_distribution));
	//qary_distribution* y = new qary_distribution[N];
	for (int i = 0; i < N; i++)
	{
		new (y + i)qary_distribution(m);	// initialize, using "placement new".
	}
	return y;
}

void qary_distribution::destroyqd(qary_distribution* pqd, int N)
{
	for (int i = 0; i < N; i++)
	{
		(pqd + i)->~qary_distribution();
	}
	operator delete ((void*)pqd);		// operator new.
}
