#include "GF.h"
#include <iostream>
using namespace std;

vector<GFTable*> GF::GFT_vec;

// Constructor.
GF::GF(int m, short x) : x(x), m(m)
{
	// Find the corresponding polynomial.
	GF_ASSERT(m >= GF_M_MIN && m <= GF_M_MAX);

	// deg(poly) = m.
	GF_ASSERT(0 <= x && x < (0x1 << m));
}

void GF::construct_GFTable(int m)
{
	if (GFT_vec.size() == 0)
	{
		for (int _m = 0; _m <= 10; _m++)
			GFT_vec.push_back(NULL);
	}
	if (GFT_vec[m] == NULL) GFT_vec[m] = new GFTable(m);
}

void GF::destroy_GFTable()
{
	for (GFTable* &p : GFT_vec)
	{
		if (p) delete p;
		p = NULL;
	}
}


GF operator+(const GF& a, const GF& b)
{
	GF_ASSERT(a.m == b.m && a.m != 0);
	return GF(a.m, a.x ^ b.x);
}

GF& operator+=(GF& a, const GF& b)
{
	a = a + b;
	return a;
}

static GF _mult(const GF& a, const GF& b)
{
	// Perform standard multiplication on GF(2^m).
	GF_ASSERT(a.m == b.m && a.m != 0);
	int m = a.m;
	short poly = GF::get_poly(m);
	short ans = 0;
	short ppa = a.x;			// partial poly-product of a.
	short ppb = b.x;

	for (int i = 0; i < m; i++)
	{
		if (ppb & 0x1)
		{
			ans ^= ppa;
		}
		ppa <<= 1;
		if (ppa & (0x1 << m))
		{
			ppa ^= poly;
		}
		ppb >>= 1;
	}
	return GF(m, ans);
}


GF operator*(const GF& a, const GF& b)
{
	GF_ASSERT(a.m == b.m);
	int m = a.m;
	GF::construct_GFTable(m);
	if (a.x == 0 || b.x == 0)
		return GF(m, 0);
	else if (a.x == 1)
		return b;
	else if (b.x == 1)
		return a;

	GFTable* gft = a.GFT_vec[m];				// Access GFTable.
	short log_a = gft->discrete_log[a.x];
	short log_b = gft->discrete_log[b.x];
	short log_ans = (log_a + log_b) % (gft->q - 1);
	return gft->primitive_chain[log_ans];
}

GF inv(const GF& a)
{
	GF_ASSERT(!a.is_zero());
	int m = a.m;
	if (a.is_one())return GF(m, 1);

	GF::construct_GFTable(m);
	GFTable* gft = a.GFT_vec[m];
	short log_a = gft->discrete_log[a.x];
	short log_b = (gft->q - 1) - log_a;
	return gft->primitive_chain[log_b];
}

GF operator/(const GF& a, const GF& b)
{
	GF_ASSERT(a.m == b.m);
	return a * inv(b);
}

ostream& operator<<(ostream& ostr, const GF& a)
{
	short t = a.x;
	for (int i = a.m-1; i >= 0; i--)
	{
		ostr << ((t & (0x1 << i)) != 0);
	}
	return ostr;
}

bool operator==(const GF& a, const GF& b)
{
	GF_ASSERT(a.m == b.m);
	return (a.x == b.x);
}

GFTable::GFTable(int m)
{
	GF_ASSERT(m >= 1 && m <= 10);
	poly = GF::get_poly(m);
	q = (0x1) << m;

	// Step1: Construct mult-table.
	GFMT = new GF[q*q];
	for (int i = 0; i < q; i++)
	{
		for (int j = i; j < q; j++)
		{
			GFMT[q*i + j] = _mult(GF(m, i), GF(m, j));		// stored by row.
			if (i != j) GFMT[q*j + i] = GFMT[q*i + j];
		}
	}

	// Step2: Find the primitive element.
	
	int i = 0;
	for (i = 2; i < q; i++)
	{
		// evaluate g(GFMT[i]).
		GF a(m, i);
		GF pow(m, 1);
		GF g_val(m, 0);
		short poly_temp = poly;
		for (int k = 0; k <= m; k++)
		{
			if (poly_temp & 0x1)
				g_val += pow;
			pow = _mult(pow, a);
			poly_temp >>= 1;
		}
		if (g_val.is_zero())
			break;
	}
	GF_ASSERT(i != q);

	primitive_chain = new GF[q - 1];
	primitive_chain[0] = GF(m, 1);
	primitive_chain[1] = GF(m, i);

	discrete_log = new short[q];	// discrete_log[x] means log(GF(x)) w.r.t primitive element alpha = primitive_chain[1].
	discrete_log[0] = -1;				// log(-1) is not defined.
	discrete_log[1] = 0;					

	for (int k = 2; k < q - 1; k++)
	{
		short t = primitive_chain[k-1].x;
		primitive_chain[k] = GFMT[t*q + i];
		discrete_log[primitive_chain[k].x] = k;
	}

	GF_ASSERT((_mult(primitive_chain[q - 2], primitive_chain[1])).is_one());
}