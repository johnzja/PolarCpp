#pragma once
#include <cassert>
#include <vector>
using namespace std;

#define GF_ASSERT assert
#define GF_M_MAX (10)
#define GF_M_MIN (1)

class GFTable;
class GF;

// GF(2^m) field:  C++ Implementation.
class GF
{
public:
	GF()
	{
		m = 0;	// m=0 means GF is uninitialized.
	}

	GF(int m, short x);
	bool is_one() const
	{
		return x == 1;
	}

	bool is_zero() const
	{
		return x == 0;
	}

	static short get_poly(int m)
	{
		int poly = 0;
		switch (m)
		{
		case 1:
			poly = 3;
			break;
		case 2:
			poly = 7;
			break;
		case 3:
			poly = 11;
			break;
		case 4:
			poly = 19;
			break;
		case 5:
			poly = 37;
			break;
		case 6:
			poly = 67;
			break;
		case 7:
			poly = 137;
			break;
		case 8:
			poly = 285;
			break;
		case 9:
			poly = 529;
			break;
		case 10:
			poly = 1033;
			break;
		}
		return poly;
	}

	static void construct_GFTable(int m);

	static void destroy_GFTable();

	GF get_prim();

	int m;
	short x;
	static vector<GFTable*> GFT_vec;
};

class GFTable
{
public:
	GFTable(int m);

	~GFTable()
	{
		delete[] GFMT;
		delete[] primitive_chain;
		delete[] discrete_log;
	}

	//private:
	int m;
	short poly;
	int q;

	GF* GFMT;	// GF multiplication table.
	GF* primitive_chain;	// vector of length q-1 that starts with the primitive element.
	short* discrete_log;
};

/* GF functions */
GF operator+(const GF& a, const GF& b);
GF& operator+=(GF& a, const GF& b);

GF operator*(const GF& a, const GF& b);
GF inv(const GF& a);
GF operator/(const GF& a, const GF& b);
ostream& operator<<(ostream& ostr, const GF& a);
bool operator==(const GF& a, const GF& b);
