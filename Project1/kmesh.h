#pragma once
#ifndef _KMESH_
#define _KMESH_
using namespace std;
typedef vector<double> dvalarray;
class dvector3 :public dvalarray
{
public:
	dvector3(double x, double y, double z) :dvalarray(3)
	{
		(*this)[0] = x;
		(*this)[1] = y;
		(*this)[2] = z;
	}
	dvector3() :dvalarray(3) {};
	double operator *(const dvector3& a) const
	{
		return (*this)[0] * a[0] + (*this)[1] * a[1] + (*this)[2] * a[2];
	}
	double length() const
	{
		return sqrt(sqr((*this)[0]) + sqr((*this)[1]) + sqr((*this)[2]));
	}
	bool operator==(const dvector3& a)
	{
		return ((fabs((*this)[0] - a[0]) < 1e-10) && (fabs((*this)[1] - a[1]) < 1e-10) && (fabs((*this)[2] - a[2]) < 1e-10));
	}
	friend dvector3 operator-(const dvector3& a, const dvector3& b);
	friend dvector3 operator+(const dvector3& a, const dvector3& b);
	void Sum(const dvector3& a, const dvector3& b)
	{
		(*this)[0] = a[0] + b[0];
		(*this)[1] = a[1] + b[1];
		(*this)[2] = a[2] + b[2];
	}
};

inline dvector3 operator-(const dvector3& a, const dvector3& b)
{
	return dvector3(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

inline dvector3 operator+(const dvector3& a, const dvector3& b)
{
	return dvector3(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}

void ChooseIrreducible(list<dvector3>& kp, vector<dvector3>& irkp, vector<double>& wkp)
{
	long unsigned int kpsize = kp.size();
	list<dvector3> irkp0;
	list<int> wkp0;
	while (kp.size() > 0)
	{
		dvector3 tk(*kp.begin());
		irkp0.push_back(tk);
		wkp0.push_back(0);
		int& w0 = wkp0.back();
		for (int ix = 1; ix >= -1; ix -= 2)
		{
			for (int iy = 1; iy >= -1; iy -= 2)
			{
				for (int iz = 1; iz >= -1; iz -= 2)
				{
					dvector3 nk(ix * tk[0], iy * tk[1], iz * tk[2]);
					sort(nk.begin(), nk.end());
					do
					{
						list<dvector3>::iterator iv = find(kp.begin(), kp.end(), nk);
						if (iv != kp.end())
						{
							w0++;
							kp.erase(iv);
						}
					} while (next_permutation(nk.begin(), nk.end()));
				}
			}
		}
	}
	irkp.resize(irkp0.size());
	wkp.resize(wkp0.size());
	int j = 0;
	for (list<dvector3>::const_iterator ik = irkp0.begin(); ik != irkp0.end(); ik++)
	{
		irkp[j++] = *ik;
	}
	j = 0;
	for (list<int>::const_iterator iw = wkp0.begin(); iw != wkp0.end(); iw++)
	{
		wkp[j++] = (*iw) / static_cast<double>(kpsize);
	}
}

inline double kv0(int iq, int q)
{
	return (iq - static_cast<int>((q + 1.5) / 2) + 1) / static_cast<double>(q);
}

inline dvector3 operator *(const dvector3& a, double x)
{
	return dvector3(a[0] * x, a[1] * x, a[2] * x);
}

void GenerateAllKPoints(int q, const dvector3& b0, const dvector3& b1, const dvector3& b2, list<dvector3>& kp, int type = 0)
{
	if (type == 0)
	{
		for (int ip = 0; ip < q; ip++)
		{
			double p = kv0(ip, q);
			for (int ir = 0; ir < q; ir++)
			{
				double r = kv0(ir, q);
				for (int is = 0; is < q; is++)
				{
					double s = kv0(is, q);
					dvector3 k = b0 * p + b1 * r + b2 * s;
					kp.push_back(k);
				}
			}
		}
	}
	else
	{
		for (int ip = 0; ip < q; ip++)
		{
			double p = (2 * ip - q + 1.) / (2 * q);
			for (int ir = 0; ir < q; ir++)
			{
				double r = (2 * ir - q + 1.) / (2 * q);
				for (int is = 0; is < q; is++)
				{
					double s = (2 * is - q + 1.) / (2 * q);
					dvector3 k = b0 * p + b1 * r + b2 * s;
					kp.push_back(k);
				}
			}
		}
	}
}
inline dvector3 cross(const dvector3& a, const dvector3& b)
{
	return dvector3(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}
inline double VProduct(const dvector3& a, const dvector3& b, const dvector3& c)
{
	return cross(a, b) * c;
}
inline dvector3 operator *(double x, const dvector3& a)
{
	return dvector3(a[0] * x, a[1] * x, a[2] * x);
}
int cmp(const dvector3& a, const dvector3& b)
{
	return a.length() < b.length();
}
inline dvector3 operator /(const dvector3& a, double x)
{
	return dvector3(a[0] / x, a[1] / x, a[2] / x);
}
#endif // !_KMESH_