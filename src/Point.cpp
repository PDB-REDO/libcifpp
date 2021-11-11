/*-
 * SPDX-License-Identifier: BSD-2-Clause
 * 
 * Copyright (c) 2020 NKI/AVL, Netherlands Cancer Institute
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <random>
#include <valarray>

#include "cif++/Point.hpp"

namespace mmcif
{

// --------------------------------------------------------------------
// uBlas compatible matrix types

// matrix is m x n, addressing i,j is 0 <= i < m and 0 <= j < n
// element m i,j is mapped to [i * n + j] and thus storage is row major

template <typename T>
class MatrixBase
{
  public:
	using value_type = T;

	virtual ~MatrixBase() {}

	virtual uint32_t dim_m() const = 0;
	virtual uint32_t dim_n() const = 0;

	virtual value_type &operator()(uint32_t i, uint32_t j) { throw std::runtime_error("unimplemented method"); }
	virtual value_type operator()(uint32_t i, uint32_t j) const = 0;

	MatrixBase &operator*=(const value_type &rhs);

	MatrixBase &operator-=(const value_type &rhs);
};

template <typename T>
MatrixBase<T> &MatrixBase<T>::operator*=(const T &rhs)
{
	for (uint32_t i = 0; i < dim_m(); ++i)
	{
		for (uint32_t j = 0; j < dim_n(); ++j)
		{
			operator()(i, j) *= rhs;
		}
	}

	return *this;
}

template <typename T>
MatrixBase<T> &MatrixBase<T>::operator-=(const T &rhs)
{
	for (uint32_t i = 0; i < dim_m(); ++i)
	{
		for (uint32_t j = 0; j < dim_n(); ++j)
		{
			operator()(i, j) -= rhs;
		}
	}

	return *this;
}

template <typename T>
class Matrix : public MatrixBase<T>
{
  public:
	using value_type = T;

	template <typename T2>
	Matrix(const MatrixBase<T2> &m)
		: m_m(m.dim_m())
		, m_n(m.dim_n())
	{
		m_data = new value_type[m_m * m_n];
		for (uint32_t i = 0; i < m_m; ++i)
		{
			for (uint32_t j = 0; j < m_n; ++j)
				operator()(i, j) = m(i, j);
		}
	}

	Matrix()
		: m_data(nullptr)
		, m_m(0)
		, m_n(0)
	{
	}

	Matrix(const Matrix &m)
		: m_m(m.m_m)
		, m_n(m.m_n)
	{
		m_data = new value_type[m_m * m_n];
		std::copy(m.m_data, m.m_data + (m_m * m_n), m_data);
	}

	Matrix &operator=(const Matrix &m)
	{
		value_type *t = new value_type[m.m_m * m.m_n];
		std::copy(m.m_data, m.m_data + (m.m_m * m.m_n), t);

		delete[] m_data;
		m_data = t;
		m_m = m.m_m;
		m_n = m.m_n;

		return *this;
	}

	Matrix(uint32_t m, uint32_t n, T v = T())
		: m_m(m)
		, m_n(n)
	{
		m_data = new value_type[m_m * m_n];
		std::fill(m_data, m_data + (m_m * m_n), v);
	}

	virtual ~Matrix()
	{
		delete[] m_data;
	}

	virtual uint32_t dim_m() const { return m_m; }
	virtual uint32_t dim_n() const { return m_n; }

	virtual value_type operator()(uint32_t i, uint32_t j) const
	{
		assert(i < m_m);
		assert(j < m_n);
		return m_data[i * m_n + j];
	}

	virtual value_type &operator()(uint32_t i, uint32_t j)
	{
		assert(i < m_m);
		assert(j < m_n);
		return m_data[i * m_n + j];
	}

	template <typename Func>
	void each(Func f)
	{
		for (uint32_t i = 0; i < m_m * m_n; ++i)
			f(m_data[i]);
	}

	template <typename U>
	Matrix &operator/=(U v)
	{
		for (uint32_t i = 0; i < m_m * m_n; ++i)
			m_data[i] /= v;

		return *this;
	}

  private:
	value_type *m_data;
	uint32_t m_m, m_n;
};

// --------------------------------------------------------------------

template <typename T>
class SymmetricMatrix : public MatrixBase<T>
{
  public:
	typedef typename MatrixBase<T>::value_type value_type;

	SymmetricMatrix(uint32_t n, T v = T())
		: m_owner(true)
		, m_n(n)
	{
		uint32_t N = (m_n * (m_n + 1)) / 2;
		m_data = new value_type[N];
		std::fill(m_data, m_data + N, v);
	}

	SymmetricMatrix(const T *data, uint32_t n)
		: m_owner(false)
		, m_data(const_cast<T *>(data))
		, m_n(n)
	{
	}

	virtual ~SymmetricMatrix()
	{
		if (m_owner)
			delete[] m_data;
	}

	virtual uint32_t dim_m() const { return m_n; }
	virtual uint32_t dim_n() const { return m_n; }

	T operator()(uint32_t i, uint32_t j) const;
	virtual T &operator()(uint32_t i, uint32_t j);

	// erase two rows, add one at the end (for neighbour joining)
	void erase_2(uint32_t i, uint32_t j);

	template <typename Func>
	void each(Func f)
	{
		uint32_t N = (m_n * (m_n + 1)) / 2;

		for (uint32_t i = 0; i < N; ++i)
			f(m_data[i]);
	}

	template <typename U>
	SymmetricMatrix &operator/=(U v)
	{
		uint32_t N = (m_n * (m_n + 1)) / 2;

		for (uint32_t i = 0; i < N; ++i)
			m_data[i] /= v;

		return *this;
	}

  private:
	bool m_owner;
	value_type *m_data;
	uint32_t m_n;
};

template <typename T>
inline T SymmetricMatrix<T>::operator()(uint32_t i, uint32_t j) const
{
	return i < j
	           ? m_data[(j * (j + 1)) / 2 + i]
	           : m_data[(i * (i + 1)) / 2 + j];
}

template <typename T>
inline T &SymmetricMatrix<T>::operator()(uint32_t i, uint32_t j)
{
	if (i > j)
		std::swap(i, j);
	assert(j < m_n);
	return m_data[(j * (j + 1)) / 2 + i];
}

template <typename T>
void SymmetricMatrix<T>::erase_2(uint32_t di, uint32_t dj)
{
	uint32_t s = 0, d = 0;
	for (uint32_t i = 0; i < m_n; ++i)
	{
		for (uint32_t j = 0; j < i; ++j)
		{
			if (i != di and j != dj and i != dj and j != di)
			{
				if (s != d)
					m_data[d] = m_data[s];
				++d;
			}

			++s;
		}
	}

	--m_n;
}

template <typename T>
class IdentityMatrix : public MatrixBase<T>
{
  public:
	typedef typename MatrixBase<T>::value_type value_type;

	IdentityMatrix(uint32_t n)
		: m_n(n)
	{
	}

	virtual uint32_t dim_m() const { return m_n; }
	virtual uint32_t dim_n() const { return m_n; }

	virtual value_type operator()(uint32_t i, uint32_t j) const
	{
		value_type result = 0;
		if (i == j)
			result = 1;
		return result;
	}

  private:
	uint32_t m_n;
};

// --------------------------------------------------------------------
// matrix functions

template <typename T>
Matrix<T> operator*(const MatrixBase<T> &lhs, const MatrixBase<T> &rhs)
{
	Matrix<T> result(std::min(lhs.dim_m(), rhs.dim_m()), std::min(lhs.dim_n(), rhs.dim_n()));

	for (uint32_t i = 0; i < result.dim_m(); ++i)
	{
		for (uint32_t j = 0; j < result.dim_n(); ++j)
		{
			for (uint32_t li = 0, rj = 0; li < lhs.dim_m() and rj < rhs.dim_n(); ++li, ++rj)
				result(i, j) += lhs(li, j) * rhs(i, rj);
		}
	}

	return result;
}

template <typename T>
Matrix<T> operator*(const MatrixBase<T> &lhs, T rhs)
{
	Matrix<T> result(lhs);
	result *= rhs;

	return result;
}

template <typename T>
Matrix<T> operator-(const MatrixBase<T> &lhs, const MatrixBase<T> &rhs)
{
	Matrix<T> result(std::min(lhs.dim_m(), rhs.dim_m()), std::min(lhs.dim_n(), rhs.dim_n()));

	for (uint32_t i = 0; i < result.dim_m(); ++i)
	{
		for (uint32_t j = 0; j < result.dim_n(); ++j)
		{
			result(i, j) = lhs(i, j) - rhs(i, j);
		}
	}

	return result;
}

template <typename T>
Matrix<T> operator-(const MatrixBase<T> &lhs, T rhs)
{
	Matrix<T> result(lhs.dim_m(), lhs.dim_n());
	result -= rhs;
	return result;
}

template<class M1, typename T>
void cofactors(const M1& m, SymmetricMatrix<T>& cf)
{
    const size_t ixs[4][3] =
    {
        { 1, 2, 3 },
        { 0, 2, 3 },
        { 0, 1, 3 },
        { 0, 1, 2 }
    };

    for (size_t x = 0; x < 4; ++x)
    {
        const size_t* ix = ixs[x];

        for (size_t y = x; y < 4; ++y)
        {
            const size_t* iy = ixs[y];

            cf(x, y) =
                m(ix[0], iy[0]) * m(ix[1], iy[1]) * m(ix[2], iy[2]) +
                m(ix[0], iy[1]) * m(ix[1], iy[2]) * m(ix[2], iy[0]) +
                m(ix[0], iy[2]) * m(ix[1], iy[0]) * m(ix[2], iy[1]) -
                m(ix[0], iy[2]) * m(ix[1], iy[1]) * m(ix[2], iy[0]) -
                m(ix[0], iy[1]) * m(ix[1], iy[0]) * m(ix[2], iy[2]) -
                m(ix[0], iy[0]) * m(ix[1], iy[2]) * m(ix[2], iy[1]);
			
			if ((x + y) % 2 == 1)
				cf(x, y) *= -1;
        }
    }
}

// --------------------------------------------------------------------

Quaternion Normalize(Quaternion q)
{
	std::valarray<double> t(4);
	
	t[0] = q.R_component_1();
	t[1] = q.R_component_2();
	t[2] = q.R_component_3();
	t[3] = q.R_component_4();
	
	t *= t;
	
	double length = std::sqrt(t.sum());

	if (length > 0.001)
		q /= static_cast<Quaternion::value_type>(length);
	else
		q = Quaternion(1, 0, 0, 0);

	return q;
}

// --------------------------------------------------------------------

std::tuple<double,Point> QuaternionToAngleAxis(Quaternion q)
{
	if (q.R_component_1() > 1)
		q = Normalize(q);

	// angle:
	double angle = 2 * acos(q.R_component_1());
	angle = angle * 180 / kPI;

	// axis:
	float s = std::sqrt(1 - q.R_component_1() * q.R_component_1());
	if (s < 0.001)
		s = 1;
	
	Point axis(q.R_component_2() / s, q.R_component_3() / s, q.R_component_4() / s);

	return std::make_tuple(angle, axis);
}

Point CenterPoints(std::vector<Point>& Points)
{
	Point t;
	
	for (Point& pt : Points)
	{
		t.mX += pt.mX;
		t.mY += pt.mY;
		t.mZ += pt.mZ;
	}
	
	t.mX /= Points.size();
	t.mY /= Points.size();
	t.mZ /= Points.size();
	
	for (Point& pt : Points)
	{
		pt.mX -= t.mX;
		pt.mY -= t.mY;
		pt.mZ -= t.mZ;
	}
	
	return t;
}

Point Centroid(std::vector<Point>& Points)
{
	Point result;
	
	for (Point& pt : Points)
		result += pt;
	
	result /= static_cast<float>(Points.size());
	
	return result;
}

double RMSd(const std::vector<Point>& a, const std::vector<Point>& b)
{
	double sum = 0;
	for (uint32_t i = 0; i < a.size(); ++i)
	{
		std::valarray<double> d(3);
		
		d[0] = b[i].mX - a[i].mX;
		d[1] = b[i].mY - a[i].mY;
		d[2] = b[i].mZ - a[i].mZ;

		d *= d;
		
		sum += d.sum();
	}
	
	return std::sqrt(sum / a.size());
}

// The next function returns the largest solution for a quartic equation
// based on Ferrari's algorithm.
// A depressed quartic is of the form:
//
//   x^4 + ax^2 + bx + c = 0
//
// (since I'm too lazy to find out a better way, I've implemented the
//  routine using complex values to avoid nan's as a result of taking
//  sqrt of a negative number)
double LargestDepressedQuarticSolution(double a, double b, double c)
{
	std::complex<double> P = - (a * a) / 12 - c;
	std::complex<double> Q = - (a * a * a) / 108 + (a * c) / 3 - (b * b) / 8;
	std::complex<double> R = - Q / 2.0 + std::sqrt((Q * Q) / 4.0 + (P * P * P) / 27.0);
	
	std::complex<double> U = std::pow(R, 1 / 3.0);
	
	std::complex<double> y;
	if (U == 0.0)
		y = -5.0 * a / 6.0 + U - std::pow(Q, 1.0 / 3.0);
	else
		y = -5.0 * a / 6.0 + U - P / (3.0 * U);

	std::complex<double> W = std::sqrt(a + 2.0 * y);
	
	// And to get the final result:
	// result = (±W + std::sqrt(-(3 * alpha + 2 * y ± 2 * beta / W))) / 2;
	// We want the largest result, so:

	std::valarray<double> t(4);

	t[0] = (( W + std::sqrt(-(3.0 * a + 2.0 * y + 2.0 * b / W))) / 2.0).real();
	t[1] = (( W + std::sqrt(-(3.0 * a + 2.0 * y - 2.0 * b / W))) / 2.0).real();
	t[2] = ((-W + std::sqrt(-(3.0 * a + 2.0 * y + 2.0 * b / W))) / 2.0).real();
	t[3] = ((-W + std::sqrt(-(3.0 * a + 2.0 * y - 2.0 * b / W))) / 2.0).real();

	return t.max();
}

Quaternion AlignPoints(const std::vector<Point>& pa, const std::vector<Point>& pb)
{
	// First calculate M, a 3x3 Matrix containing the sums of products of the coordinates of A and B
	Matrix<double> M(3, 3, 0);

	for (uint32_t i = 0; i < pa.size(); ++i)
	{
		const Point& a = pa[i];
		const Point& b = pb[i];
		
		M(0, 0) += a.mX * b.mX;	M(0, 1) += a.mX * b.mY;	M(0, 2) += a.mX * b.mZ;
		M(1, 0) += a.mY * b.mX;	M(1, 1) += a.mY * b.mY;	M(1, 2) += a.mY * b.mZ;
		M(2, 0) += a.mZ * b.mX;	M(2, 1) += a.mZ * b.mY;	M(2, 2) += a.mZ * b.mZ;
	}
	
	// Now calculate N, a symmetric 4x4 Matrix
	SymmetricMatrix<double> N(4);
	
	N(0, 0) =  M(0, 0) + M(1, 1) + M(2, 2);
	N(0, 1) =  M(1, 2) - M(2, 1);
	N(0, 2) =  M(2, 0) - M(0, 2);
	N(0, 3) =  M(0, 1) - M(1, 0);
	
	N(1, 1) =  M(0, 0) - M(1, 1) - M(2, 2);
	N(1, 2) =  M(0, 1) + M(1, 0);
	N(1, 3) =  M(0, 2) + M(2, 0);
	
	N(2, 2) = -M(0, 0) + M(1, 1) - M(2, 2);
	N(2, 3) =  M(1, 2) + M(2, 1);
	
	N(3, 3) = -M(0, 0) - M(1, 1) + M(2, 2);

	// det(N - λI) = 0
	// find the largest λ (λm)
	//
	// Aλ4 + Bλ3 + Cλ2 + Dλ + E = 0
	// A = 1
	// B = 0
	// and so this is a so-called depressed quartic
	// solve it using Ferrari's algorithm
	
	double C = -2 * (
		M(0, 0) * M(0, 0) + M(0, 1) * M(0, 1) + M(0, 2) * M(0, 2) +
		M(1, 0) * M(1, 0) + M(1, 1) * M(1, 1) + M(1, 2) * M(1, 2) +
		M(2, 0) * M(2, 0) + M(2, 1) * M(2, 1) + M(2, 2) * M(2, 2));
	
	double D = 8 * (M(0, 0) * M(1, 2) * M(2, 1) +
					M(1, 1) * M(2, 0) * M(0, 2) +
					M(2, 2) * M(0, 1) * M(1, 0)) -
			   8 * (M(0, 0) * M(1, 1) * M(2, 2) +
					M(1, 2) * M(2, 0) * M(0, 1) +
					M(2, 1) * M(1, 0) * M(0, 2));
	
	// E is the determinant of N:
	double E = 
		(N(0,0) * N(1,1) - N(0,1) * N(0,1)) * (N(2,2) * N(3,3) - N(2,3) * N(2,3)) +
		(N(0,1) * N(0,2) - N(0,0) * N(2,1)) * (N(2,1) * N(3,3) - N(2,3) * N(1,3)) +
		(N(0,0) * N(1,3) - N(0,1) * N(0,3)) * (N(2,1) * N(2,3) - N(2,2) * N(1,3)) +
		(N(0,1) * N(2,1) - N(1,1) * N(0,2)) * (N(0,2) * N(3,3) - N(2,3) * N(0,3)) +
		(N(1,1) * N(0,3) - N(0,1) * N(1,3)) * (N(0,2) * N(2,3) - N(2,2) * N(0,3)) +
		(N(0,2) * N(1,3) - N(2,1) * N(0,3)) * (N(0,2) * N(1,3) - N(2,1) * N(0,3));
	
	// solve quartic
	double lm = LargestDepressedQuarticSolution(C, D, E);
	
	// calculate t = (N - λI)
	Matrix<double> li = IdentityMatrix<double>(4) * lm;
	Matrix<double> t = N - li;
	
	// calculate a Matrix of cofactors for t, since N is symmetric, t must be symmetric as well and so will be cf
	SymmetricMatrix<double> cf(4);
	cofactors(t, cf);

	int maxR = 0;
	for (int r = 1; r < 4; ++r)
	{
		if (cf(r, 0) > cf(maxR, 0))
			maxR = r;
	}
	
	Quaternion q(cf(maxR, 0), cf(maxR, 1), cf(maxR, 2), cf(maxR, 3));
	q = Normalize(q);
	
	return q;
}

// --------------------------------------------------------------------

Point Nudge(Point p, float offset)
{
	static std::random_device rd;
	static std::mt19937_64 rng(rd());

	std::uniform_real_distribution<> randomAngle(0, 2 * kPI);
	std::normal_distribution<> randomOffset(0, offset);

	float theta = static_cast<float>(randomAngle(rng));
	float phi1 = static_cast<float>(randomAngle(rng) - kPI);
	float phi2 = static_cast<float>(randomAngle(rng) - kPI);
		
	Quaternion q = boost::math::spherical(1.0f, theta, phi1, phi2);

	Point r{ 0, 0, 1 };
	r.rotate(q);
	r *= static_cast<float>(randomOffset(rng));
	
	return p + r;
}

}
