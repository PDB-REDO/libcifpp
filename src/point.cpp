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

#include <cif++/point.hpp>

namespace mmcif
{

// --------------------------------------------------------------------
// We're using expression templates here

template <typename M>
class MatrixExpression
{
  public:
	uint32_t dim_m() const { return static_cast<const M&>(*this).dim_m(); }
	uint32_t dim_n() const { return static_cast<const M&>(*this).dim_n(); }

	double &operator()(uint32_t i, uint32_t j)
	{
		return static_cast<M&>(*this).operator()(i, j);
	}

	double operator()(uint32_t i, uint32_t j) const
	{
		return static_cast<const M&>(*this).operator()(i, j);
	}
};

// --------------------------------------------------------------------
// matrix is m x n, addressing i,j is 0 <= i < m and 0 <= j < n
// element m i,j is mapped to [i * n + j] and thus storage is row major

class Matrix : public MatrixExpression<Matrix>
{
  public:
	template <typename M2>
	Matrix(const MatrixExpression<M2> &m)
		: m_m(m.dim_m())
		, m_n(m.dim_n())
		, m_data(m_m * m_n)
	{
		for (uint32_t i = 0; i < m_m; ++i)
		{
			for (uint32_t j = 0; j < m_n; ++j)
				operator()(i, j) = m(i, j);
		}
	}

	Matrix(size_t m, size_t n, double v = 0)
		: m_m(m)
		, m_n(n)
		, m_data(m_m * m_n)
	{
		std::fill(m_data.begin(), m_data.end(), v);
	}

	Matrix() = default;
	Matrix(Matrix &&m) = default;
	Matrix(const Matrix &m) = default;
	Matrix &operator=(Matrix &&m) = default;
	Matrix &operator=(const Matrix &m) = default;

	uint32_t dim_m() const { return m_m; }
	uint32_t dim_n() const { return m_n; }

	double operator()(uint32_t i, uint32_t j) const
	{
		assert(i < m_m);
		assert(j < m_n);
		return m_data[i * m_n + j];
	}

	double &operator()(uint32_t i, uint32_t j)
	{
		assert(i < m_m);
		assert(j < m_n);
		return m_data[i * m_n + j];
	}

  private:
	uint32_t m_m = 0, m_n = 0;
	std::vector<double> m_data;
};

// --------------------------------------------------------------------

class SymmetricMatrix : public MatrixExpression<SymmetricMatrix>
{
  public:
	SymmetricMatrix(uint32_t n, double v = 0)
		: m_n(n)
		, m_data((m_n * (m_n + 1)) / 2)
	{
		std::fill(m_data.begin(), m_data.end(), v);
	}

	SymmetricMatrix() = default;
	SymmetricMatrix(SymmetricMatrix &&m) = default;
	SymmetricMatrix(const SymmetricMatrix &m) = default;
	SymmetricMatrix &operator=(SymmetricMatrix &&m) = default;
	SymmetricMatrix &operator=(const SymmetricMatrix &m) = default;

	uint32_t dim_m() const { return m_n; }
	uint32_t dim_n() const { return m_n; }

	double operator()(uint32_t i, uint32_t j) const
	{
		return i < j
				? m_data[(j * (j + 1)) / 2 + i]
				: m_data[(i * (i + 1)) / 2 + j];
	}

	double &operator()(uint32_t i, uint32_t j)
	{
		if (i > j)
			std::swap(i, j);
		assert(j < m_n);
		return m_data[(j * (j + 1)) / 2 + i];
	}

  private:
	uint32_t m_n;
	std::vector<double> m_data;
};

class IdentityMatrix : public MatrixExpression<IdentityMatrix>
{
  public:
	IdentityMatrix(uint32_t n)
		: m_n(n)
	{
	}

	uint32_t dim_m() const { return m_n; }
	uint32_t dim_n() const { return m_n; }

	double operator()(uint32_t i, uint32_t j) const
	{
		return i == j ? 1 : 0;
	}

  private:
	uint32_t m_n;
};

// --------------------------------------------------------------------
// matrix functions, implemented as expression templates

template<typename M1, typename M2>
class MatrixSubtraction : public MatrixExpression<MatrixSubtraction<M1, M2>>
{
  public:
	MatrixSubtraction(const M1 &m1, const M2 &m2)
		: m_m1(m1), m_m2(m2)
	{
		assert(m_m1.dim_m() == m_m2.dim_m());
		assert(m_m1.dim_n() == m_m2.dim_n());
	}

	uint32_t dim_m() const { return m_m1.dim_m(); }
	uint32_t dim_n() const { return m_m1.dim_n(); }

	double operator()(uint32_t i, uint32_t j) const
	{
		return m_m1(i, j) - m_m2(i, j);
	}

  private:
	const M1 &m_m1;
	const M2 &m_m2;
};

template<typename M1, typename M2>
MatrixSubtraction<M1, M2> operator-(const MatrixExpression<M1> &m1, const MatrixExpression<M2> &m2)
{
	return MatrixSubtraction(*static_cast<const M1*>(&m1), *static_cast<const M2*>(&m2));
}

template<typename M>
class MatrixMultiplication : public MatrixExpression<MatrixMultiplication<M>>
{
  public:
	MatrixMultiplication(const M &m, double v)
		: m_m(m), m_v(v)
	{
	}

	uint32_t dim_m() const { return m_m.dim_m(); }
	uint32_t dim_n() const { return m_m.dim_n(); }

	double operator()(uint32_t i, uint32_t j) const
	{
		return m_m(i, j) * m_v;
	}

  private:
	const M &m_m;
	double m_v;
};

template<typename M>
MatrixMultiplication<M> operator*(const MatrixExpression<M> &m, double v)
{
	return MatrixMultiplication(*static_cast<const M*>(&m), v);
}

// --------------------------------------------------------------------

template<class M1>
Matrix Cofactors(const M1& m)
{
	Matrix cf(m.dim_m(), m.dim_m());

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

        for (size_t y = 0; y < 4; ++y)
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

	return cf;
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

Quaternion ConstructFromAngleAxis(float angle, Point axis)
{
	auto q = std::cos((angle * mmcif::kPI / 180) / 2);
	auto s = std::sqrt(1 - q * q);

	axis.normalize();

	return Normalize(Quaternion{
		static_cast<float>(q),
		static_cast<float>(s * axis.mX),
		static_cast<float>(s * axis.mY),
		static_cast<float>(s * axis.mZ)});
}

std::tuple<double,Point> QuaternionToAngleAxis(Quaternion q)
{
	if (q.R_component_1() > 1)
		q = Normalize(q);

	// angle:
	double angle = 2 * std::acos(q.R_component_1());
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

Point Centroid(const std::vector<Point>& pts)
{
	Point result;
	
	for (auto &pt : pts)
		result += pt;
	
	result /= static_cast<float>(pts.size());
	
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
	Matrix M(3, 3, 0);

	for (uint32_t i = 0; i < pa.size(); ++i)
	{
		const Point& a = pa[i];
		const Point& b = pb[i];
		
		M(0, 0) += a.mX * b.mX;	M(0, 1) += a.mX * b.mY;	M(0, 2) += a.mX * b.mZ;
		M(1, 0) += a.mY * b.mX;	M(1, 1) += a.mY * b.mY;	M(1, 2) += a.mY * b.mZ;
		M(2, 0) += a.mZ * b.mX;	M(2, 1) += a.mZ * b.mY;	M(2, 2) += a.mZ * b.mZ;
	}
	
	// Now calculate N, a symmetric 4x4 Matrix
	SymmetricMatrix N(4);
	
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
	double lambda = LargestDepressedQuarticSolution(C, D, E);
	
	// calculate t = (N - λI)
	Matrix t = N - IdentityMatrix(4) * lambda;
	
	// calculate a Matrix of cofactors for t
	Matrix cf = Cofactors(t);

	int maxR = 0;
	for (int r = 1; r < 4; ++r)
	{
		if (std::abs(cf(r, 0)) > std::abs(cf(maxR, 0)))
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
