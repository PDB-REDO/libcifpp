/*-
 * SPDX-License-Identifier: BSD-2-Clause
 *
 * Copyright (c) 2023 NKI/AVL, Netherlands Cancer Institute
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

#pragma once

#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <ostream>
#include <tuple>
#include <type_traits>
#include <vector>

namespace cif
{
// --------------------------------------------------------------------
// We're using expression templates here

template <typename M>
class matrix_expression
{
  public:
	constexpr uint32_t dim_m() const { return static_cast<const M &>(*this).dim_m(); }
	constexpr uint32_t dim_n() const { return static_cast<const M &>(*this).dim_n(); }

	constexpr auto &operator()(uint32_t i, uint32_t j)
	{
		return static_cast<M &>(*this).operator()(i, j);
	}

	constexpr auto operator()(uint32_t i, uint32_t j) const
	{
		return static_cast<const M &>(*this).operator()(i, j);
	}

	friend std::ostream &operator<<(std::ostream &os, const matrix_expression &m)
	{
		os << '[';

		for (size_t i = 0; i < m.dim_m(); ++i)
		{
			os << '[';

			for (size_t j = 0; j < m.dim_n(); ++j)
			{
				os << m(i, j);
				if (j + 1 < m.dim_n())
					os << ", ";
			}

			if (i + 1 < m.dim_m())
				os << ", ";

			os << ']';
		}

		os << ']';

		return os;
	}
};

// --------------------------------------------------------------------
// matrix is m x n, addressing i,j is 0 <= i < m and 0 <= j < n
// element m i,j is mapped to [i * n + j] and thus storage is row major

template <typename F = float>
class matrix : public matrix_expression<matrix<F>>
{
  public:
	using value_type = F;

	template <typename M2>
	matrix(const matrix_expression<M2> &m)
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

	matrix(size_t m, size_t n, value_type v = 0)
		: m_m(m)
		, m_n(n)
		, m_data(m_m * m_n)
	{
		std::fill(m_data.begin(), m_data.end(), v);
	}

	matrix() = default;
	matrix(matrix &&m) = default;
	matrix(const matrix &m) = default;
	matrix &operator=(matrix &&m) = default;
	matrix &operator=(const matrix &m) = default;

	constexpr size_t dim_m() const { return m_m; }
	constexpr size_t dim_n() const { return m_n; }

	constexpr value_type operator()(size_t i, size_t j) const
	{
		assert(i < m_m);
		assert(j < m_n);
		return m_data[i * m_n + j];
	}

	constexpr value_type &operator()(size_t i, size_t j)
	{
		assert(i < m_m);
		assert(j < m_n);
		return m_data[i * m_n + j];
	}

  private:
	size_t m_m = 0, m_n = 0;
	std::vector<value_type> m_data;
};

// --------------------------------------------------------------------
// special case, 3x3 matrix

template <typename F, size_t M, size_t N>
class matrix_fixed : public matrix_expression<matrix_fixed<F, M, N>>
{
  public:
	using value_type = F;

	template <typename M2>
	matrix_fixed(const M2 &m)
	{
		assert(M == m.dim_m() and N == m.dim_n());
		for (uint32_t i = 0; i < M; ++i)
		{
			for (uint32_t j = 0; j < N; ++j)
				operator()(i, j) = m(i, j);
		}
	}

	matrix_fixed(value_type v = 0)
	{
		std::fill(m_data.begin(), m_data.end(), v);
	}

	matrix_fixed(matrix_fixed &&m) = default;
	matrix_fixed(const matrix_fixed &m) = default;
	matrix_fixed &operator=(matrix_fixed &&m) = default;
	matrix_fixed &operator=(const matrix_fixed &m) = default;

	constexpr size_t dim_m() const { return M; }
	constexpr size_t dim_n() const { return N; }

	constexpr value_type operator()(size_t i, size_t j) const
	{
		assert(i < M);
		assert(j < N);
		return m_data[i * N + j];
	}

	constexpr value_type &operator()(size_t i, size_t j)
	{
		assert(i < M);
		assert(j < N);
		return m_data[i * N + j];
	}

  private:
	std::array<value_type, M * N> m_data;
};

template <typename F>
using matrix3x3 = matrix_fixed<F, 3, 3>;

template <typename F>
using matrix4x4 = matrix_fixed<F, 4, 4>;

// --------------------------------------------------------------------

template <typename F = float>
class symmetric_matrix : public matrix_expression<symmetric_matrix<F>>
{
  public:
	using value_type = F;

	symmetric_matrix(uint32_t n, value_type v = 0)
		: m_n(n)
		, m_data((m_n * (m_n + 1)) / 2)
	{
		std::fill(m_data.begin(), m_data.end(), v);
	}

	symmetric_matrix() = default;
	symmetric_matrix(symmetric_matrix &&m) = default;
	symmetric_matrix(const symmetric_matrix &m) = default;
	symmetric_matrix &operator=(symmetric_matrix &&m) = default;
	symmetric_matrix &operator=(const symmetric_matrix &m) = default;

	constexpr uint32_t dim_m() const { return m_n; }
	constexpr uint32_t dim_n() const { return m_n; }

	constexpr value_type operator()(uint32_t i, uint32_t j) const
	{
		return i < j
		           ? m_data[(j * (j + 1)) / 2 + i]
		           : m_data[(i * (i + 1)) / 2 + j];
	}

	constexpr value_type &operator()(uint32_t i, uint32_t j)
	{
		if (i > j)
			std::swap(i, j);
		assert(j < m_n);
		return m_data[(j * (j + 1)) / 2 + i];
	}

  private:
	uint32_t m_n;
	std::vector<value_type> m_data;
};

// --------------------------------------------------------------------

template <typename F, size_t M>
class symmetric_matrix_fixed : public matrix_expression<symmetric_matrix_fixed<F, M>>
{
  public:
	using value_type = F;

	symmetric_matrix_fixed(value_type v = 0)
	{
		std::fill(m_data.begin(), m_data.end(), v);
	}

	symmetric_matrix_fixed(symmetric_matrix_fixed &&m) = default;
	symmetric_matrix_fixed(const symmetric_matrix_fixed &m) = default;
	symmetric_matrix_fixed &operator=(symmetric_matrix_fixed &&m) = default;
	symmetric_matrix_fixed &operator=(const symmetric_matrix_fixed &m) = default;

	constexpr uint32_t dim_m() const { return M; }
	constexpr uint32_t dim_n() const { return M; }

	constexpr value_type operator()(uint32_t i, uint32_t j) const
	{
		return i < j
		           ? m_data[(j * (j + 1)) / 2 + i]
		           : m_data[(i * (i + 1)) / 2 + j];
	}

	constexpr value_type &operator()(uint32_t i, uint32_t j)
	{
		if (i > j)
			std::swap(i, j);
		assert(j < M);
		return m_data[(j * (j + 1)) / 2 + i];
	}

  private:
	std::array<value_type, (M * (M + 1)) / 2> m_data;
};

template <typename F>
using symmetric_matrix3x3 = symmetric_matrix_fixed<F, 3>;

template <typename F>
using symmetric_matrix4x4 = symmetric_matrix_fixed<F, 4>;

// --------------------------------------------------------------------

template <typename F = float>
class identity_matrix : public matrix_expression<identity_matrix<F>>
{
  public:
	using value_type = F;

	identity_matrix(uint32_t n)
		: m_n(n)
	{
	}

	constexpr uint32_t dim_m() const { return m_n; }
	constexpr uint32_t dim_n() const { return m_n; }

	constexpr value_type operator()(uint32_t i, uint32_t j) const
	{
		return i == j ? 1 : 0;
	}

  private:
	uint32_t m_n;
};

// --------------------------------------------------------------------
// matrix functions, implemented as expression templates

template <typename M1, typename M2>
class matrix_subtraction : public matrix_expression<matrix_subtraction<M1, M2>>
{
  public:
	matrix_subtraction(const M1 &m1, const M2 &m2)
		: m_m1(m1)
		, m_m2(m2)
	{
		assert(m_m1.dim_m() == m_m2.dim_m());
		assert(m_m1.dim_n() == m_m2.dim_n());
	}

	constexpr uint32_t dim_m() const { return m_m1.dim_m(); }
	constexpr uint32_t dim_n() const { return m_m1.dim_n(); }

	constexpr auto operator()(uint32_t i, uint32_t j) const
	{
		return m_m1(i, j) - m_m2(i, j);
	}

  private:
	const M1 &m_m1;
	const M2 &m_m2;
};

template <typename M1, typename M2>
auto operator-(const matrix_expression<M1> &m1, const matrix_expression<M2> &m2)
{
	return matrix_subtraction(m1, m2);
}

template <typename M1, typename M2>
class matrix_matrix_multiplication : public matrix_expression<matrix_matrix_multiplication<M1, M2>>
{
  public:
	matrix_matrix_multiplication(const M1 &m1, const M2 &m2)
		: m_m1(m1)
		, m_m2(m2)
	{
		assert(m1.dim_m() == m2.dim_n());
	}

	constexpr uint32_t dim_m() const { return m_m1.dim_m(); }
	constexpr uint32_t dim_n() const { return m_m1.dim_n(); }

	constexpr auto operator()(uint32_t i, uint32_t j) const
	{
		using value_type = decltype(m_m1(0, 0));

		value_type result = {};

		for (uint32_t k = 0; k < m_m1.dim_m(); ++k)
			result += m_m1(i, k) * m_m2(k, j);

		return result;
	}

  private:
	const M1 &m_m1;
	const M2 &m_m2;
};

template <typename M, typename T>
class matrix_scalar_multiplication : public matrix_expression<matrix_scalar_multiplication<M, T>>
{
  public:
	using value_type = T;

	matrix_scalar_multiplication(const M &m, value_type v)
		: m_m(m)
		, m_v(v)
	{
	}

	constexpr uint32_t dim_m() const { return m_m.dim_m(); }
	constexpr uint32_t dim_n() const { return m_m.dim_n(); }

	constexpr auto operator()(uint32_t i, uint32_t j) const
	{
		return m_m(i, j) * m_v;
	}

  private:
	const M &m_m;
	value_type m_v;
};

template <typename M1, typename T, std::enable_if_t<std::is_floating_point_v<T>, int> = 0>
auto operator*(const matrix_expression<M1> &m, T v)
{
	return matrix_scalar_multiplication(m, v);
}

template <typename M1, typename M2, std::enable_if_t<not std::is_floating_point_v<M2>, int> = 0>
auto operator*(const matrix_expression<M1> &m1, const matrix_expression<M2> &m2)
{
	return matrix_matrix_multiplication(m1, m2);
}

// --------------------------------------------------------------------

template <typename M>
auto determinant(const M &m);

template <typename F = float>
auto determinant(const matrix3x3<F> &m)
{
	return (m(0, 0) * (m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1)) +
			m(0, 1) * (m(1, 2) * m(2, 0) - m(1, 0) * m(2, 2)) +
			m(0, 2) * (m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0)));
}

template <typename M>
M inverse(const M &m);

template <typename F = float>
matrix3x3<F> inverse(const matrix3x3<F> &m)
{
	F det = determinant(m);

	matrix3x3<F> result;

	result(0, 0) = (m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1)) / det;
	result(1, 0) = (m(1, 2) * m(2, 0) - m(1, 0) * m(2, 2)) / det;
	result(2, 0) = (m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0)) / det;
	result(0, 1) = (m(2, 1) * m(0, 2) - m(2, 2) * m(0, 1)) / det;
	result(1, 1) = (m(2, 2) * m(0, 0) - m(2, 0) * m(0, 2)) / det;
	result(2, 1) = (m(2, 0) * m(0, 1) - m(2, 1) * m(0, 0)) / det;
	result(0, 2) = (m(0, 1) * m(1, 2) - m(0, 2) * m(1, 1)) / det;
	result(1, 2) = (m(0, 2) * m(1, 0) - m(0, 0) * m(1, 2)) / det;
	result(2, 2) = (m(0, 0) * m(1, 1) - m(0, 1) * m(1, 0)) / det;

	return result;
}

template <typename M>
auto eigen(const matrix_expression<M> &mat)
{
	using value_type = decltype(mat(0, 0));

	assert(mat.dim_m() == mat.dim_n());

	const size_t N = mat.dim_m();

	matrix<float> m = mat;

	matrix<value_type> em = identity_matrix(N);
	std::vector<value_type> ev(N);
	std::vector<value_type> b(N);

	// Set ev & b to diagonal.
	for (size_t i = 0; i < N; ++i)
		ev[i] = b[i] = m(i, i);

	for (int cyc = 1; cyc <= 50; ++cyc)
	{
		// calc sum of diagonal, off-diagonal
		value_type spp = 0, spq = 0;

		for (size_t i = 0; i < N - 1; i++)
		{
			for (size_t j = i + 1; j < N; j++)
				spq += std::fabs(m(i, j));
			spp += std::fabs(m(i, i));
		}

		if (spq <= 1.0e-12 * spp)
			break;

		std::vector<value_type> z(N);

		// now try and reduce each off-diagonal element in turn
		for (size_t i = 0; i < N - 1; i++)
		{
			for (size_t j = i + 1; j < N; j++)
			{
				value_type a_ij = m(i, j);
				value_type h = ev[j] - ev[i];
				value_type t;

				if (std::fabs(a_ij) > 1.0e-12 * std::fabs(h))
				{
					value_type theta = 0.5 * h / a_ij;
					t = 1.0 / (std::fabs(theta) + std::sqrt(1 + theta * theta));
					if (theta < 0)
						t = -t;
				}
				else
					t = a_ij / h;

				// calc trig properties
				value_type c = 1.f / std::sqrt(1 + t * t);
				value_type s = t * c;
				value_type tau = s / (1 + c);
				h = t * a_ij;

				// update eigenvalues
				z[i] -= h;
				z[j] += h;
				ev[i] -= h;
				ev[j] += h;

				// rotate the upper diagonal of the matrix
				m(i, j) = 0;

				for (size_t k = 0; k < i; k++)
				{
					value_type ai = m(k, i);
					value_type aj = m(k, j);
					m(k, i) = ai - s * (aj + ai * tau);
					m(k, j) = aj + s * (ai - aj * tau);
				}

				for (size_t k = i + 1; k < j; k++)
				{
					value_type ai = m(i, k);
					value_type aj = m(k, j);
					m(i, k) = ai - s * (aj + ai * tau);
					m(k, j) = aj + s * (ai - aj * tau);
				}

				for (size_t k = j + 1; k < N; k++)
				{
					value_type ai = m(i, k);
					value_type aj = m(j, k);
					m(i, k) = ai - s * (aj + ai * tau);
					m(j, k) = aj + s * (ai - aj * tau);
				}

				// apply corresponding rotation to result
				for (size_t k = 0; k < N; k++)
				{
					value_type ai = em(k, i);
					value_type aj = em(k, j);
					em(k, i) = ai - s * (aj + ai * tau);
					em(k, j) = aj + s * (ai - aj * tau);
				}
			}
		}

		for (size_t p = 0; p < N; p++)
		{
			b[p] += z[p];
			ev[p] = b[p];
		}
	}

	return std::make_tuple(ev, em);
}

// --------------------------------------------------------------------

template <typename M>
class matrix_cofactors : public matrix_expression<matrix_cofactors<M>>
{
  public:
	matrix_cofactors(const M &m)
		: m_m(m)
	{
	}

	constexpr uint32_t dim_m() const { return m_m.dim_m(); }
	constexpr uint32_t dim_n() const { return m_m.dim_n(); }

	constexpr auto operator()(uint32_t i, uint32_t j) const
	{
		const size_t ixs[4][3] = {
			{ 1, 2, 3 },
			{ 0, 2, 3 },
			{ 0, 1, 3 },
			{ 0, 1, 2 }
		};

		const size_t *ix = ixs[i];
		const size_t *iy = ixs[j];

		auto result =
			m_m(ix[0], iy[0]) * m_m(ix[1], iy[1]) * m_m(ix[2], iy[2]) +
			m_m(ix[0], iy[1]) * m_m(ix[1], iy[2]) * m_m(ix[2], iy[0]) +
			m_m(ix[0], iy[2]) * m_m(ix[1], iy[0]) * m_m(ix[2], iy[1]) -
			m_m(ix[0], iy[2]) * m_m(ix[1], iy[1]) * m_m(ix[2], iy[0]) -
			m_m(ix[0], iy[1]) * m_m(ix[1], iy[0]) * m_m(ix[2], iy[2]) -
			m_m(ix[0], iy[0]) * m_m(ix[1], iy[2]) * m_m(ix[2], iy[1]);

		return i + j % 2 == 1 ? -result : result;
	}

  private:
	const M &m_m;
};

} // namespace cif
