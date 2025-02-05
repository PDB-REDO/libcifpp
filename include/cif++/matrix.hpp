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

/**
 * @file matrix.hpp
 * 
 * Some basic matrix operations and classes to hold matrices.
 * 
 * We're using expression templates for optimal performance.
 * 
 */

namespace cif
{
// --------------------------------------------------------------------
// We're using expression templates here

/**
 * @brief Base for the matrix expression templates
 * This all uses the Curiously recurring template pattern
 * 
 * @tparam M The type of the derived class
 */
template <typename M>
class matrix_expression
{
  public:
	constexpr std::size_t dim_m() const { return static_cast<const M &>(*this).dim_m(); } ///< Return the size (dimension) in direction m
	constexpr std::size_t dim_n() const { return static_cast<const M &>(*this).dim_n(); } ///< Return the size (dimension) in direction n

	constexpr bool empty() const { return dim_m() == 0 or dim_n() == 0; } ///< Convenient way to test for empty matrices

	/** Return a reference to element [ @a i, @a j ] */
	constexpr auto &operator()(std::size_t i, std::size_t j)
	{
		return static_cast<M &>(*this).operator()(i, j);
	}

	/** Return the value of element [ @a i, @a j ] */
	constexpr auto operator()(std::size_t i, std::size_t j) const
	{
		return static_cast<const M &>(*this).operator()(i, j);
	}

	/** Swap the contents of rows @a r1 and @a r2 */
	void swap_row(std::size_t r1, std::size_t r2)
	{
		for (std::size_t c = 0; c < dim_m(); ++c)
		{
			auto v = operator()(r1, c);
			operator()(r1, c) = operator()(r2, c);
			operator()(r2, c) = v;
		}
	}

	/** Swap the contents of columns @a c1 and @a c2 */
	void swap_col(std::size_t c1, std::size_t c2)
	{
		for (std::size_t r = 0; r < dim_n(); ++r)
		{
			auto &a = operator()(r, c1);
			auto &b = operator()(r, c2);
			std::swap(a, b);
		}
	}

	/** write the matrix @a m to std::ostream @a os */
	friend std::ostream &operator<<(std::ostream &os, const matrix_expression &m)
	{
		os << '[';

		for (std::size_t i = 0; i < m.dim_m(); ++i)
		{
			os << '[';

			for (std::size_t j = 0; j < m.dim_n(); ++j)
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

/**
 * @brief Storage class implementation of matrix_expression.
 * 
 * @tparam F The type of the stored values
 *  
 * matrix is m x n, addressing i,j is 0 <= i < m and 0 <= j < n
 * element m i,j is mapped to [i * n + j] and thus storage is row major
 */

template <typename F = float>
class matrix : public matrix_expression<matrix<F>>
{
  public:
	/** The value type */
	using value_type = F;

	/**
	 * @brief Copy construct a new matrix object using @a m
	 * 
	 * @tparam M2 Type of @a m
	 * @param m The matrix expression to copy values from
	 */
	template <typename M2>
	matrix(const matrix_expression<M2> &m)
		: m_m(m.dim_m())
		, m_n(m.dim_n())
		, m_data(m_m * m_n)
	{
		for (std::size_t i = 0; i < m_m; ++i)
		{
			for (std::size_t j = 0; j < m_n; ++j)
				operator()(i, j) = m(i, j);
		}
	}

	/**
	 * @brief Construct a new matrix object with dimension @a m and @a n
	 * setting the values to @a v
	 * 
	 * @param m Requested dimension M
	 * @param n Requested dimension N
	 * @param v Value to store in each element
	 */
	matrix(std::size_t m, std::size_t n, value_type v = 0)
		: m_m(m)
		, m_n(n)
		, m_data(m_m * m_n)
	{
		std::fill(m_data.begin(), m_data.end(), v);
	}

	/** @cond */
	matrix() = default;
	matrix(matrix &&m) = default;
	matrix(const matrix &m) = default;
	matrix &operator=(matrix &&m) = default;
	matrix &operator=(const matrix &m) = default;
	/** @endcond */

	constexpr std::size_t dim_m() const { return m_m; } ///< Return dimension m
	constexpr std::size_t dim_n() const { return m_n; } ///< Return dimension n

	/** Return the value of element [ @a i, @a j ] */
	constexpr value_type operator()(std::size_t i, std::size_t j) const
	{
		assert(i < m_m);
		assert(j < m_n);
		return m_data[i * m_n + j];
	}

	/** Return a reference to element [ @a i, @a j ] */
	constexpr value_type &operator()(std::size_t i, std::size_t j)
	{
		assert(i < m_m);
		assert(j < m_n);
		return m_data[i * m_n + j];
	}

  private:
	std::size_t m_m = 0, m_n = 0;
	std::vector<value_type> m_data;
};

// --------------------------------------------------------------------
// special case, 3x3 matrix

/**
 * @brief Storage class implementation of matrix_expression
 * with compile time fixed size.
 * 
 * @tparam F The type of the stored values
 *  
 * matrix is m x n, addressing i,j is 0 <= i < m and 0 <= j < n
 * element m i,j is mapped to [i * n + j] and thus storage is row major
 */

template <typename F, std::size_t M, std::size_t N>
class matrix_fixed : public matrix_expression<matrix_fixed<F, M, N>>
{
  public:
	/** The value type */
	using value_type = F;

	/** The storage size */
	static constexpr std::size_t kSize = M * N;

	/** Copy constructor */
	template <typename M2>
	matrix_fixed(const M2 &m)
	{
		assert(M == m.dim_m() and N == m.dim_n());
		for (std::size_t i = 0; i < M; ++i)
		{
			for (std::size_t j = 0; j < N; ++j)
				operator()(i, j) = m(i, j);
		}
	}

	/** default constructor */
	matrix_fixed(value_type v = 0)
	{
		m_data.fill(v);
	}

	/** Alternate constructor taking an array of values to store */
	matrix_fixed(const F (&v)[kSize])
	{
		fill(v, std::make_index_sequence<kSize>{});
	}

	/** @cond */
	matrix_fixed(matrix_fixed &&m) = default;
	matrix_fixed(const matrix_fixed &m) = default;
	matrix_fixed &operator=(matrix_fixed &&m) = default;
	matrix_fixed &operator=(const matrix_fixed &m) = default;
	/** @endcond */

	/** Store the values in @a a in the matrix */
	template<std::size_t... Ixs>
	matrix_fixed& fill(const F (&a)[kSize], std::index_sequence<Ixs...>)
	{
		m_data = { a[Ixs]... };
		return *this;
	}

	constexpr std::size_t dim_m() const { return M; } ///< Return dimension m
	constexpr std::size_t dim_n() const { return N; } ///< Return dimension n

	/** Return the value of element [ @a i, @a j ] */
	constexpr value_type operator()(std::size_t i, std::size_t j) const
	{
		assert(i < M);
		assert(j < N);
		return m_data[i * N + j];
	}

	/** Return a reference to element [ @a i, @a j ] */
	constexpr value_type &operator()(std::size_t i, std::size_t j)
	{
		assert(i < M);
		assert(j < N);
		return m_data[i * N + j];
	}

  private:
	std::array<value_type, M * N> m_data;
};

/** typedef of a fixed matrix of size 3x3 */
template <typename F>
using matrix3x3 = matrix_fixed<F, 3, 3>;

/** typedef of a fixed matrix of size 4x4 */
template <typename F>
using matrix4x4 = matrix_fixed<F, 4, 4>;

// --------------------------------------------------------------------

/**
 * @brief Storage class implementation of symmetric matrix_expression
 * 
 * @tparam F The type of the stored values
 *  
 * matrix is m x n, addressing i,j is 0 <= i < m and 0 <= j < n
 * element m i,j is mapped to [i * n + j] and thus storage is row major
 */
template <typename F = float>
class symmetric_matrix : public matrix_expression<symmetric_matrix<F>>
{
  public:
	/** The value type */
	using value_type = F;

	/** constructor for a matrix of size @a n x @a n elements with value @a v */
	symmetric_matrix(std::size_t n, value_type v = 0)
		: m_n(n)
		, m_data((m_n * (m_n + 1)) / 2)
	{
		std::fill(m_data.begin(), m_data.end(), v);
	}

	/** @cond */
	symmetric_matrix() = default;
	symmetric_matrix(symmetric_matrix &&m) = default;
	symmetric_matrix(const symmetric_matrix &m) = default;
	symmetric_matrix &operator=(symmetric_matrix &&m) = default;
	symmetric_matrix &operator=(const symmetric_matrix &m) = default;
	/** @endcond */

	constexpr std::size_t dim_m() const { return m_n; } ///< Return dimension m
	constexpr std::size_t dim_n() const { return m_n; } ///< Return dimension n

	/** Return the value of element [ @a i, @a j ] */
	constexpr value_type operator()(std::size_t i, std::size_t j) const
	{
		return i < j
		           ? m_data[(j * (j + 1)) / 2 + i]
		           : m_data[(i * (i + 1)) / 2 + j];
	}

	/** Return a reference to element [ @a i, @a j ] */
	constexpr value_type &operator()(std::size_t i, std::size_t j)
	{
		if (i > j)
			std::swap(i, j);
		assert(j < m_n);
		return m_data[(j * (j + 1)) / 2 + i];
	}

  private:
	std::size_t m_n;
	std::vector<value_type> m_data;
};

// --------------------------------------------------------------------

/**
 * @brief Storage class implementation of symmetric matrix_expression
 * with compile time fixed size.
 * 
 * @tparam F The type of the stored values
 *  
 * matrix is m x n, addressing i,j is 0 <= i < m and 0 <= j < n
 * element m i,j is mapped to [i * n + j] and thus storage is row major
 */
template <typename F, std::size_t M>
class symmetric_matrix_fixed : public matrix_expression<symmetric_matrix_fixed<F, M>>
{
  public:
	/** The value type */
	using value_type = F;

	/** constructor with all elements set to value @a v */
	symmetric_matrix_fixed(value_type v = 0)
	{
		std::fill(m_data.begin(), m_data.end(), v);
	}

	/** @cond */
	symmetric_matrix_fixed(symmetric_matrix_fixed &&m) = default;
	symmetric_matrix_fixed(const symmetric_matrix_fixed &m) = default;
	symmetric_matrix_fixed &operator=(symmetric_matrix_fixed &&m) = default;
	symmetric_matrix_fixed &operator=(const symmetric_matrix_fixed &m) = default;
	/** @endcond */

	constexpr std::size_t dim_m() const { return M; } ///< Return dimension m
	constexpr std::size_t dim_n() const { return M; } ///< Return dimension n

	/** Return the value of element [ @a i, @a j ] */
	constexpr value_type operator()(std::size_t i, std::size_t j) const
	{
		return i < j
		           ? m_data[(j * (j + 1)) / 2 + i]
		           : m_data[(i * (i + 1)) / 2 + j];
	}

	/** Return a reference to element [ @a i, @a j ] */
	constexpr value_type &operator()(std::size_t i, std::size_t j)
	{
		if (i > j)
			std::swap(i, j);
		assert(j < M);
		return m_data[(j * (j + 1)) / 2 + i];
	}

  private:
	std::array<value_type, (M * (M + 1)) / 2> m_data;
};

/** typedef of a fixed symmetric matrix of size 3x3 */
template <typename F>
using symmetric_matrix3x3 = symmetric_matrix_fixed<F, 3>;

/** typedef of a fixed symmetric matrix of size 4x4 */
template <typename F>
using symmetric_matrix4x4 = symmetric_matrix_fixed<F, 4>;

// --------------------------------------------------------------------

/**
 * @brief implementation of symmetric matrix_expression with a value
 * of 1 for the diagonal values and 0 for all the others.
 *  
 * @tparam F The type of the stored values
 *  
 * matrix is m x n, addressing i,j is 0 <= i < m and 0 <= j < n
 * element m i,j is mapped to [i * n + j] and thus storage is row major
 */
template <typename F = float>
class identity_matrix : public matrix_expression<identity_matrix<F>>
{
  public:
	/** the value type */
	using value_type = F;

	/** constructor taking a dimension @a n */
	identity_matrix(std::size_t n)
		: m_n(n)
	{
	}

	constexpr std::size_t dim_m() const { return m_n; } ///< Return dimension m
	constexpr std::size_t dim_n() const { return m_n; } ///< Return dimension n

	/** Return the value of element [ @a i, @a j ] */
	constexpr value_type operator()(std::size_t i, std::size_t j) const
	{
		return static_cast<value_type>(i == j ? 1 : 0);
	}

  private:
	std::size_t m_n;
};

// --------------------------------------------------------------------
// matrix functions, implemented as expression templates

/**
 * @brief Implementation of a substraction operation as a matrix expression
 * 
 * @tparam M1 Type of matrix 1
 * @tparam M2 Type of matrix 2
 */
template <typename M1, typename M2>
class matrix_subtraction : public matrix_expression<matrix_subtraction<M1, M2>>
{
  public:
	/** constructor */
	matrix_subtraction(const M1 &m1, const M2 &m2)
		: m_m1(m1)
		, m_m2(m2)
	{
		assert(m_m1.dim_m() == m_m2.dim_m());
		assert(m_m1.dim_n() == m_m2.dim_n());
	}

	constexpr std::size_t dim_m() const { return m_m1.dim_m(); } ///< Return dimension m
	constexpr std::size_t dim_n() const { return m_m1.dim_n(); } ///< Return dimension n

	/** Access to the value of element [ @a i, @a j ] */
	constexpr auto operator()(std::size_t i, std::size_t j) const
	{
		return m_m1(i, j) - m_m2(i, j);
	}

  private:
	const M1 &m_m1;
	const M2 &m_m2;
};

/** operator to subtract two matrices and return a matrix expression */
template <typename M1, typename M2>
auto operator-(const matrix_expression<M1> &m1, const matrix_expression<M2> &m2)
{
	return matrix_subtraction(m1, m2);
}

/**
 * @brief Implementation of a multiplication operation as a matrix expression
 * 
 * @tparam M1 Type of matrix 1
 * @tparam M2 Type of matrix 2
 */
template <typename M1, typename M2>
class matrix_matrix_multiplication : public matrix_expression<matrix_matrix_multiplication<M1, M2>>
{
  public:
	/** constructor */
	matrix_matrix_multiplication(const M1 &m1, const M2 &m2)
		: m_m1(m1)
		, m_m2(m2)
	{
		assert(m1.dim_m() == m2.dim_n());
	}

	constexpr std::size_t dim_m() const { return m_m1.dim_m(); } ///< Return dimension m
	constexpr std::size_t dim_n() const { return m_m1.dim_n(); } ///< Return dimension n

	/** Access to the value of element [ @a i, @a j ] */
	constexpr auto operator()(std::size_t i, std::size_t j) const
	{
		using value_type = decltype(m_m1(0, 0));

		value_type result = {};

		for (std::size_t k = 0; k < m_m1.dim_m(); ++k)
			result += m_m1(i, k) * m_m2(k, j);

		return result;
	}

  private:
	const M1 &m_m1;
	const M2 &m_m2;
};

/**
 * @brief Implementation of a multiplication operation of a matrix and a scalar value as a matrix expression
 * 
 * @tparam M1 Type of matrix
 * @tparam M2 Type of scalar value
 */
template <typename M, typename T>
class matrix_scalar_multiplication : public matrix_expression<matrix_scalar_multiplication<M, T>>
{
  public:
	/** value type */
	using value_type = T;

	/** constructor */
	matrix_scalar_multiplication(const M &m, value_type v)
		: m_m(m)
		, m_v(v)
	{
	}

	constexpr std::size_t dim_m() const { return m_m.dim_m(); } ///< Return dimension m
	constexpr std::size_t dim_n() const { return m_m.dim_n(); } ///< Return dimension n

	/** Access to the value of element [ @a i, @a j ] */
	constexpr auto operator()(std::size_t i, std::size_t j) const
	{
		return m_m(i, j) * m_v;
	}

  private:
	const M &m_m;
	value_type m_v;
};

/** First implementation of operator*, enabled if the second parameter is a scalar */
template <typename M1, typename T, std::enable_if_t<std::is_floating_point_v<T>, int> = 0>
auto operator*(const matrix_expression<M1> &m, T v)
{
	return matrix_scalar_multiplication(m, v);
}

/** First implementation of operator*, enabled if the second parameter is not a scalar and thus must be a matrix, right? */
template <typename M1, typename M2, std::enable_if_t<not std::is_floating_point_v<M2>, int> = 0>
auto operator*(const matrix_expression<M1> &m1, const matrix_expression<M2> &m2)
{
	return matrix_matrix_multiplication(m1, m2);
}

// --------------------------------------------------------------------

/** Generic routine to calculate the determinant of a matrix
 * 
 * @note This is currently only implemented for fixed matrices of size 3x3
 */
template <typename M>
auto determinant(const M &m);

/** Implementation of the determinant function for fixed size matrices of size 3x3 */
template <typename F = float>
auto determinant(const matrix3x3<F> &m)
{
	return (m(0, 0) * (m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1)) +
			m(0, 1) * (m(1, 2) * m(2, 0) - m(1, 0) * m(2, 2)) +
			m(0, 2) * (m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0)));
}

/** Generic routine to calculate the inverse of a matrix
 * 
 * @note This is currently only implemented for fixed matrices of size 3x3
 */
template <typename M>
M inverse(const M &m);

/** Implementation of the inverse function for fixed size matrices of size 3x3 */
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

// --------------------------------------------------------------------

/**
 * @brief Implementation of a cofactor calculation as a matrix expression
 * 
 * @tparam M Type of matrix
 */
template <typename M>
class matrix_cofactors : public matrix_expression<matrix_cofactors<M>>
{
  public:
	/** constructor */
	matrix_cofactors(const M &m)
		: m_m(m)
	{
	}

	constexpr std::size_t dim_m() const { return m_m.dim_m(); } ///< Return dimension m
	constexpr std::size_t dim_n() const { return m_m.dim_n(); } ///< Return dimension n

	/** Access to the value of element [ @a i, @a j ] */
	constexpr auto operator()(std::size_t i, std::size_t j) const
	{
		const std::size_t ixs[4][3] = {
			{ 1, 2, 3 },
			{ 0, 2, 3 },
			{ 0, 1, 3 },
			{ 0, 1, 2 }
		};

		const std::size_t *ix = ixs[i];
		const std::size_t *iy = ixs[j];

		auto result =
			m_m(ix[0], iy[0]) * m_m(ix[1], iy[1]) * m_m(ix[2], iy[2]) +
			m_m(ix[0], iy[1]) * m_m(ix[1], iy[2]) * m_m(ix[2], iy[0]) +
			m_m(ix[0], iy[2]) * m_m(ix[1], iy[0]) * m_m(ix[2], iy[1]) -
			m_m(ix[0], iy[2]) * m_m(ix[1], iy[1]) * m_m(ix[2], iy[0]) -
			m_m(ix[0], iy[1]) * m_m(ix[1], iy[0]) * m_m(ix[2], iy[2]) -
			m_m(ix[0], iy[0]) * m_m(ix[1], iy[2]) * m_m(ix[2], iy[1]);

		return (i + j) % 2 == 1 ? -result : result;
	}

  private:
	const M &m_m;
};

} // namespace cif
