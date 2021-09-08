/*-
 * SPDX-License-Identifier: BSD-2-Clause
 * 
 * Copyright Maarten L. Hekkelman, Radboud University 2008-2011.
 * Copyright (c) 2021 NKI/AVL, Netherlands Cancer Institute
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

// --------------------------------------------------------------------
// uBlas compatible matrix types

#pragma once

#include <iostream>
#include <vector>

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
std::ostream &operator<<(std::ostream &lhs, const MatrixBase<T> &rhs)
{
	lhs << '[' << rhs.dim_m() << ',' << rhs.dim_n() << ']' << '(';
	for (uint32_t i = 0; i < rhs.dim_m(); ++i)
	{
		lhs << '(';
		for (uint32_t j = 0; j < rhs.dim_n(); ++j)
		{
			if (j > 0)
				lhs << ',';
			lhs << rhs(i, j);
		}
		lhs << ')';
	}
	lhs << ')';

	return lhs;
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

// template <typename T>
// symmetric_matrix<T> hammingDistance(const MatrixBase<T> &lhs, T rhs);

// template <typename T>
// std::vector<T> sum(const MatrixBase<T> &m);
