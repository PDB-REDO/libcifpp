/*-
 * SPDX-License-Identifier: BSD-2-Clause
 *
 * Copyright (c) 2022 NKI/AVL, Netherlands Cancer Institute
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

#include <list>

#include "cif++/datablock.hpp"
#include "cif++/parser.hpp"

/** \file file.hpp
 * 
 * The file class defined here encapsulates the contents of an mmCIF file
 * It is mainly a list of @ref cif::datablock objects
 * 
 * The class file has methods to load dictionaries. These dictionaries are
 * loaded from resources (if available) or from disk from several locations.
 * 
 * See the documentation on load_resource() in file: utilities.hpp for more
 * information on how data is loaded. 
 */

namespace cif
{

// --------------------------------------------------------------------

/**
 * @brief The class file is actually a list of datablock objects
 * 
 */

class file : public std::list<datablock>
{
  public:
	file() = default;

	/**
	 * @brief Construct a new file object using the data in the file @a p as content
	 * 
	 * @param p Path to a file containing the data to load
	 */
	explicit file(const std::filesystem::path &p)
	{
		load(p);
	}

	/**
	 * @brief Construct a new file object using the data in the std::istream @a is
	 * 
	 * @param is The istream containing the data to load
	 */
	explicit file(std::istream &is)
	{
		load(is);
	}

	/**
	 * @brief Construct a new file object with data in the constant string defined
	 * by @a data and @a length
	 * 
	 * @param data The pointer to the character string with data to load
	 * @param length The length of the data
	 */
	explicit file(const char *data, size_t length)
	{
		struct membuf : public std::streambuf
		{
			membuf(char *text, size_t length)
			{
				this->setg(text, text, text + length);
			}
		} buffer(const_cast<char *>(data), length);

		std::istream is(&buffer);
		load(is);
	}

	/** @cond */
	file(const file &rhs)
		: std::list<datablock>(rhs)
	{
	}

	file(file &&rhs)
	{
		this->swap(rhs);
	}

	file &operator=(file f)
	{
		this->swap(f);
		return *this;
	}

	/** @endcond */

	/**
	 * @brief Set the validator object to @a v
	 */
	void set_validator(const validator *v);

	/**
	 * @brief Get the validator object
	 */
	const validator *get_validator() const
	{
		return m_validator;
	}

	/**
	 * @brief Validate the content and return true if everything was valid.
	 * 
	 * Will throw an exception if there is no validator defined.
	 * 
	 * If each category was valid, validate_links will also be called.
	 * 
	 * @return true If the content is valid
	 * @return false If the content is not valid
	 */
	bool is_valid() const;

	/**
	 * @brief Validate the content and return true if everything was valid.
	 * 
	 * Will attempt to load the referenced dictionary if none was specified.
	 * 
	 * If each category was valid, validate_links will also be called.
	 * 
	 * @return true If the content is valid
	 * @return false If the content is not valid
	 */
	bool is_valid();

	/**
	 * @brief Validate the links for all datablocks contained.
	 * 
	 * Will throw an exception if no validator was specified.
	 * 
	 * @return true If all links were valid
	 * @return false If all links were not valid
	 */
	bool validate_links() const;

	/**
	 * @brief Attempt to load a dictionary (validator) based on
	 * the contents of the *audit_conform* category, if available.
	 */
	void load_dictionary();


	/**
	 * @brief Attempt to load the named dictionary @a name and
	 * create a validator based on it.
	 * 
	 * @param name The name of the dictionary to load
	 */
	void load_dictionary(std::string_view name);

	/**
	 * @brief Return true if a datablock with the name @a name is part of this file
	 */
	bool contains(std::string_view name) const;

	/**
	 * @brief return a reference to the first datablock in the file
	 */
	datablock &front()
	{
		assert(not empty());
		return std::list<datablock>::front();
	}

	/**
	 * @brief return a const reference to the first datablock in the file
	 */
	const datablock &front() const
	{
		assert(not empty());
		return std::list<datablock>::front();
	}

	/**
	 * @brief return a reference to the datablock named @a name
	 */
	datablock &operator[](std::string_view name);

	/**
	 * @brief return a const reference to the datablock named @a name
	 */
	const datablock &operator[](std::string_view name) const;

	/**
	 * @brief Tries to find a datablock with name @a name and will create a
	 * new one if it is not found. The result is a tuple of an iterator
	 * pointing to the datablock and a boolean indicating whether the datablock
	 * was created or not.
	 * 
	 * @param name The name for the datablock
	 * @return std::tuple<iterator, bool> A tuple containing an iterator pointing
	 * at the datablock and a boolean indicating whether the datablock was newly
	 * created.
	 */
	std::tuple<iterator, bool> emplace(std::string_view name);

	/** Load the data from the file specified by @a p */
	void load(const std::filesystem::path &p);

	/** Load the data from @a is */
	void load(std::istream &is);

	/** Save the data to the file specified by @a p */
	void save(const std::filesystem::path &p) const;

	/** Save the data to @a is */
	void save(std::ostream &os) const;

	/**
	 * @brief Friend operator<< to write file @a f to std::ostream @a os
	 */
	friend std::ostream &operator<<(std::ostream &os, const file &f)
	{
		f.save(os);
		return os;
	}

  private:
	const validator *m_validator = nullptr;
};

} // namespace cif