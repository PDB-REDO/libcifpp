//          Copyright Maarten L. Hekkelman, 2022
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#pragma once

#include <array>
#include <filesystem>
#include <fstream>
#include <memory>
#include <utility>

#include <zlib.h>

/** \file gzio.hpp
 * 
 * Single header file for the implementation of stream classes
 * that can transparently read and write compressed files.
 * 
 * The gzio::istream_buf class sniffs the input and decides whether to use
 * a decompressor if a signature was recognized.
 * 
 * There's also an ifstream and ofstream class here that can
 * read and write compressed files. In this case the decission
 * whether to use a compressions/decompression algorithm is
 * based on the extension of the \a filename argument.
 *
 * This is a stripped down version of the gxrio library from
 * https://github.com/mhekkel/gxrio.git
 * Most notably, the lzma support has been removed since getting
 * that to work in Windows proved to be too much work.
 */

namespace cif::gzio
{

/** The default buffer size to use */
const std::size_t kDefaultBufferSize = 256;

// --------------------------------------------------------------------

/// \brief A base class for the streambuf classes in gzio
///
/// \tparam CharT Type of the character stream.
/// \tparam Traits Traits for character type, defaults to char_traits<_CharT>.
///
/// The base class for all streambuf classes in this library.
/// It maintains the pointer to the upstream streambuf.

template <typename CharT, typename Traits>
class basic_streambuf : public std::basic_streambuf<CharT, Traits>
{
  public:
	/** @cond */

	using char_type = CharT;
	using traits_type = Traits;

	using int_type = typename traits_type::int_type;
	using pos_type = typename traits_type::pos_type;
	using off_type = typename traits_type::off_type;

	using streambuf_type = std::basic_streambuf<CharT, Traits>;

	basic_streambuf() = default;

	basic_streambuf(const basic_streambuf &) = delete;

	basic_streambuf(basic_streambuf &&rhs)
		: streambuf_type(std::move(rhs))
	{
		m_upstream = std::exchange(rhs.m_upstream, nullptr);
	}

	basic_streambuf &operator=(const basic_streambuf &) = delete;

	basic_streambuf &operator=(basic_streambuf &&rhs)
	{
		m_upstream = std::exchange(rhs.m_upstream, nullptr);
		return *this;
	}

	/** @endcond */

	/** Set the upstream streambuf to @a upstream */
	void set_upstream(streambuf_type *upstream)
	{
		m_upstream = upstream;
	}

	/** @cond */

	virtual basic_streambuf *init(streambuf_type *sb) = 0;
	virtual basic_streambuf *close() = 0;

	/** @endcond */

  protected:
	/// \brief The upstream streambuf object, usually this is a basic_filebuf
	streambuf_type *m_upstream = nullptr;
};

// --------------------------------------------------------------------

/// \brief A streambuf class that can be used to decompress gzipped data
///
/// \tparam CharT		Type of the character stream.
/// \tparam Traits		Traits for character type, defaults to char_traits<_CharT>.
/// \tparam BufferSize	The size of the internal buffers.
///
/// This implementation of streambuf can decompress (inflate) data compressed
/// using zlib.

template <typename CharT, typename Traits, std::size_t BufferSize = kDefaultBufferSize>
class basic_igzip_streambuf : public basic_streambuf<CharT, Traits>
{
  public:
	/** @cond */

	static_assert(sizeof(CharT) == 1, "Unfortunately, support for wide characters is not implemented yet.");

	using char_type = CharT;
	using traits_type = Traits;

	using streambuf_type = std::basic_streambuf<char_type, traits_type>;
	using base_type = basic_streambuf<CharT, Traits>;

	using int_type = typename traits_type::int_type;
	using pos_type = typename traits_type::pos_type;
	using off_type = typename traits_type::off_type;

	basic_igzip_streambuf() = default;

	basic_igzip_streambuf(const basic_igzip_streambuf &) = delete;

	/** @endcond */

	/// \brief Move constructor
	basic_igzip_streambuf(basic_igzip_streambuf &&rhs)
		: base_type(std::move(rhs))
	{
		std::swap(m_zstream, rhs.m_zstream);
		std::swap(m_gzheader, rhs.m_gzheader);

		auto p = std::copy(rhs.gptr(), rhs.egptr(), m_out_buffer.data());
		this->setg(m_out_buffer.data(), m_out_buffer.data(), p);

		if (m_zstream and m_zstream->avail_in > 0)
		{
			auto next_in_offset = m_zstream->next_in - rhs.m_in_buffer.data();
			std::copy(rhs.m_in_buffer.data() + next_in_offset,
				rhs.m_in_buffer.data() + next_in_offset + m_zstream->avail_in,
				m_in_buffer.data());
			m_zstream->next_in = m_in_buffer.data();
		}
	}

	/** @cond */

	basic_igzip_streambuf &operator=(const basic_igzip_streambuf &) = delete;

	/// \brief Move operator= implementation
	basic_igzip_streambuf &operator=(basic_igzip_streambuf &&rhs)
	{
		base_type::operator=(std::move(rhs));

		std::swap(m_zstream, rhs.m_zstream);
		std::swap(m_gzheader, rhs.m_gzheader);

		auto p = std::copy(rhs.gptr(), rhs.egptr(), m_out_buffer.data());
		this->setg(m_out_buffer.data(), m_out_buffer.data(), p);

		if (m_zstream and m_zstream->avail_in > 0)
		{
			auto next_in_offset = m_zstream->next_in - reinterpret_cast<unsigned char *>(rhs.m_in_buffer.data());
			std::copy(rhs.m_in_buffer.data() + next_in_offset,
				rhs.m_in_buffer.data() + next_in_offset + m_zstream->avail_in,
				m_in_buffer.data());
			m_zstream->next_in = reinterpret_cast<unsigned char *>(m_in_buffer.data());
		}

		return *this;
	}

	~basic_igzip_streambuf()
	{
		close();
	}

	/** @endcond */

	/// \brief This closes the zlib stream and sets the get pointers to null.
	base_type *close() override
	{
		if (m_zstream)
		{
			::inflateEnd(m_zstream.get());

			m_zstream.reset(nullptr);
			m_gzheader.reset(nullptr);
		}

		this->setg(nullptr, nullptr, nullptr);

		return this;
	}

	/// \brief Initialize a zlib stream and set the upstream.
	///
	/// \param upstream The upstream streambuf
	///
	/// The zstream is constructed and an optional header is
	/// read from upstream. The contents of the header are ignored
	/// but we must maintain that structure.
	base_type *init(streambuf_type *upstream) override
	{
		this->set_upstream(upstream);

		close();

		m_zstream.reset(new z_stream_s);
		m_gzheader.reset(new gz_header_s);

		auto &zstream = *m_zstream.get();
		zstream = z_stream_s{};
		auto &header = *m_gzheader.get();
		header = gz_header_s{};

		int err = ::inflateInit2(&zstream, 47);
		if (err == Z_OK)
		{
			zstream.next_in = reinterpret_cast<unsigned char *>(m_in_buffer.data());
			zstream.avail_in = static_cast<uInt>(this->m_upstream->sgetn(m_in_buffer.data(), m_in_buffer.size()));

			err = ::inflateGetHeader(&zstream, &header);

			if (err != Z_OK)
				::inflateEnd(&zstream);
		}

		if (err != Z_OK)
			zstream = z_stream_s{};

		return err == Z_OK ? this : nullptr;
	}

  private:
	/// \brief The actual work is done here.
	int_type underflow() override
	{
		if (m_zstream and this->m_upstream)
		{
			auto &zstream = *m_zstream.get();
			const std::streamsize kBufferByteSize = m_out_buffer.size();

			while (this->gptr() == this->egptr())
			{
				zstream.next_out = reinterpret_cast<unsigned char *>(m_out_buffer.data());
				zstream.avail_out = static_cast<uInt>(kBufferByteSize);

				if (zstream.avail_in == 0)
				{
					zstream.next_in = reinterpret_cast<unsigned char *>(m_in_buffer.data());
					zstream.avail_in = static_cast<uInt>(this->m_upstream->sgetn(m_in_buffer.data(), m_in_buffer.size()));
				}

				if (zstream.avail_in == 0)
					break;

				int err = ::inflate(&zstream, Z_SYNC_FLUSH);
				std::streamsize n = kBufferByteSize - zstream.avail_out;

				if (n > 0)
				{
					this->setg(
						m_out_buffer.data(),
						m_out_buffer.data(),
						m_out_buffer.data() + n);
					break;
				}

				if (err == Z_STREAM_END and zstream.avail_in > 0)
					err = ::inflateReset2(&zstream, 47);

				if (err < Z_OK)
					break;
			}
		}

		return this->gptr() != this->egptr() ? traits_type::to_int_type(*this->gptr()) : traits_type::eof();
	}

  private:
	/// \brief The zlib internal structures are mainained as pointers to avoid having
	/// to copy their content in move constructors.
	std::unique_ptr<z_stream_s> m_zstream;

	/// \brief The zlib internal structures are mainained as pointers to avoid having
	/// to copy their content in move constructors.
	std::unique_ptr<gz_header> m_gzheader;

	/// \brief Input buffer, this is the input for zlib
	std::array<char_type, BufferSize> m_in_buffer;

	/// \brief Output buffer, where the ostream finds the data
	std::array<char_type, BufferSize> m_out_buffer;
};

// --------------------------------------------------------------------

/// \brief A streambuf class that can be used to compress data using zlib
///
/// \tparam CharT		Type of the character stream.
/// \tparam Traits		Traits for character type, defaults to char_traits<_CharT>.
/// \tparam BufferSize	The size of the internal buffers.
///
/// This implementation of streambuf can compress (deflate) data using zlib.

template <typename CharT, typename Traits, std::size_t BufferSize = kDefaultBufferSize>
class basic_ogzip_streambuf : public basic_streambuf<CharT, Traits>
{
  public:
	/** @cond */

	static_assert(sizeof(CharT) == 1, "Unfortunately, support for wide characters is not implemented yet.");

	using char_type = CharT;
	using traits_type = Traits;

	using streambuf_type = std::basic_streambuf<char_type, traits_type>;
	using base_type = basic_streambuf<CharT, Traits>;

	using int_type = typename traits_type::int_type;
	using pos_type = typename traits_type::pos_type;
	using off_type = typename traits_type::off_type;

	basic_ogzip_streambuf() = default;

	basic_ogzip_streambuf(const basic_ogzip_streambuf &) = delete;

	/// \brief Move constructor
	basic_ogzip_streambuf(basic_ogzip_streambuf &&rhs)
		: base_type(std::move(rhs))
	{
		std::swap(m_zstream, rhs.m_zstream);
		std::swap(m_gzheader, rhs.m_gzheader);

		this->setp(m_in_buffer.data(), m_in_buffer.data() + m_in_buffer.size());
		this->sputn(rhs.pbase(), rhs.pptr() - rhs.pbase());
		rhs.setp(nullptr, nullptr);
	}

	basic_ogzip_streambuf &operator=(const basic_ogzip_streambuf &) = delete;

	/** @endcond */

	/// \brief Move operator=
	basic_ogzip_streambuf &operator=(basic_ogzip_streambuf &&rhs)
	{
		base_type::operator=(std::move(rhs));

		std::swap(m_zstream, rhs.m_zstream);
		std::swap(m_gzheader, rhs.m_gzheader);

		this->setp(m_in_buffer.data(), m_in_buffer.data() + m_in_buffer.size());
		this->sputn(rhs.pbase(), rhs.pptr() - rhs.pbase());
		rhs.setp(nullptr, nullptr);

		return *this;
	}

	~basic_ogzip_streambuf()
	{
		close();
	}

	/// \brief This closes the zlib stream and sets the put pointers to null.
	base_type *close() override
	{
		if (m_zstream)
		{
			overflow(traits_type::eof());

			::deflateEnd(m_zstream.get());

			m_zstream.reset(nullptr);
			m_gzheader.reset(nullptr);
		}

		this->setp(nullptr, nullptr);

		return this;
	}

	/// \brief Initialize the internal zlib structures
	///
	/// \param upstream The upstream streambuf
	///
	/// The zlib stream is initialized as one that can accept
	/// a gzip header.
	base_type *init(streambuf_type *upstream) override
	{
		this->set_upstream(upstream);

		close();

		m_zstream.reset(new z_stream_s);
		m_gzheader.reset(new gz_header_s);

		auto &zstream = *m_zstream.get();
		zstream = z_stream_s{};
		auto &header = *m_gzheader.get();
		header = gz_header_s{};

		const int WINDOW_BITS = 15, GZIP_ENCODING = 16;

		int err = deflateInit2(&zstream, Z_BEST_COMPRESSION, Z_DEFLATED,
			WINDOW_BITS | GZIP_ENCODING, Z_DEFLATED, Z_DEFAULT_STRATEGY);

		if (err == Z_OK)
			err = ::deflateSetHeader(&zstream, &header);

		if (err == Z_OK)
			this->setp(this->m_in_buffer.data(), this->m_in_buffer.data() + this->m_in_buffer.size());
		else
			zstream = z_stream_s{};

		return err == Z_OK ? this : nullptr;
	}

  private:
	/// \brief The actual work is done here
	///
	/// \param ch The character that did not fit, in case it is eof we need to flush
	///
	int_type overflow(int_type ch) override
	{
		if (not m_zstream)
			return traits_type::eof();

		auto &zstream = *m_zstream;

		zstream.next_in = reinterpret_cast<unsigned char *>(this->pbase());
		zstream.avail_in = static_cast<uInt>(this->pptr() - this->pbase());

		char_type buffer[BufferSize];

		for (;;)
		{
			zstream.next_out = reinterpret_cast<unsigned char *>(buffer);
			zstream.avail_out = sizeof(buffer);

			int err = ::deflate(&zstream, ch == traits_type::eof() ? Z_FINISH : Z_NO_FLUSH);

			std::streamsize n = sizeof(buffer) - zstream.avail_out;
			if (n > 0)
			{
				auto r = this->m_upstream->sputn(reinterpret_cast<char_type *>(buffer), n);

				if (r != n)
					return traits_type::eof();
			}

			if (zstream.avail_out == 0)
				continue;

			if (err == Z_OK and ch == traits_type::eof())
				continue;

			break;
		}

		this->setp(this->m_in_buffer.data(), this->m_in_buffer.data() + this->m_in_buffer.size());

		if (not traits_type::eq_int_type(ch, traits_type::eof()))
		{
			*this->pptr() = traits_type::to_char_type(ch);
			this->pbump(1);
		}

		return ch;
	}

  private:
	/// \brief The zlib internal structures are mainained as pointers to avoid having
	/// to copy their content in move constructors.
	std::unique_ptr<z_stream_s> m_zstream;

	/// \brief The zlib internal structures are mainained as pointers to avoid having
	/// to copy their content in move constructors.
	std::unique_ptr<gz_header> m_gzheader;

	/// \brief Input buffer, this is the input for zlib
	std::array<char_type, BufferSize> m_in_buffer;
};

// --------------------------------------------------------------------

/// \brief An istream implementation that wraps a streambuf with a decompressing streambuf
///
/// \tparam CharT		Type of the character stream.
/// \tparam Traits		Traits for character type, defaults to char_traits<_CharT>.
///
/// This is an istream implementation that can take a source streambuf and then wraps
/// this streambuf with a decompressing streambuf class defined above.
/// The class inherits from std::basic_istream and offers all the associated functionality.

template <typename CharT, typename Traits>
class basic_istream : public std::basic_istream<CharT, Traits>
{
  public:
	/** @cond */

	using base_type = std::basic_istream<CharT, Traits>;

	using traits_type = Traits;
	using char_type = typename traits_type::char_type;
	using int_type = typename traits_type::int_type;

	using z_streambuf_type = basic_streambuf<char_type, traits_type>;
	using upstreambuf_type = std::basic_streambuf<char_type, traits_type>;

	using gzip_streambuf_type = basic_igzip_streambuf<char_type, traits_type>;

	/** @endcond */

	/// \brief Regular move constructor
	basic_istream(basic_istream &&rhs)
		: base_type(std::move(rhs))
	{
		m_gziobuf = std::move(rhs.m_gziobuf);

		if (m_gziobuf)
			this->rdbuf(m_gziobuf.get());
		else
			this->rdbuf(nullptr);
	}

	/// \brief Regular move operator=
	basic_istream &operator=(basic_istream &&rhs)
	{
		base_type::operator=(std::move(rhs));
		m_gziobuf = std::move(rhs.m_gziobuf);

		if (m_gziobuf)
			this->rdbuf(m_gziobuf.get());
		else
			this->rdbuf(nullptr);

		return *this;
	}

	/// \brief Construct an istream with the passed in streambuf \a buf
	///
	/// \param buf The streambuf that provides the compressed data
	///
	/// This constructor will initialize the zlib code with the \a buf streambuf.

	explicit basic_istream(upstreambuf_type *buf)
		: base_type(nullptr)
	{
		init_z(buf);
	}

  protected:
	basic_istream()
		: base_type(nullptr) {}

	/// \brief Initialise internals with streambuf \a sb
	/// \param sb The upstream streambuf class
	///
	/// This will sniff the content in \a sb and decide upon what is found
	/// what implementation is used. If it doesn't look like compressed data
	/// the \a sb streambuf is used without any decompression being done.

	void init_z(upstreambuf_type *sb)
	{
		int_type ch = sb->sgetc();
		if (ch == 0x1f)
		{
			sb->sbumpc();
			ch = sb->sgetc();
			sb->sungetc();

			if (ch == 0x8b) // Read gzip header
				m_gziobuf.reset(new gzip_streambuf_type);
		}

		if (m_gziobuf)
		{
			if (not m_gziobuf->init(sb))
				this->setstate(std::ios_base::failbit);
			else
				this->init(m_gziobuf.get());
		}
		else
			this->init(sb);
	}

  protected:
	/// \brief Our streambuf class
	std::unique_ptr<z_streambuf_type> m_gziobuf;
};

// --------------------------------------------------------------------

/// \brief Control input from files compressed with gzip.
///
/// \tparam CharT		Type of the character stream.
/// \tparam Traits		Traits for character type, defaults to char_traits<_CharT>.
///
/// This is an ifstream implementation that can read from named files compressed with
/// gzip directly. The class inherits from std::basic_istream and offers all the
/// associated functionality.

template <typename CharT, typename Traits>
class basic_ifstream : public basic_istream<CharT, Traits>
{
  public:
	/** @cond */

	using base_type = basic_istream<CharT, Traits>;

	using char_type = CharT;
	using traits_type = Traits;

	using filebuf_type = std::basic_filebuf<char_type, traits_type>;

	using gzip_streambuf_type = typename base_type::gzip_streambuf_type;

	/// \brief Default constructor, does not open a file since none is specified
	basic_ifstream() = default;

	~basic_ifstream()
	{
		close();
	}

	/** @endcond */

	/// \brief Construct an ifstream
	/// \param filename Null terminated string specifying the file to open
	/// \param mode The mode in which to open the file

	explicit basic_ifstream(const char *filename, std::ios_base::openmode mode = std::ios_base::in)
	{
		open(filename, mode);
	}

	/// \brief Construct an ifstream
	/// \param filename std::string specifying the file to open
	/// \param mode The mode in which to open the file

	explicit basic_ifstream(const std::string &filename, std::ios_base::openmode mode = std::ios_base::in)
	{
		open(filename, mode);
	}

	/// \brief Construct an ifstream
	/// \param filename std::filesystem::path specifying the file to open
	/// \param mode The mode in which to open the file

	explicit basic_ifstream(const std::filesystem::path &filename, std::ios_base::openmode mode = std::ios_base::in)
	{
		open(filename, mode);
	}

	/// \brief Move constructor
	basic_ifstream(basic_ifstream &&rhs)
		: base_type(std::move(rhs))
	{
		m_filebuf = std::move(rhs.m_filebuf);

		if (this->m_gziobuf)
			this->m_gziobuf->set_upstream(&m_filebuf);
		else
			this->rdbuf(&m_filebuf);
	}

	/** @cond */
	basic_ifstream(const basic_ifstream &) = delete;

	basic_ifstream &operator=(const basic_ifstream &) = delete;

	/** @endcond */

	/// \brief Move version of operator=
	basic_ifstream &operator=(basic_ifstream &&rhs)
	{
		base_type::operator=(std::move(rhs));

		m_filebuf = std::move(rhs.m_filebuf);
		if (this->m_gziobuf)
			this->m_gziobuf->set_upstream(&m_filebuf);
		else
			this->rdbuf(&m_filebuf);

		return *this;
	}

	/// \brief Open the file \a filename with mode \a mode
	/// \param filename std::filesystem::path specifying the file to open
	/// \param mode The mode in which to open the file

	void open(const std::filesystem::path &filename, std::ios_base::openmode mode = std::ios_base::in)
	{
		if (not m_filebuf.open(filename, mode | std::ios::binary))
			this->setstate(std::ios_base::failbit);
		else
		{
			if (filename.extension() == ".gz")
				this->m_gziobuf.reset(new gzip_streambuf_type);

			if (not this->m_gziobuf)
			{
				this->rdbuf(&m_filebuf);
				this->clear();
			}
			else if (not this->m_gziobuf->init(&m_filebuf))
				this->setstate(std::ios_base::failbit);
			else
			{
				this->rdbuf(this->m_gziobuf.get());
				this->clear();
			}
		}
	}

	/// \brief Open the file \a filename with mode \a mode
	/// \param filename std::string specifying the file to open
	/// \param mode The mode in which to open the file

	void open(const std::string &filename, std::ios_base::openmode mode = std::ios_base::in)
	{
		open(std::filesystem::path{filename}, mode);
	}

	/// \brief Open the file \a filename with mode \a mode
	/// \param filename Null terminated string specifying the file to open
	/// \param mode The mode in which to open the file

	void open(const char *filename, std::ios_base::openmode mode = std::ios_base::in)
	{
		open(std::filesystem::path{filename}, mode);
	}

	/// \brief Return true if the file is open
	/// \return m_filebuf.is_open()

	bool is_open() const
	{
		return m_filebuf.is_open();
	}

	/// \brief Close the file
	///
	/// Calls m_filebuf.close(). If that fails, the failbit is set.

	void close()
	{
		if (this->m_gziobuf and not this->m_gziobuf->close())
			this->setstate(std::ios_base::failbit);

		if (not m_filebuf.close())
			this->setstate(std::ios_base::failbit);
	}

	/// \brief Swap the contents with those of \a rhs
	/// \param rhs The ifstream to swap with

	void swap(basic_ifstream &rhs)
	{
		base_type::swap(rhs);
		m_filebuf.swap(rhs.m_filebuf);

		if (this->m_gziobuf)
		{
			this->m_gziobuf.set_upstream(&m_filebuf);
			this->rdbuf(this->m_gziobuf.get());
		}
		else
			this->rdbuf(&m_filebuf);

		if (rhs.m_gziobuf)
		{
			rhs.m_gziobuf.set_upstream(&rhs.m_filebuf);
			rhs.rdbuf(rhs.m_gziobuf.get());
		}
		else
			rhs.rdbuf(&rhs.m_filebuf);
	}

  private:
	/// \brief The filebuf
	filebuf_type m_filebuf;
};

// --------------------------------------------------------------------

/// \brief An ostream implementation that wraps a streambuf with a compressing streambuf
///
/// \tparam CharT		Type of the character stream.
/// \tparam Traits		Traits for character type, defaults to char_traits<_CharT>.
///
/// This is an ostream implementation that can take an upstream streambuf and then wraps
/// this streambuf with a compressing streambuf class defined above.
/// The class inherits from std::basic_ostream and offers all the associated functionality.

template <typename CharT, typename Traits>
class basic_ostream : public std::basic_ostream<CharT, Traits>
{
  public:
	/** @cond */

	using base_type = std::basic_ostream<CharT, Traits>;

	using char_type = CharT;
	using traits_type = Traits;

	using z_streambuf_type = basic_streambuf<char_type, traits_type>;
	using upstreambuf_type = std::basic_streambuf<char_type, traits_type>;

	/** @endcond */

	/// \brief Regular move constructor
	basic_ostream(basic_ostream &&rhs)
		: base_type(std::move(rhs))
	{
		m_gziobuf = std::move(rhs.m_gziobuf);
		this->rdbuf(m_gziobuf.get());
	}

	/// \brief Regular move operator=
	basic_ostream &operator=(basic_ostream &&rhs)
	{
		base_type::operator=(std::move(rhs));
		m_gziobuf = std::move(rhs.m_gziobuf);

		this->rdbuf(m_gziobuf.get());

		return *this;
	}

	// One might expect a constructor taking a streambuf pointer
	// as the regular std::ostream class does. However, that's not
	// going to work since no information is known at this time
	// what compression to use.
	//
	// explicit basic_ostream(upstreambuf_type *buf)
	// {
	// 	init_z(buf);
	// 	this->init(m_gziobuf.get());
	// }

  protected:
	basic_ostream()
		: base_type(nullptr) {}

	/// \brief Initialise internals with streambuf \a sb
	void init_z(std::streambuf *sb)
	{
		if (not m_gziobuf->init(sb))
			this->setstate(std::ios_base::failbit);
	}

  protected:
	/// \brief Our streambuf class
	std::unique_ptr<z_streambuf_type> m_gziobuf;
};

// --------------------------------------------------------------------

/// \brief Control output to files compressing the contents with gzip.
///
/// \tparam CharT		Type of the character stream.
/// \tparam Traits		Traits for character type, defaults to char_traits<_CharT>.
///
/// This is an ofstream implementation that can writeto named files compressing the content
/// with gzip directly. The class inherits from std::basic_ostream and offers all the
/// associated functionality.

template <typename CharT, typename Traits>
class basic_ofstream : public basic_ostream<CharT, Traits>
{
  public:
	/** @cond */

	using base_type = basic_ostream<CharT, Traits>;

	using char_type = CharT;
	using traits_type = Traits;

	using filebuf_type = std::basic_filebuf<char_type, traits_type>;
	using gzip_streambuf_type = basic_ogzip_streambuf<char_type, traits_type>;

	basic_ofstream() = default;

	~basic_ofstream()
	{
		close();
	}

	/** @endcond */

	/// \brief Construct an ofstream
	/// \param filename Null terminated string specifying the file to open
	/// \param mode The mode in which to open the file

	explicit basic_ofstream(const char *filename, std::ios_base::openmode mode = std::ios_base::out)
	{
		open(filename, mode);
	}

	/// \brief Construct an ofstream
	/// \param filename std::string specifying the file to open
	/// \param mode The mode in which to open the file

	explicit basic_ofstream(const std::string &filename, std::ios_base::openmode mode = std::ios_base::out)
	{
		open(filename, mode);
	}

	/// \brief Construct an ofstream
	/// \param filename std::filesystem::path specifying the file to open
	/// \param mode The mode in which to open the file

	explicit basic_ofstream(const std::filesystem::path &filename, std::ios_base::openmode mode = std::ios_base::out)
	{
		open(filename, mode);
	}

	/// \brief Move constructor
	basic_ofstream(basic_ofstream &&rhs)
		: base_type(std::move(rhs))
	{
		m_filebuf = std::move(rhs.m_filebuf);
		if (this->m_gziobuf)
			this->m_gziobuf->set_upstream(&m_filebuf);
		else
			this->rdbuf(&m_filebuf);
	}

	/** @cond */

	basic_ofstream(const basic_ofstream &) = delete;

	basic_ofstream &operator=(const basic_ofstream &) = delete;

	/** @endcond */

	/// \brief Move operator=
	basic_ofstream &operator=(basic_ofstream &&rhs)
	{
		base_type::operator=(std::move(rhs));
		m_filebuf = std::move(rhs.m_filebuf);
		if (this->m_gziobuf)
			this->m_gziobuf->set_upstream(&m_filebuf);
		else
			this->rdbuf(&m_filebuf);

		return *this;
	}

	/// \brief Open the file \a filename with mode \a mode
	/// \param filename std::filesystem::path specifying the file to open
	/// \param mode The mode in which to open the file
	///
	/// A compression algorithm is chosen upon the contents of the
	/// extension() of \a filename with .gz mapping to gzip compression
	/// and .xz to xz compression.

	void open(const std::filesystem::path &filename, std::ios_base::openmode mode = std::ios_base::out)
	{
		if (not m_filebuf.open(filename, mode | std::ios::binary))
			this->setstate(std::ios_base::failbit);
		else
		{
			if (filename.extension() == ".gz")
				this->m_gziobuf.reset(new gzip_streambuf_type);
			else
				this->m_gziobuf.reset(nullptr);

			if (this->m_gziobuf)
			{
				if (not this->m_gziobuf->init(&m_filebuf))
					this->setstate(std::ios_base::failbit);
				else
				{
					this->rdbuf(this->m_gziobuf.get());
					this->clear();
				}
			}
			else
			{
				this->rdbuf(&m_filebuf);
				this->clear();
			}
		}
	}

	/// \brief Open the file \a filename with mode \a mode
	/// \param filename std::string specifying the file to open
	/// \param mode The mode in which to open the file
	///
	/// A compression algorithm is chosen upon the contents of the
	/// extension of \a filename with .gz mapping to gzip compression
	/// and .xz to xz compression.

	void open(const std::string &filename, std::ios_base::openmode mode = std::ios_base::out)
	{
		open(std::filesystem::path{filename}, mode);
	}

	/// \brief Open the file \a filename with mode \a mode
	/// \param filename Null terminated string specifying the file to open
	/// \param mode The mode in which to open the file
	///
	/// A compression algorithm is chosen upon the contents of the
	/// extension of \a filename with .gz mapping to gzip compression
	/// and .xz to xz compression.

	void open(const char *filename, std::ios_base::openmode mode = std::ios_base::out)
	{
		open(std::filesystem::path{filename}, mode);
	}

	/// \brief Return true if the file is open
	/// \return m_filebuf.is_open()

	bool is_open() const
	{
		return m_filebuf.is_open();
	}

	/// \brief Close the file
	///
	/// Calls m_filebuf.close(). If that fails, the failbit is set.

	void close()
	{
		if (this->m_gziobuf and not this->m_gziobuf->close())
			this->setstate(std::ios_base::failbit);

		if (not m_filebuf.close())
			this->setstate(std::ios_base::failbit);
	}

	/// \brief Swap the contents with those of \a rhs
	/// \param rhs The ifstream to swap with

	void swap(basic_ofstream &rhs)
	{
		base_type::swap(rhs);
		m_filebuf.swap(rhs.m_filebuf);

		if (this->m_gziobuf)
		{
			this->m_gziobuf.set_upstream(&m_filebuf);
			this->rdbuf(this->m_gziobuf.get());
		}
		else
			this->rdbuf(&m_filebuf);

		if (rhs.m_gziobuf)
		{
			rhs.m_gziobuf.set_upstream(&rhs.m_filebuf);
			rhs.rdbuf(rhs.m_gziobuf.get());
		}
		else
			rhs.rdbuf(&rhs.m_filebuf);
	}

  private:
	/// \brief The filebuf
	filebuf_type m_filebuf;
};

// --------------------------------------------------------------------

/// \brief Convenience typedef for a regular istream
using istream = basic_istream<char, std::char_traits<char>>;

/// \brief Convenience typedef for a file ifstream
using ifstream = basic_ifstream<char, std::char_traits<char>>;

/// \brief Convenience typedef for a file ofstream
using ofstream = basic_ofstream<char, std::char_traits<char>>;

} // namespace gzio
