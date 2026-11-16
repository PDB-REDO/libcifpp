// Copyright Maarten L. Hekkelman, 2022
// SPDX-FileCopyrightText: 2026 Maarten L. Hekkelman
//
// SPDX-License-Identifier: BSL-1.0

#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <ios>
#include <istream>
#include <memory>
#include <ostream>
#include <streambuf>
#include <string>
#include <type_traits>
#include <utility>

// Allow standalone use of this header: unless the build system configured
// HAVE_LibLZMA explicitly, detect liblzma by the presence of its header.
#if !defined(HAVE_LibLZMA)
# if defined(__has_include)
#  if __has_include(<lzma.h>)
#   define HAVE_LibLZMA 1
#  else
#   define HAVE_LibLZMA 0
#  endif
# else
#  define HAVE_LibLZMA 0
# endif
#endif

#include <zlib.h>
#if HAVE_LibLZMA
#include <lzma.h>
#endif

/// \file gxrio.hpp
///
/// Single header file for the implementation of stream classes
/// that can transparently read and write compressed files.
///
/// The gxrio::istream_buf class sniffs the input and decides whether to use
/// a decompressor if a signature was recognized.
///
/// There's also an ifstream and ofstream class here that can
/// read and write compressed files. In this case the decision
/// whether to use a compressions/decompression algorithm is
/// based on the extension of the \a filename argument.

namespace gxrio
{

constexpr size_t kDefaultBufferSize = 8192;
constexpr int kDefaultCompressionLevel = 9;

/// \brief The longest signature of the supported compression formats (the xz magic)
constexpr std::streamsize kMaxSignatureLength = 6;

namespace detail
{

/// \brief Read from an upstream streambuf into \a buffer, clamping a negative
/// result to zero so it can safely be used as an unsigned byte count.
template <typename Streambuf, typename CharT, size_t N>
std::streamsize sgetn(Streambuf *sb, std::array<CharT, N> &buffer)
{
	std::streamsize n = sb->sgetn(buffer.data(), static_cast<std::streamsize>(buffer.size()));
	return n < 0 ? 0 : n;
}

} // namespace detail

// --------------------------------------------------------------------

/// \brief A base class for the streambuf classes in gxrio
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
	using char_type = CharT;
	using traits_type = Traits;

	using int_type = typename traits_type::int_type;
	using pos_type = typename traits_type::pos_type;
	using off_type = typename traits_type::off_type;

	using streambuf_type = std::basic_streambuf<CharT, Traits>;

	basic_streambuf() = default;

	basic_streambuf(const basic_streambuf &) = delete;

	basic_streambuf(basic_streambuf &&rhs) noexcept
		: streambuf_type(static_cast<streambuf_type &&>(rhs))
	{
		m_upstream = std::exchange(rhs.m_upstream, nullptr);
		m_ios = std::exchange(rhs.m_ios, nullptr);
		m_has_error = rhs.m_has_error;
		rhs.m_has_error = false;
	}

	basic_streambuf &operator=(const basic_streambuf &) = delete;

	basic_streambuf &operator=(basic_streambuf &&rhs) noexcept
	{
		m_upstream = std::exchange(rhs.m_upstream, nullptr);
		m_ios = std::exchange(rhs.m_ios, nullptr);
		m_has_error = rhs.m_has_error;
		rhs.m_has_error = false;
		return *this;
	}

	void set_upstream(streambuf_type *upstream) noexcept
	{
		m_upstream = upstream;
	}

	/// \brief Register the stream that should be notified of (de)compression errors
	///
	/// The stream is expected to be the stream that owns this streambuf as its
	/// rdbuf. It receives badbit when an error is encountered.
	void set_ios(std::basic_ios<CharT, Traits> *ios) noexcept
	{
		m_ios = ios;
	}

	/// \brief Return true if a (de)compression error was encountered
	[[nodiscard]] bool error() const noexcept
	{
		return m_has_error;
	}

	virtual basic_streambuf *init(streambuf_type *sb) = 0;
	virtual basic_streambuf *close() = 0;

	/// \brief Set the compression level, only meaningful for the compressing streambufs
	virtual void set_compression_level(int /*level*/) noexcept {}

  protected:
	/// \brief Mark the streambuf as failed and set badbit on the owning stream
	///
	/// Error reporting is noexcept on purpose: this method is reached from
	/// close(), which in turn is called from destructors, and from underflow()
	/// and overflow(). Throwing here would either terminate the process or
	/// interrupt an in-progress close(). The error remains visible through
	/// error() and through the stream's badbit.
	void set_error() noexcept
	{
		m_has_error = true;

		if (m_ios != nullptr)
		{
			try
			{
				m_ios->setstate(std::ios_base::badbit);
			}
			catch (...)
			{
				// badbit has been set; swallow the notification so error
				// reporting (and close from a destructor) cannot throw.
				m_has_error = true;
			}
		}
	}

	/// \brief The upstream streambuf object, usually this is a basic_filebuf
	streambuf_type *m_upstream = nullptr;

	/// \brief The stream to notify when an error occurs, may be null
	std::basic_ios<CharT, Traits> *m_ios = nullptr;

	/// \brief Set when a (de)compression error was encountered
	bool m_has_error = false;
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

template <typename CharT, typename Traits, size_t BufferSize = kDefaultBufferSize, std::enable_if_t<sizeof(CharT) == 1, int> = 0>
class basic_igzip_streambuf : public basic_streambuf<CharT, Traits>
{
  public:
	using char_type = CharT;
	using traits_type = Traits;

	using streambuf_type = std::basic_streambuf<char_type, traits_type>;
	using base_type = basic_streambuf<CharT, Traits>;

	using int_type = typename traits_type::int_type;
	using pos_type = typename traits_type::pos_type;
	using off_type = typename traits_type::off_type;

	basic_igzip_streambuf() = default;

	basic_igzip_streambuf(const basic_igzip_streambuf &) = delete;

	/// \brief Move constructor
	basic_igzip_streambuf(basic_igzip_streambuf &&rhs) noexcept
		: base_type(static_cast<base_type &&>(rhs))
	{
		std::swap(m_zstream, rhs.m_zstream);
		std::swap(m_gzheader, rhs.m_gzheader);
		m_finished = rhs.m_finished;

		auto p = m_out_buffer.data();
		if (rhs.gptr() != nullptr and rhs.egptr() != nullptr)
			p = std::copy(rhs.gptr(), rhs.egptr(), m_out_buffer.data());
		this->setg(m_out_buffer.data(), m_out_buffer.data(), p);

		if (m_zstream and m_zstream->avail_in > 0)
		{
			auto next_in_offset = m_zstream->next_in - reinterpret_cast<unsigned char *>(rhs.m_in_buffer.data());
			std::copy(rhs.m_in_buffer.data() + next_in_offset,
				rhs.m_in_buffer.data() + next_in_offset + m_zstream->avail_in,
				m_in_buffer.data());
			m_zstream->next_in = reinterpret_cast<unsigned char *>(m_in_buffer.data());
		}
	}

	basic_igzip_streambuf &operator=(const basic_igzip_streambuf &) = delete;

	/// \brief Move operator= implementation
	basic_igzip_streambuf &operator=(basic_igzip_streambuf &&rhs) noexcept
	{
		base_type::operator=(static_cast<base_type &&>(rhs));

		std::swap(m_zstream, rhs.m_zstream);
		std::swap(m_gzheader, rhs.m_gzheader);
		m_finished = rhs.m_finished;

		auto p = m_out_buffer.data();
		if (rhs.gptr() != nullptr and rhs.egptr() != nullptr)
			p = std::copy(rhs.gptr(), rhs.egptr(), m_out_buffer.data());
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

	~basic_igzip_streambuf() override
	{
		close();
	}

	/// \brief This closes the zlib stream and sets the get pointers to null.
	base_type *close() noexcept override
	{
		if (m_zstream)
		{
			::inflateEnd(m_zstream.get());

			m_zstream.reset(nullptr);
			m_gzheader.reset(nullptr);
		}

		m_finished = false;
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
		if (upstream == nullptr)
			return nullptr;

		this->set_upstream(upstream);

		close();

		m_finished = false;
		this->m_has_error = false;

		m_zstream = std::make_unique<z_stream_s>();
		m_gzheader = std::make_unique<gz_header_s>();

		auto &zstream = *m_zstream.get();
		zstream = z_stream_s{};
		auto &header = *m_gzheader.get();
		header = gz_header_s{};

		int err = ::inflateInit2(&zstream, 47);
		if (err == Z_OK)
		{
			zstream.next_in = reinterpret_cast<unsigned char *>(m_in_buffer.data());
			zstream.avail_in = static_cast<uInt>(detail::sgetn(this->m_upstream, m_in_buffer));

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
		if (this->m_has_error)
			return traits_type::eof();

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
					zstream.avail_in = static_cast<uInt>(detail::sgetn(this->m_upstream, m_in_buffer));
				}

				if (zstream.avail_in == 0)
				{
					// The upstream is exhausted. If the compressed stream did not
					// end cleanly, the input was truncated.
					if (not m_finished)
						this->set_error();
					break;
				}

				// A previous call saw the end of a gzip member while there was
				// still input available, i.e. a concatenated member. Reset for it.
				if (m_finished)
				{
					if (::inflateReset2(&zstream, 47) != Z_OK)
					{
						this->set_error();
						break;
					}
					m_finished = false;
				}

				int err = ::inflate(&zstream, Z_NO_FLUSH);
				std::streamsize n = kBufferByteSize - zstream.avail_out;

				if (err == Z_STREAM_END)
				{
					if (zstream.avail_in > 0)
					{
						// Another gzip member follows in this input buffer.
						if (::inflateReset2(&zstream, 47) != Z_OK)
						{
							this->set_error();
							break;
						}
						err = Z_OK;
					}
					else
						m_finished = true;
				}
				else if (err < Z_OK)
				{
					this->set_error();
					break;
				}

				if (n > 0)
				{
					this->setg(
						m_out_buffer.data(),
						m_out_buffer.data(),
						m_out_buffer.data() + n);
					break;
				}

				// No output was produced. If this was not a clean end of the
				// stream, try again with more input or output space.
				if (err == Z_STREAM_END)
					continue;
			}
		}

		return this->gptr() != this->egptr() ? traits_type::to_int_type(*this->gptr()) : traits_type::eof();
	}

  private:
	/// \brief The zlib internal structures are maintained as pointers to avoid having
	/// to copy their content in move constructors.
	std::unique_ptr<z_stream_s> m_zstream;

	/// \brief The zlib internal structures are maintained as pointers to avoid having
	/// to copy their content in move constructors.
	std::unique_ptr<gz_header> m_gzheader;

	/// \brief Set when the end of a gzip member was reached
	bool m_finished = false;

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

template <typename CharT, typename Traits, size_t BufferSize = kDefaultBufferSize, std::enable_if_t<sizeof(CharT) == 1, int> = 0>
class basic_ogzip_streambuf : public basic_streambuf<CharT, Traits>
{
  public:
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
	basic_ogzip_streambuf(basic_ogzip_streambuf &&rhs) noexcept
		: base_type(static_cast<base_type &&>(rhs))
	{
		std::swap(m_zstream, rhs.m_zstream);
		std::swap(m_gzheader, rhs.m_gzheader);

		m_compression_level = rhs.m_compression_level;

		this->setp(m_in_buffer.data(), m_in_buffer.data() + m_in_buffer.size());
		if (rhs.pbase() != nullptr and rhs.pptr() != nullptr)
			this->sputn(rhs.pbase(), rhs.pptr() - rhs.pbase());
		rhs.setp(nullptr, nullptr);
	}

	basic_ogzip_streambuf &operator=(const basic_ogzip_streambuf &) = delete;

	/// \brief Move operator=
	basic_ogzip_streambuf &operator=(basic_ogzip_streambuf &&rhs) noexcept
	{
		if (this == &rhs)
			return *this;

		// Flush the old state before the base move discards our upstream.
		close();

		base_type::operator=(static_cast<base_type &&>(rhs));

		std::swap(m_zstream, rhs.m_zstream);
		std::swap(m_gzheader, rhs.m_gzheader);

		m_compression_level = rhs.m_compression_level;

		this->setp(m_in_buffer.data(), m_in_buffer.data() + m_in_buffer.size());
		if (rhs.pbase() != nullptr and rhs.pptr() != nullptr)
			this->sputn(rhs.pbase(), rhs.pptr() - rhs.pbase());
		rhs.setp(nullptr, nullptr);

		return *this;
	}

	~basic_ogzip_streambuf() override
	{
		close();
	}

	/// \brief This closes the zlib stream and sets the put pointers to null.
	///
	/// Returns nullptr when flushing the remaining data failed, so that
	/// callers can propagate the error.
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

		return this->m_has_error ? nullptr : this;
	}

	/// \brief Set the compression level, 0 (no compression) to 9 (maximum compression)
	///
	/// When the stream is already initialized the new level is applied
	/// immediately using deflateParams().
	void set_compression_level(int level) noexcept override
	{
		m_compression_level = std::clamp(level, 0, 9);

		if (m_zstream)
			::deflateParams(m_zstream.get(), m_compression_level, Z_DEFAULT_STRATEGY);
	}

	/// \brief Initialize the internal zlib structures
	///
	/// \param upstream The upstream streambuf
	///
	/// The zlib stream is initialized as one that can accept
	/// a gzip header.
	base_type *init(streambuf_type *upstream) override
	{
		if (upstream == nullptr)
			return nullptr;

		// Finish and tear down any existing stream while the old upstream is
		// still installed. Otherwise buffered data would be flushed to the new
		// upstream. Only then start a clean session on the new upstream.
		close();

		this->set_upstream(upstream);
		this->m_has_error = false;

		m_zstream = std::make_unique<z_stream_s>();
		m_gzheader = std::make_unique<gz_header_s>();

		auto &zstream = *m_zstream.get();
		zstream = z_stream_s{};
		auto &header = *m_gzheader.get();
		header = gz_header_s{};

		const int WINDOW_BITS = 15, GZIP_ENCODING = 16;

		int err = deflateInit2(&zstream, m_compression_level, Z_DEFLATED,
			WINDOW_BITS | GZIP_ENCODING, 8 /* memLevel */, Z_DEFAULT_STRATEGY);

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
		if (not m_zstream or this->m_has_error or this->pbase() == nullptr)
			return traits_type::eof();

		auto &zstream = *m_zstream;

		zstream.next_in = reinterpret_cast<unsigned char *>(this->pbase());
		zstream.avail_in = static_cast<uInt>(this->pptr() - this->pbase());

		char_type buffer[BufferSize];

		bool failed = false;

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
				{
					this->set_error();
					return traits_type::eof();
				}
			}

			if (err == Z_OK)
			{
				if (zstream.avail_out == 0)
					continue;

				if (ch == traits_type::eof())
					continue;
			}
			else if (err != Z_STREAM_END)
			{
				this->set_error();
				failed = true;
			}

			break;
		}

		this->setp(this->m_in_buffer.data(), this->m_in_buffer.data() + this->m_in_buffer.size());

		if (failed or this->m_has_error)
			return traits_type::eof();

		if (not traits_type::eq_int_type(ch, traits_type::eof()))
		{
			*this->pptr() = traits_type::to_char_type(ch);
			this->pbump(1);
		}

		return ch;
	}

	/// \brief Flush the pending input through zlib and out to the upstream
	///
	/// Unlike close(), this does not finish the compressed stream; writing can
	/// continue afterwards.
	int sync() override
	{
		if (this->m_has_error)
			return -1;

		if (not m_zstream)
			return 0;

		if (this->pbase() == nullptr or this->pptr() == nullptr)
			return -1;

		auto &zstream = *m_zstream;

		zstream.next_in = reinterpret_cast<unsigned char *>(this->pbase());
		zstream.avail_in = static_cast<uInt>(this->pptr() - this->pbase());

		char_type buffer[BufferSize];

		for (;;)
		{
			zstream.next_out = reinterpret_cast<unsigned char *>(buffer);
			zstream.avail_out = sizeof(buffer);

			int err = ::deflate(&zstream, Z_SYNC_FLUSH);

			std::streamsize n = sizeof(buffer) - zstream.avail_out;
			if (n > 0 and this->m_upstream->sputn(reinterpret_cast<char_type *>(buffer), n) != n)
			{
				this->set_error();
				return -1;
			}

			if (err != Z_OK)
			{
				this->set_error();
				return -1;
			}

			if (zstream.avail_out != 0 and zstream.avail_in == 0)
				break;
		}

		this->setp(this->m_in_buffer.data(), this->m_in_buffer.data() + this->m_in_buffer.size());

		if (this->m_upstream->pubsync() != 0)
		{
			this->set_error();
			return -1;
		}

		return 0;
	}

  private:
	/// \brief The zlib internal structures are maintained as pointers to avoid having
	/// to copy their content in move constructors.
	std::unique_ptr<z_stream_s> m_zstream;

	/// \brief The zlib internal structures are maintained as pointers to avoid having
	/// to copy their content in move constructors.
	std::unique_ptr<gz_header> m_gzheader;

	/// \brief The compression level used, 0 (none) to 9 (maximum)
	int m_compression_level = kDefaultCompressionLevel;


	/// \brief Input buffer, this is the input for zlib
	std::array<char_type, BufferSize> m_in_buffer;
};

// --------------------------------------------------------------------
#if HAVE_LibLZMA

/// \brief A streambuf class that can be used to decompress xz data
///
/// \tparam CharT		Type of the character stream.
/// \tparam Traits		Traits for character type, defaults to char_traits<_CharT>.
/// \tparam BufferSize	The size of the internal buffers.
///
/// This implementation of streambuf can decompress data compressed
/// using xz or lzma.

template <typename CharT, typename Traits, size_t BufferSize = kDefaultBufferSize, std::enable_if_t<sizeof(CharT) == 1, int> = 0>
class basic_ixz_streambuf : public basic_streambuf<CharT, Traits>
{
  public:
	using char_type = CharT;
	using traits_type = Traits;

	using streambuf_type = std::basic_streambuf<char_type, traits_type>;
	using base_type = basic_streambuf<CharT, Traits>;

	using int_type = typename traits_type::int_type;
	using pos_type = typename traits_type::pos_type;
	using off_type = typename traits_type::off_type;

	basic_ixz_streambuf() = default;

	basic_ixz_streambuf(const basic_ixz_streambuf &) = delete;

	/// \brief Move constructor
	basic_ixz_streambuf(basic_ixz_streambuf &&rhs) noexcept
		: base_type(static_cast<base_type &&>(rhs))
	{
		std::swap(m_xzstream, rhs.m_xzstream);
		m_finished = rhs.m_finished;

		auto p = m_out_buffer.data();
		if (rhs.gptr() != nullptr and rhs.egptr() != nullptr)
			p = std::copy(rhs.gptr(), rhs.egptr(), m_out_buffer.data());
		this->setg(m_out_buffer.data(), m_out_buffer.data(), p);

		if (m_xzstream and m_xzstream->avail_in > 0)
		{
			auto next_in_offset = m_xzstream->next_in - reinterpret_cast<const unsigned char *>(rhs.m_in_buffer.data());
			std::copy(rhs.m_in_buffer.data() + next_in_offset,
				rhs.m_in_buffer.data() + next_in_offset + m_xzstream->avail_in,
				m_in_buffer.data());
			m_xzstream->next_in = reinterpret_cast<const unsigned char *>(m_in_buffer.data());
		}
	}

	basic_ixz_streambuf &operator=(const basic_ixz_streambuf &) = delete;

	/// \brief Move operator= implementation
	basic_ixz_streambuf &operator=(basic_ixz_streambuf &&rhs) noexcept
	{
		base_type::operator=(static_cast<base_type &&>(rhs));
		std::swap(m_xzstream, rhs.m_xzstream);
		m_finished = rhs.m_finished;

		auto p = m_out_buffer.data();
		if (rhs.gptr() != nullptr and rhs.egptr() != nullptr)
			p = std::copy(rhs.gptr(), rhs.egptr(), m_out_buffer.data());
		this->setg(m_out_buffer.data(), m_out_buffer.data(), p);

		if (m_xzstream and m_xzstream->avail_in > 0)
		{
			auto next_in_offset = m_xzstream->next_in - reinterpret_cast<unsigned char *>(rhs.m_in_buffer.data());
			std::copy(rhs.m_in_buffer.data() + next_in_offset,
				rhs.m_in_buffer.data() + next_in_offset + m_xzstream->avail_in,
				m_in_buffer.data());
			m_xzstream->next_in = reinterpret_cast<unsigned char *>(m_in_buffer.data());
		}

		return *this;
	}

	~basic_ixz_streambuf() override
	{
		close();
	}

	/// \brief This closes the xz stream and sets the get pointers to null.
	base_type *close() noexcept override
	{
		if (m_xzstream)
		{
			::lzma_end(m_xzstream.get());

			m_xzstream.reset(nullptr);
		}

		m_finished = false;
		this->setg(nullptr, nullptr, nullptr);

		return this;
	}

	/// \brief Initialize a xz stream and set the upstream.
	///
	/// \param upstream The upstream streambuf
	///
	/// The zstream is constructed and an optional header is
	/// read from upstream. The contents of the header are ignored
	/// but we must maintain that structure.
	base_type *init(streambuf_type *upstream) override
	{
		if (upstream == nullptr)
			return nullptr;

		this->set_upstream(upstream);

		close();

		m_finished = false;
		this->m_has_error = false;

		m_xzstream = std::make_unique<lzma_stream>();

		auto &xzstream = *m_xzstream.get();
		xzstream = LZMA_STREAM_INIT;

		int err = lzma_stream_decoder(&xzstream, UINT64_MAX, LZMA_TELL_NO_CHECK | LZMA_CONCATENATED);

		return err == LZMA_OK ? this : nullptr;
	}

  private:
	/// \brief The actual work is done here.
	int_type underflow() override
	{
		if (this->m_has_error)
			return traits_type::eof();

		if (m_xzstream and this->m_upstream)
		{
			auto &zstream = *m_xzstream.get();
			const std::streamsize kBufferByteSize = m_out_buffer.size();

			while (not m_finished and this->gptr() == this->egptr())
			{
				zstream.next_out = reinterpret_cast<unsigned char *>(m_out_buffer.data());
				zstream.avail_out = kBufferByteSize;

				bool at_eof = false;
				if (zstream.avail_in == 0)
				{
					zstream.next_in = reinterpret_cast<unsigned char *>(m_in_buffer.data());
					zstream.avail_in = static_cast<size_t>(detail::sgetn(this->m_upstream, m_in_buffer));
					at_eof = zstream.avail_in == 0;
				}

				// Because LZMA_CONCATENATED is used, the decoder must be told
				// that the input has ended. This is also how truncation is
				// detected.
				const lzma_action action = at_eof ? LZMA_FINISH : LZMA_RUN;
				int err = ::lzma_code(&zstream, action);
				std::streamsize n = kBufferByteSize - zstream.avail_out;

				// LZMA_NO_CHECK is an informational return code that is only
				// produced because we enabled LZMA_TELL_NO_CHECK. It means the
				// stream has no integrity check; decoding must simply continue.
				if (err == LZMA_NO_CHECK)
					err = LZMA_OK;

				if (err == LZMA_STREAM_END)
				{
					m_finished = true;

					if (n > 0)
					{
						this->setg(
							m_out_buffer.data(),
							m_out_buffer.data(),
							m_out_buffer.data() + n);
					}

					break;
				}

				if (err != LZMA_OK)
				{
					this->set_error();
					break;
				}

				if (n > 0)
				{
					this->setg(
						m_out_buffer.data(),
						m_out_buffer.data(),
						m_out_buffer.data() + n);
					break;
				}
			}
		}

		return this->gptr() != this->egptr() ? traits_type::to_int_type(*this->gptr()) : traits_type::eof();
	}

  private:
	/// \brief The lzma internal structures are maintained as pointers to avoid having
	/// to copy their content in move constructors.
	std::unique_ptr<lzma_stream> m_xzstream;

	/// \brief Set when the end of the xz stream was reached
	bool m_finished = false;

	/// \brief Input buffer, this is the input for xz
	std::array<char_type, BufferSize> m_in_buffer;

	/// \brief Output buffer, where the ostream finds the data
	std::array<char_type, BufferSize> m_out_buffer;
};

// --------------------------------------------------------------------

/// \brief A streambuf class that can be used to compress data using xz
///
/// \tparam CharT		Type of the character stream.
/// \tparam Traits		Traits for character type, defaults to char_traits<_CharT>.
/// \tparam BufferSize	The size of the internal buffers.
///
/// This implementation of streambuf can compress (deflate) data using xz.

template <typename CharT, typename Traits, size_t BufferSize = kDefaultBufferSize, std::enable_if_t<sizeof(CharT) == 1, int> = 0>
class basic_oxz_streambuf : public basic_streambuf<CharT, Traits>
{
  public:
	using char_type = CharT;
	using traits_type = Traits;

	using streambuf_type = std::basic_streambuf<char_type, traits_type>;
	using base_type = basic_streambuf<CharT, Traits>;

	using int_type = typename traits_type::int_type;
	using pos_type = typename traits_type::pos_type;
	using off_type = typename traits_type::off_type;

	basic_oxz_streambuf() = default;

	basic_oxz_streambuf(const basic_oxz_streambuf &) = delete;

	/// \brief Move constructor
	basic_oxz_streambuf(basic_oxz_streambuf &&rhs) noexcept
		: base_type(static_cast<base_type &&>(rhs))
	{
		std::swap(m_xzstream, rhs.m_xzstream);

		m_compression_level = rhs.m_compression_level;

		this->setp(m_in_buffer.data(), m_in_buffer.data() + m_in_buffer.size());
		if (rhs.pbase() != nullptr and rhs.pptr() != nullptr)
			this->sputn(rhs.pbase(), rhs.pptr() - rhs.pbase());
		rhs.setp(nullptr, nullptr);
	}

	basic_oxz_streambuf &operator=(const basic_oxz_streambuf &) = delete;

	/// \brief Move operator=
	basic_oxz_streambuf &operator=(basic_oxz_streambuf &&rhs) noexcept
	{
		if (this == &rhs)
			return *this;

		// Flush the old state before the base move discards our upstream.
		close();

		base_type::operator=(static_cast<base_type &&>(rhs));

		std::swap(m_xzstream, rhs.m_xzstream);

		m_compression_level = rhs.m_compression_level;

		this->setp(m_in_buffer.data(), m_in_buffer.data() + m_in_buffer.size());
		if (rhs.pbase() != nullptr and rhs.pptr() != nullptr)
			this->sputn(rhs.pbase(), rhs.pptr() - rhs.pbase());
		rhs.setp(nullptr, nullptr);

		return *this;
	}

	~basic_oxz_streambuf() override
	{
		close();
	}

	/// \brief This closes the xz stream and sets the put pointers to null.
	///
	/// Returns nullptr when flushing the remaining data failed, so that
	/// callers can propagate the error.
	base_type *close() override
	{
		if (m_xzstream)
		{
			overflow(traits_type::eof());

			::lzma_end(m_xzstream.get());

			m_xzstream.reset(nullptr);
		}

		this->setp(nullptr, nullptr);

		return this->m_has_error ? nullptr : this;
	}

	/// \brief Set the compression level, 0 (no compression) to 9 (maximum compression)
	void set_compression_level(int level) noexcept override
	{
		m_compression_level = std::clamp(level, 0, 9);
	}

	/// \brief Initialize the internal lzma structures
	///
	/// \param upstream The upstream streambuf
	///
	/// The zlib stream is initialized as one that can accept
	/// a gzip header.
	base_type *init(streambuf_type *upstream) override
	{
		if (upstream == nullptr)
			return nullptr;

		// Finish and tear down any existing stream while the old upstream is
		// still installed. Otherwise buffered data would be flushed to the new
		// upstream. Only then start a clean session on the new upstream.
		close();

		this->set_upstream(upstream);
		this->m_has_error = false;

		m_xzstream = std::make_unique<lzma_stream>();

		auto &zstream = *m_xzstream.get();
		zstream = LZMA_STREAM_INIT;

		int err = lzma_easy_encoder(&zstream, m_compression_level, LZMA_CHECK_CRC64);

		if (err == LZMA_OK)
			this->setp(this->m_in_buffer.data(), this->m_in_buffer.data() + this->m_in_buffer.size());

		return err == LZMA_OK ? this : nullptr;
	}

  private:
	/// \brief The actual work is done here
	///
	/// \param ch The character that did not fit, in case it is eof we need to flush
	///
	int_type overflow(int_type ch) override
	{
		if (not m_xzstream or this->m_has_error or this->pbase() == nullptr)
			return traits_type::eof();

		auto &zstream = *m_xzstream;

		zstream.next_in = reinterpret_cast<unsigned char *>(this->pbase());
		zstream.avail_in = this->pptr() - this->pbase();

		char_type buffer[BufferSize];

		bool failed = false;

		for (;;)
		{
			zstream.next_out = reinterpret_cast<unsigned char *>(buffer);
			zstream.avail_out = sizeof(buffer);

			int err = ::lzma_code(&zstream, ch == traits_type::eof() ? LZMA_FINISH : LZMA_RUN);

			std::streamsize n = sizeof(buffer) - zstream.avail_out;
			if (n > 0)
			{
				auto r = this->m_upstream->sputn(reinterpret_cast<char_type *>(buffer), n);

				if (r != n)
				{
					this->set_error();
					return traits_type::eof();
				}
			}

			if (err == LZMA_OK)
			{
				if (zstream.avail_out == 0)
					continue;

				if (ch == traits_type::eof())
					continue;
			}
			else if (err != LZMA_STREAM_END)
			{
				this->set_error();
				failed = true;
			}

			break;
		}

		this->setp(this->m_in_buffer.data(), this->m_in_buffer.data() + this->m_in_buffer.size());

		if (failed or this->m_has_error)
			return traits_type::eof();

		if (not traits_type::eq_int_type(ch, traits_type::eof()))
		{
			*this->pptr() = traits_type::to_char_type(ch);
			this->pbump(1);
		}

		return ch;
	}

	/// \brief Flush the pending input through liblzma and out to the upstream
	///
	/// Unlike close(), this does not finish the compressed stream; writing can
	/// continue afterwards.
	int sync() override
	{
		if (this->m_has_error)
			return -1;

		if (not m_xzstream)
			return 0;

		if (this->pbase() == nullptr or this->pptr() == nullptr)
			return -1;

		auto &zstream = *m_xzstream;

		zstream.next_in = reinterpret_cast<unsigned char *>(this->pbase());
		zstream.avail_in = this->pptr() - this->pbase();

		char_type buffer[BufferSize];

		for (;;)
		{
			zstream.next_out = reinterpret_cast<unsigned char *>(buffer);
			zstream.avail_out = sizeof(buffer);

			int err = ::lzma_code(&zstream, LZMA_SYNC_FLUSH);

			std::streamsize n = sizeof(buffer) - zstream.avail_out;
			if (n > 0 and this->m_upstream->sputn(reinterpret_cast<char_type *>(buffer), n) != n)
			{
				this->set_error();
				return -1;
			}

			// For an encoder, LZMA_STREAM_END after LZMA_SYNC_FLUSH means the
			// flush request completed; the stream itself is still alive.
			if (err == LZMA_STREAM_END)
				break;

			if (err != LZMA_OK)
			{
				this->set_error();
				return -1;
			}
		}

		this->setp(this->m_in_buffer.data(), this->m_in_buffer.data() + this->m_in_buffer.size());

		if (this->m_upstream->pubsync() != 0)
		{
			this->set_error();
			return -1;
		}

		return 0;
	}

  private:
	/// \brief The lzma internal structures are maintained as pointers to avoid having
	/// to copy their content in move constructors.
	std::unique_ptr<lzma_stream> m_xzstream;

	/// \brief The compression level used, 0 (none) to 9 (maximum)
	int m_compression_level = kDefaultCompressionLevel;


	/// \brief Input buffer, this is the input for lzma
	std::array<char_type, BufferSize> m_in_buffer;
};

#endif

// --------------------------------------------------------------------

/// \brief A streambuf that replays a number of bytes before delegating to another streambuf
///
/// \tparam CharT		Type of the character stream.
/// \tparam Traits		Traits for character type, defaults to char_traits<_CharT>.
///
/// This streambuf is used while sniffing the format of an input stream. The signature
/// bytes that were already read from the stream are played back from an internal buffer
/// before the actual upstream streambuf is consulted. This avoids the need to put back
/// bytes into the upstream streambuf, which not every streambuf implementation supports.

template <typename CharT, typename Traits>
class basic_prefix_streambuf : public std::basic_streambuf<CharT, Traits>
{
  public:
	using char_type = CharT;
	using traits_type = Traits;

	using int_type = typename traits_type::int_type;
	using streambuf_type = std::basic_streambuf<CharT, Traits>;

	basic_prefix_streambuf(streambuf_type *upstream, const char_type *prefix, std::streamsize length) noexcept
		: m_upstream(upstream)
	{
		const auto max_length = static_cast<std::streamsize>(m_prefix.size());
		const auto count = std::clamp(length, std::streamsize{ 0 }, max_length);
		std::copy(prefix, prefix + count, m_prefix.begin());
		this->setg(m_prefix.data(), m_prefix.data(), m_prefix.data() + count);
	}

  protected:
	int_type underflow() override
	{
		if (this->gptr() < this->egptr())
			return traits_type::to_int_type(*this->gptr());

		if (m_upstream == nullptr)
			return traits_type::eof();

		std::streamsize n = detail::sgetn(m_upstream, m_buffer);

		if (n == 0)
			return traits_type::eof();

		this->setg(m_buffer.data(), m_buffer.data(), m_buffer.data() + n);

		return traits_type::to_int_type(*this->gptr());
	}

  private:
	/// \brief The streambuf we read from once the prefix is exhausted
	streambuf_type *m_upstream;

	/// \brief The signature bytes that were already consumed from the upstream
	std::array<char_type, kMaxSignatureLength> m_prefix;

	/// \brief A buffer to read the remaining data into
	std::array<char_type, kDefaultBufferSize> m_buffer;
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
	using base_type = std::basic_istream<CharT, Traits>;

	using traits_type = Traits;
	using char_type = typename traits_type::char_type;
	using int_type = typename traits_type::int_type;

	using z_streambuf_type = basic_streambuf<char_type, traits_type>;
	using upstreambuf_type = std::basic_streambuf<char_type, traits_type>;

	using gzip_streambuf_type = basic_igzip_streambuf<char_type, traits_type>;
#if HAVE_LibLZMA
	using xz_streambuf_type = basic_ixz_streambuf<char_type, traits_type>;
#endif
	using prefix_streambuf_type = basic_prefix_streambuf<char_type, traits_type>;

	/// \brief Regular move constructor
	basic_istream(basic_istream &&rhs) noexcept
		: base_type(static_cast<base_type &&>(rhs))
	{
		m_gxriobuf = std::move(rhs.m_gxriobuf);
		m_upstream_buf = std::move(rhs.m_upstream_buf);
		m_upstream = rhs.m_upstream;

		if (m_gxriobuf)
			m_gxriobuf->set_ios(this);

		this->rdbuf(m_gxriobuf ? m_gxriobuf.get() : m_upstream);
	}

	/// \brief Regular move operator=
	basic_istream &operator=(basic_istream &&rhs) noexcept
	{
		if (this == &rhs)
			return *this;

		base_type::operator=(static_cast<base_type &&>(rhs));
		m_gxriobuf = std::move(rhs.m_gxriobuf);
		m_upstream_buf = std::move(rhs.m_upstream_buf);
		m_upstream = rhs.m_upstream;

		if (m_gxriobuf)
			m_gxriobuf->set_ios(this);

		this->rdbuf(m_gxriobuf ? m_gxriobuf.get() : m_upstream);

		return *this;
	}

	/// \brief Return true if a (de)compression error was encountered
	///
	/// This reports errors from the underlying compression streambuf. It is
	/// false for streams that read uncompressed data.
	[[nodiscard]] bool error() const noexcept
	{
		return m_gxriobuf != nullptr and m_gxriobuf->error();
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
	///
	/// The signature bytes that are consumed while sniffing are replayed to
	/// the reader through a basic_prefix_streambuf, so the upstream does not
	/// need to support putting characters back.

	void init_z(upstreambuf_type *sb)
	{
		if (sb == nullptr)
		{
			this->clear();
			this->setstate(std::ios_base::failbit);
			return;
		}

		m_gxriobuf.reset(nullptr);
		m_upstream_buf.reset(nullptr);
		m_upstream = nullptr;

		char_type signature[kMaxSignatureLength];
		std::streamsize signature_length = 0;

		int_type ch = sb->sgetc();
		if (ch == 0x1f)
		{
			signature[signature_length++] = traits_type::to_char_type(ch);
			sb->sbumpc();
			ch = sb->sgetc();

			if (ch == 0x8b) // gzip magic
			{
				signature[signature_length++] = traits_type::to_char_type(ch);
				sb->sbumpc();
				ch = sb->sgetc();

				if (ch == 0x08) // compression method must be deflate
				{
					signature[signature_length++] = traits_type::to_char_type(ch);
					sb->sbumpc();
					m_gxriobuf.reset(new gzip_streambuf_type);
				}
			}
		}
		else if (ch == 0xfd)
		{
			// xz magic: 0xfd '7' 'z' 'X' 'Z' 0x00
			constexpr std::streamsize kXZMagicLength = 6;
			constexpr int kXZMagic[kXZMagicLength] = { 0xfd, 0x37, 0x7a, 0x58, 0x5a, 0x00 };

			int_type next = ch;
			for (std::streamsize i = 0; i < kXZMagicLength and next == kXZMagic[i]; ++i)
			{
				signature[signature_length++] = traits_type::to_char_type(next);
				sb->sbumpc();
				next = sb->sgetc();
			}

			if (signature_length == kXZMagicLength)
			{
#if HAVE_LibLZMA
				m_gxriobuf.reset(new xz_streambuf_type);
#else
				// xz data cannot be decoded in this build. Detaching rdbuf sets
				// badbit instead of returning the compressed bytes, and clear()
				// cannot make a stream without an rdbuf readable again.
				this->rdbuf(nullptr);
				return;
#endif
			}
		}

		if (m_gxriobuf)
		{
			m_upstream_buf.reset(new prefix_streambuf_type(sb, signature, signature_length));
			m_upstream = m_upstream_buf.get();

			if (not m_gxriobuf->init(m_upstream))
				this->setstate(std::ios_base::failbit);
			else
			{
				m_gxriobuf->set_ios(this);
				this->init(m_gxriobuf.get());
			}
		}
		else if (signature_length > 0)
		{
			m_upstream_buf.reset(new prefix_streambuf_type(sb, signature, signature_length));
			m_upstream = m_upstream_buf.get();
			this->init(m_upstream);
		}
		else
		{
			m_upstream = sb;
			this->init(m_upstream);
		}
	}

  protected:
	/// \brief Our streambuf class
	std::unique_ptr<z_streambuf_type> m_gxriobuf;

	/// \brief Replays the sniffed signature bytes, only used when a basic_istream
	/// is constructed directly with an upstream streambuf.
	std::unique_ptr<upstreambuf_type> m_upstream_buf;

	/// \brief The streambuf that provides the data to this stream, either the
	/// upstream directly or the prefix streambuf wrapping it.
	upstreambuf_type *m_upstream = nullptr;
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
	using base_type = basic_istream<CharT, Traits>;

	using char_type = CharT;
	using traits_type = Traits;

	using filebuf_type = std::basic_filebuf<char_type, traits_type>;

	using gzip_streambuf_type = typename base_type::gzip_streambuf_type;
#if HAVE_LibLZMA
	using xz_streambuf_type = typename base_type::xz_streambuf_type;
#endif

	/// \brief Default constructor, does not open a file since none is specified
	basic_ifstream()
	{
		this->rdbuf(&m_filebuf);
		this->m_upstream = &m_filebuf;
		this->clear();
	}

	~basic_ifstream() override
	{
		// Do not report close failures through setstate() here: that can throw
		// for a stream with exceptions enabled, and destructors must not throw.
		if (is_open())
			close_impl(false);
	}

	/// \brief Construct an ifstream
	/// \param filename Null terminated string specifying the file to open
	/// \param mode The mode in which to open the file

	explicit basic_ifstream(const char *filename, std::ios_base::openmode mode = std::ios_base::in)
		: basic_ifstream()
	{
		open(filename, mode);
	}

	/// \brief Construct an ifstream
	/// \param filename std::string specifying the file to open
	/// \param mode The mode in which to open the file

	explicit basic_ifstream(const std::string &filename, std::ios_base::openmode mode = std::ios_base::in)
		: basic_ifstream()
	{
		open(filename, mode);
	}

	/// \brief Construct an ifstream
	/// \param filename std::filesystem::path specifying the file to open
	/// \param mode The mode in which to open the file

	explicit basic_ifstream(const std::filesystem::path &filename, std::ios_base::openmode mode = std::ios_base::in)
		: basic_ifstream()
	{
		open(filename, mode);
	}

	/// \brief Move constructor
	basic_ifstream(basic_ifstream &&rhs) noexcept
		: base_type(static_cast<base_type &&>(rhs))
	{
		m_filebuf = std::move(rhs.m_filebuf);
		this->m_upstream = &m_filebuf;

		if (this->m_gxriobuf)
		{
			this->m_gxriobuf->set_upstream(&m_filebuf);
			this->m_gxriobuf->set_ios(this);
		}
		else
			this->rdbuf(&m_filebuf);
	}

	basic_ifstream(const basic_ifstream &) = delete;

	basic_ifstream &operator=(const basic_ifstream &) = delete;

	/// \brief Move version of operator=
	basic_ifstream &operator=(basic_ifstream &&rhs) noexcept
	{
		if (this == &rhs)
			return *this;

		base_type::operator=(static_cast<base_type &&>(rhs));

		m_filebuf = std::move(rhs.m_filebuf);
		this->m_upstream = &m_filebuf;

		if (this->m_gxriobuf)
		{
			this->m_gxriobuf->set_upstream(&m_filebuf);
			this->m_gxriobuf->set_ios(this);
		}
		else
			this->rdbuf(&m_filebuf);

		return *this;
	}

	/// \brief Open the file \a filename with mode \a mode
	/// \param filename std::filesystem::path specifying the file to open
	/// \param mode The mode in which to open the file

	void open(const std::filesystem::path &filename, std::ios_base::openmode mode = std::ios_base::in)
	{
#if !HAVE_LibLZMA
		if (filename.extension() == ".xz")
		{
			// Reject unsupported input before touching the file. Detaching
			// rdbuf sets badbit and stops clear() from making the stream return
			// the still-compressed bytes.
			this->m_gxriobuf.reset(nullptr);
			this->rdbuf(nullptr);
			this->m_upstream = nullptr;
			return;
		}
#endif

		if (not m_filebuf.open(filename, mode | std::ios::binary))
			this->setstate(std::ios_base::failbit);
		else
		{
			if (filename.extension() == ".gz")
				this->m_gxriobuf.reset(new gzip_streambuf_type);
#if HAVE_LibLZMA
			else if (filename.extension() == ".xz")
				this->m_gxriobuf.reset(new xz_streambuf_type);
#endif
			else
				this->m_gxriobuf.reset(nullptr);

			this->m_upstream = &m_filebuf;

			if (not this->m_gxriobuf)
			{
				this->rdbuf(&m_filebuf);
				this->clear();
			}
			else if (not this->m_gxriobuf->init(&m_filebuf))
			{
				// Do not leave rdbuf() pointing at a compressor that is about
				// to be destroyed. Close the file and fall back to the filebuf.
				this->rdbuf(&m_filebuf);
				this->m_upstream = &m_filebuf;
				m_filebuf.close();
				this->m_gxriobuf.reset(nullptr);
				this->setstate(std::ios_base::failbit);
			}
			else
			{
				this->m_gxriobuf->set_ios(this);
				this->rdbuf(this->m_gxriobuf.get());
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

	[[nodiscard]] bool is_open() const noexcept
	{
		return m_filebuf.is_open();
	}

	/// \brief Close the file
	///
	/// Calls m_filebuf.close(). If that fails, the failbit is set.

	void close()
	{
		close_impl(true);
	}

	/// \brief Swap the contents with those of \a rhs
	/// \param rhs The ifstream to swap with

	void swap(basic_ifstream &rhs) noexcept
	{
		base_type::swap(rhs);
		m_filebuf.swap(rhs.m_filebuf);
		std::swap(this->m_gxriobuf, rhs.m_gxriobuf);

		if (this->m_gxriobuf)
		{
			this->m_gxriobuf->set_upstream(&m_filebuf);
			this->m_gxriobuf->set_ios(this);
			this->rdbuf(this->m_gxriobuf.get());
		}
		else
			this->rdbuf(&m_filebuf);

		this->m_upstream = &m_filebuf;

		if (rhs.m_gxriobuf)
		{
			rhs.m_gxriobuf->set_upstream(&rhs.m_filebuf);
			rhs.m_gxriobuf->set_ios(&rhs);
			rhs.rdbuf(rhs.m_gxriobuf.get());
		}
		else
			rhs.rdbuf(&rhs.m_filebuf);

		rhs.m_upstream = &rhs.m_filebuf;
	}

  private:
	/// \brief Close the compressor and the file
	///
	/// When \a report_error is false, failures are not reported through
	/// setstate(), which keeps the destructor path from throwing.
	void close_impl(bool report_error)
	{
		bool failed = false;

		if (this->m_gxriobuf and not this->m_gxriobuf->close())
			failed = true;

		if (not m_filebuf.close())
			failed = true;

		if (report_error and failed)
			this->setstate(std::ios_base::failbit);
	}

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
	using base_type = std::basic_ostream<CharT, Traits>;

	using char_type = CharT;
	using traits_type = Traits;

	using z_streambuf_type = basic_streambuf<char_type, traits_type>;
	using upstreambuf_type = std::basic_streambuf<char_type, traits_type>;

	/// \brief Construct an ostream that writes uncompressed data to \a buf
	///
	/// No compression is applied; use gxrio::ofstream or a compressing
	/// streambuf when compressed output is required.
	explicit basic_ostream(upstreambuf_type *buf)
		: base_type(buf), m_upstream(buf)
	{
		if (buf == nullptr)
		{
			this->clear();
			this->setstate(std::ios_base::failbit);
		}
	}

	/// \brief Regular move constructor
	basic_ostream(basic_ostream &&rhs) noexcept
		: base_type(static_cast<base_type &&>(rhs))
	{
		m_gxriobuf = std::move(rhs.m_gxriobuf);
		m_upstream = rhs.m_upstream;

		this->rdbuf(m_gxriobuf ? m_gxriobuf.get() : m_upstream);
	}

	/// \brief Regular move operator=
	basic_ostream &operator=(basic_ostream &&rhs) noexcept
	{
		if (this == &rhs)
			return *this;

		base_type::operator=(static_cast<base_type &&>(rhs));
		m_gxriobuf = std::move(rhs.m_gxriobuf);
		m_upstream = rhs.m_upstream;

		this->rdbuf(m_gxriobuf ? m_gxriobuf.get() : m_upstream);

		return *this;
	}

	/// \brief Return true if a (de)compression error was encountered
	///
	/// This reports errors from the underlying compression streambuf. It is
	/// false for streams that write uncompressed data.
	[[nodiscard]] bool error() const noexcept
	{
		return m_gxriobuf != nullptr and m_gxriobuf->error();
	}

  protected:
	basic_ostream()
		: base_type(nullptr) {}

  protected:
	/// \brief Our streambuf class
	std::unique_ptr<z_streambuf_type> m_gxriobuf;

	/// \brief The streambuf that receives the output, the upstream of the compressing streambuf.
	upstreambuf_type *m_upstream = nullptr;
};

// --------------------------------------------------------------------

/// \brief Control output to files compressing the contents with gzip.
///
/// \tparam CharT		Type of the character stream.
/// \tparam Traits		Traits for character type, defaults to char_traits<_CharT>.
///
/// This is an ofstream implementation that can write to named files compressing the content
/// with gzip directly. The class inherits from std::basic_ostream and offers all the
/// associated functionality.

template <typename CharT, typename Traits>
class basic_ofstream : public basic_ostream<CharT, Traits>
{
  public:
	using base_type = basic_ostream<CharT, Traits>;

	using char_type = CharT;
	using traits_type = Traits;

	using filebuf_type = std::basic_filebuf<char_type, traits_type>;
	using gzip_streambuf_type = basic_ogzip_streambuf<char_type, traits_type>;
#if HAVE_LibLZMA
	using xz_streambuf_type = basic_oxz_streambuf<char_type, traits_type>;
#endif

	/// \brief Default constructor, does not open a file since none is specified
	basic_ofstream()
	{
		this->rdbuf(&m_filebuf);
		this->m_upstream = &m_filebuf;
		this->clear();
	}

	~basic_ofstream() override
	{
		// Do not report close failures through setstate() here: that can throw
		// for a stream with exceptions enabled, and destructors must not throw.
		if (is_open())
			close_impl(false);
	}

	/// \brief Construct an ofstream
	/// \param filename Null terminated string specifying the file to open
	/// \param mode The mode in which to open the file

	explicit basic_ofstream(const char *filename, std::ios_base::openmode mode = std::ios_base::out)
		: basic_ofstream()
	{
		open(filename, mode);
	}

	/// \brief Construct an ofstream
	/// \param filename std::string specifying the file to open
	/// \param mode The mode in which to open the file

	explicit basic_ofstream(const std::string &filename, std::ios_base::openmode mode = std::ios_base::out)
		: basic_ofstream()
	{
		open(filename, mode);
	}

	/// \brief Construct an ofstream
	/// \param filename std::filesystem::path specifying the file to open
	/// \param mode The mode in which to open the file

	explicit basic_ofstream(const std::filesystem::path &filename, std::ios_base::openmode mode = std::ios_base::out)
		: basic_ofstream()
	{
		open(filename, mode);
	}

	/// \brief Move constructor
	basic_ofstream(basic_ofstream &&rhs) noexcept
		: base_type(static_cast<base_type &&>(rhs))
	{
		m_filebuf = std::move(rhs.m_filebuf);
		this->m_upstream = &m_filebuf;
		m_compression_level = rhs.m_compression_level;

		if (this->m_gxriobuf)
		{
			this->m_gxriobuf->set_upstream(&m_filebuf);
			this->m_gxriobuf->set_ios(this);
		}
		else
			this->rdbuf(&m_filebuf);
	}

	basic_ofstream(const basic_ofstream &) = delete;

	basic_ofstream &operator=(const basic_ofstream &) = delete;

	/// \brief Move operator=
	basic_ofstream &operator=(basic_ofstream &&rhs) noexcept
	{
		if (this == &rhs)
			return *this;

		base_type::operator=(static_cast<base_type &&>(rhs));
		m_filebuf = std::move(rhs.m_filebuf);
		this->m_upstream = &m_filebuf;
		m_compression_level = rhs.m_compression_level;

		if (this->m_gxriobuf)
		{
			this->m_gxriobuf->set_upstream(&m_filebuf);
			this->m_gxriobuf->set_ios(this);
		}
		else
			this->rdbuf(&m_filebuf);

		return *this;
	}

	/// \brief Set the compression level to use, 0 (no compression) to 9 (maximum).
	///
	/// This applies to files that are opened with a .gz or .xz extension.
	/// The default is 9. When called on an already open compressed file the
	/// level is forwarded to the active compressor; for gzip it takes effect
	/// immediately, for xz it applies the next time a file is opened.

	void set_compression_level(int level)
	{
		m_compression_level = std::clamp(level, 0, 9);

		if (this->m_gxriobuf)
			this->m_gxriobuf->set_compression_level(m_compression_level);
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
#if !HAVE_LibLZMA
		if (filename.extension() == ".xz")
		{
			// Reject unsupported output before creating or truncating the
			// destination file. Detaching rdbuf sets badbit and stops clear()
			// from making the stream writable again.
			this->m_gxriobuf.reset(nullptr);
			this->rdbuf(nullptr);
			this->m_upstream = nullptr;
			return;
		}
#endif

		if (not m_filebuf.open(filename, mode | std::ios::binary))
			this->setstate(std::ios_base::failbit);
		else
		{
			if (filename.extension() == ".gz")
				this->m_gxriobuf.reset(new gzip_streambuf_type);
#if HAVE_LibLZMA
			else if (filename.extension() == ".xz")
				this->m_gxriobuf.reset(new xz_streambuf_type);
#endif
			else
				this->m_gxriobuf.reset(nullptr);

			this->m_upstream = &m_filebuf;

			if (this->m_gxriobuf)
			{
				this->m_gxriobuf->set_compression_level(m_compression_level);

				if (not this->m_gxriobuf->init(&m_filebuf))
				{
					// Do not leave rdbuf() pointing at a compressor that is
					// about to be destroyed. Close the file and fall back to
					// the filebuf.
					this->rdbuf(&m_filebuf);
					this->m_upstream = &m_filebuf;
					m_filebuf.close();
					this->m_gxriobuf.reset(nullptr);
					this->setstate(std::ios_base::failbit);
				}
				else
				{
					this->m_gxriobuf->set_ios(this);
					this->rdbuf(this->m_gxriobuf.get());
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

	[[nodiscard]] bool is_open() const noexcept
	{
		return m_filebuf.is_open();
	}

	/// \brief Close the file
	///
	/// Calls m_filebuf.close(). If that fails, the failbit is set.

	void close()
	{
		bool failed = false;

		if (this->m_gxriobuf and not this->m_gxriobuf->close())
			failed = true;

		if (not m_filebuf.close())
			failed = true;

		// Report last: setstate() may throw when the stream has exceptions
		// enabled, and by then both the compressor and the file must be closed.
		if (failed)
			this->setstate(std::ios_base::failbit);
	}

	/// \brief Swap the contents with those of \a rhs
	/// \param rhs The ifstream to swap with

	void swap(basic_ofstream &rhs) noexcept
	{
		base_type::swap(rhs);
		m_filebuf.swap(rhs.m_filebuf);
		std::swap(this->m_gxriobuf, rhs.m_gxriobuf);
		std::swap(m_compression_level, rhs.m_compression_level);

		if (this->m_gxriobuf)
		{
			this->m_gxriobuf->set_upstream(&m_filebuf);
			this->m_gxriobuf->set_ios(this);
			this->rdbuf(this->m_gxriobuf.get());
		}
		else
			this->rdbuf(&m_filebuf);

		this->m_upstream = &m_filebuf;

		if (rhs.m_gxriobuf)
		{
			rhs.m_gxriobuf->set_upstream(&rhs.m_filebuf);
			rhs.m_gxriobuf->set_ios(&rhs);
			rhs.rdbuf(rhs.m_gxriobuf.get());
		}
		else
			rhs.rdbuf(&rhs.m_filebuf);

		rhs.m_upstream = &rhs.m_filebuf;
	}

  private:
	/// \brief Close the compressor and the file
	///
	/// When \a report_error is false, failures are not reported through
	/// setstate(), which keeps the destructor path from throwing.
	void close_impl(bool report_error)
	{
		bool failed = false;

		if (this->m_gxriobuf and not this->m_gxriobuf->close())
			failed = true;

		if (not m_filebuf.close())
			failed = true;

		if (report_error and failed)
			this->setstate(std::ios_base::failbit);
	}

	/// \brief The filebuf
	filebuf_type m_filebuf;

	/// \brief The compression level used when opening compressed files
	int m_compression_level = kDefaultCompressionLevel;
};

// --------------------------------------------------------------------

/// \brief Convenience typedefs
using istream = basic_istream<char, std::char_traits<char>>;
using ostream = basic_ostream<char, std::char_traits<char>>;
using ifstream = basic_ifstream<char, std::char_traits<char>>;
using ofstream = basic_ofstream<char, std::char_traits<char>>;

} // namespace gxrio
