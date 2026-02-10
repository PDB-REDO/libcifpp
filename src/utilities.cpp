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

#include "cif++/utilities.hpp"

#include "revision.hpp"

#include <atomic>
#include <cassert>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstddef>
#include <cstdlib>
#include <deque>
#include <exception>
#include <format>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <stop_token>
#include <string>
#include <system_error>
#include <thread>
#include <utility>

#if __cpp_lib_jthread >= 201911L
# include <stop_token>
#endif

namespace fs = std::filesystem;

// --------------------------------------------------------------------

namespace cif
{

int VERBOSE = 0;

// --------------------------------------------------------------------

std::string get_version_nr()
{
	std::ostringstream s;
	write_version_string(s, false);
	return s.str();
}

// --------------------------------------------------------------------

#if defined(_WIN32) or defined(__MINGW32__)
}
// clang-format off
# include <windows.h>
# include <libloaderapi.h>
# include <wincon.h>
// clang-format on

namespace cif
{

uint32_t get_terminal_width()
{
	CONSOLE_SCREEN_BUFFER_INFO csbi;
	if (::GetConsoleScreenBufferInfo(::GetStdHandle(STD_OUTPUT_HANDLE), &csbi))
		return csbi.srWindow.Right - csbi.srWindow.Left + 1;
	return 80;
}

std::tuple<uint32_t, uint32_t> get_terminal_width_and_height()
{
	CONSOLE_SCREEN_BUFFER_INFO csbi;
	if (::GetConsoleScreenBufferInfo(::GetStdHandle(STD_OUTPUT_HANDLE), &csbi))
		return { csbi.srWindow.Right - csbi.srWindow.Left + 1, csbi.srWindow.Bottom - csbi.srWindow.Top + 1 };
	return { 80, 24 };
}

void write_to_console(const std::string &s)
{
	auto h = ::GetStdHandle(STD_OUTPUT_HANDLE);
	CONSOLE_SCREEN_BUFFER_INFO csbi;

	if (auto l = ::MultiByteToWideChar(CP_UTF8, 0, s.data(), s.length(), nullptr, 0);
		l > 0 and ::GetConsoleScreenBufferInfo(::GetStdHandle(STD_OUTPUT_HANDLE), &csbi))
	{
		std::u16string ws(l, 0);

		::MultiByteToWideChar(CP_UTF8, 0, s.data(), s.length(), (LPWSTR)ws.data(), l);

		DWORD w;
		::WriteConsoleW(h, ws.data(), ws.length(), &w, nullptr);
	}
	else
	{
		std::cout.write(s.data(), s.length());
		std::cout.flush();
	}
}

#else

# include <climits>
# include <sys/ioctl.h>
# include <termios.h>

uint32_t get_terminal_width()
{
	uint32_t width = 80;

	if (isatty(STDOUT_FILENO))
	{
		struct winsize w;
		ioctl(0, TIOCGWINSZ, &w);
		width = w.ws_col;
	}
	return width;
}

std::tuple<uint32_t, uint32_t> get_terminal_width_and_height()
{
	if (isatty(STDOUT_FILENO))
	{
		struct winsize w;
		ioctl(0, TIOCGWINSZ, &w);
		return { w.ws_col, w.ws_row };
	}
	return { 80, 24 };
}

inline void write_to_console(const std::string &s)
{
	std::cout << s << std::flush;
}

#endif

// --------------------------------------------------------------------

struct progress_bar_impl
{
	progress_bar_impl(uint64_t max_value, const std::string &message)
		: m_max_value(max_value)
		, m_consumed(0)
		, m_action(message)
		, m_message(message)
	{
	}

	virtual ~progress_bar_impl() = default;

	virtual void consumed(uint64_t n);
	virtual void progress(uint64_t p);
	virtual void message(const std::string &msg);
	virtual void print_done();

	using time_point = std::chrono::time_point<std::chrono::system_clock>;

	uint64_t m_max_value;
	std::atomic<uint64_t> m_consumed;
	std::string m_action, m_message;
	time_point m_start = std::chrono::system_clock::now();
};

void progress_bar_impl::consumed(uint64_t n)
{
	m_consumed += n;
}

void progress_bar_impl::progress(uint64_t p)
{
	m_consumed = p;
}

void progress_bar_impl::message(const std::string &msg)
{
	m_message = msg;
}

void progress_bar_impl::print_done()
{
	std::chrono::duration<double> elapsed = std::chrono::system_clock::now() - m_start;
	std::string days, hours, minutes;

	auto s = static_cast<int>(std::rint(elapsed.count()));
	if (s > 24 * 60 * 60)
	{
		days = std::format("{:d}d ", s / (24 * 60 * 60));
		s %= 24 * 60 * 60;
	}

	if (s > 60 * 60)
	{
		hours = std::format("{:2d}h ", s / (60 * 60));
		s %= 60 * 60;
	}

	if (s > 60)
	{
		minutes = std::format("{:2d}m ", s / 60);
		s %= 60;
	}

	std::string msg = std::format("{} done in {}{}{}{:.1f}s", m_action, days, hours, minutes, s + 1e-6 * static_cast<double>(elapsed.count() - s));

	uint32_t width = get_terminal_width();

	if (msg.length() < width)
		msg += std::string(width - msg.length(), ' ');

	write_to_console(msg += '\n');
}

// --------------------------------------------------------------------

struct simple_progress_bar_impl : public progress_bar_impl
{
	simple_progress_bar_impl(uint64_t max_value, const std::string &message)
		: progress_bar_impl(max_value, message)
	{
	}

	~simple_progress_bar_impl() override
	{
		try
		{
			if (m_printed_any)
				print_done();
		}
		catch (const std::exception &ex)
		{
			std::cerr << "error finishing progress bar: " << ex.what() << '\n';
		}
	}

	void consumed(uint64_t n) override
	{
		using namespace std::literals;

		progress_bar_impl::consumed(n);

		// print at most 10 steps, but only if it took long enough

		int percentile = static_cast<int>(std::floor(10.f * m_consumed / m_max_value));
		if (percentile > m_last_percentile and (m_printed_any or std::chrono::system_clock::now() - m_start >= 1s))
		{
			if (not std::exchange(m_printed_any, true))
				write_to_console(m_action + ": ");

			write_to_console(std::format("...{:d}0%", percentile));
			m_last_percentile = percentile;
		}
	}

	void progress(uint64_t p) override
	{
		consumed(p - m_consumed);
	}

	void print_done() override
	{
		if (m_printed_any)
		{
			write_to_console("\n");
			progress_bar_impl::print_done();
		}
	}

	bool m_printed_any = false;
	int m_last_percentile = 0;
};

// --------------------------------------------------------------------

struct fancy_progress_bar_impl : public progress_bar_impl
{
	fancy_progress_bar_impl(uint64_t max_value, const std::string &message)
		: progress_bar_impl(max_value, message)
		, m_thread(
#if __cpp_lib_jthread >= 201911L
			  [this](std::stop_token stoken)
			  { this->run(stoken); }
#else
			  [this]()
			  { this->run(); }
#endif
		  )
	{
	}

	~fancy_progress_bar_impl() override;

#if __cpp_lib_jthread >= 201911L
	void run(std::stop_token stoken);
#else
	void run();
#endif

	void consumed(uint64_t n) override;
	void progress(uint64_t p) override;
	void message(const std::string &msg) override;

	void print_progress();
	void print_done() override;

	std::mutex m_mutex;
	std::condition_variable m_cv;

	float m_progress;
	uint32_t m_width, m_bar_width, m_height;
	uint32_t m_steps, m_last_steps = 0;
	uint64_t m_last_consumed = 0;
#if __cpp_lib_jthread >= 201911L
	std::jthread m_thread;
#else
	std::thread m_thread;
	bool m_stop = false;
#endif
};

const char *kBlocks[] = {
	" ",
	"▏",
	"▎",
	"▍",
	"▌",
	"▋",
	"▊",
	"▉",
	"█",
};

const size_t kBlockCount = sizeof(kBlocks) / sizeof(void *) - 1;

fancy_progress_bar_impl::~fancy_progress_bar_impl()
{
	using namespace std::literals;
	assert(m_thread.joinable());

#if __cpp_lib_jthread >= 201911L
	m_thread.request_stop();
#else
	m_stop = true;
#endif
	m_thread.join();
}

#if __cpp_lib_jthread >= 201911L
void fancy_progress_bar_impl::run(std::stop_token stoken)
#else
void fancy_progress_bar_impl::run()
#endif
{
	using namespace std::literals;

	bool printedAny = false;

	try
	{
		for (;;)
		{
			std::unique_lock lock(m_mutex);

			m_cv.wait_for(lock, 25ms);

#if __cpp_lib_jthread >= 201911L
			if (stoken.stop_requested())
				break;
#else
			if (m_stop)
				break;
#endif

			auto now = std::chrono::system_clock::now();

			if (m_consumed == m_last_consumed or now - m_start < 1s)
				continue;

			m_last_consumed = m_consumed;

			// See if we need to do work
			std::tie(m_width, m_height) = get_terminal_width_and_height();
			m_progress = static_cast<float>(m_consumed) / m_max_value;
			m_bar_width = 7 * m_width / 10; // 70% of the width of the terminal
			m_steps = static_cast<uint32_t>(std::ceil(m_progress * m_bar_width * kBlockCount));

			if (m_steps == m_last_steps)
				continue;

			m_last_steps = m_steps;

			// auto [w, h] = get_terminal_width_and_height();
			if (not printedAny)
				std::cout << std::format("\n\0337\033[{};{}r\0338\033[1A", 0, m_height - 1)
						  << std::flush;

			print_progress();

			printedAny = true;
		}
	}
	catch (const std::exception &ex)
	{
		std::cerr << "error finishing progress bar: " << ex.what() << '\n';
	}

	if (printedAny)
		print_done();
}

void fancy_progress_bar_impl::consumed(uint64_t n)
{
	progress_bar_impl::consumed(n);
	// m_cv.notify_one();
}

void fancy_progress_bar_impl::progress(uint64_t p)
{
	progress_bar_impl::progress(p);
	// m_cv.notify_one();
}

void fancy_progress_bar_impl::message(const std::string &msg)
{
	std::unique_lock lock(m_mutex);
	progress_bar_impl::message(msg);
	// m_cv.notify_one();
}

const uint32_t kMinBarWidth = 40, kMinMsgWidth = 12;

void fancy_progress_bar_impl::print_progress()
{
	const uint32_t pct_width = 5;
	uint32_t msg_width = m_width - m_bar_width - pct_width - 1;

	if (msg_width < kMinMsgWidth)
	{
		m_bar_width += kMinMsgWidth - msg_width;
		msg_width = kMinMsgWidth;
	}

	std::string bar;
	bar.reserve(m_bar_width * 4UL);

	for (uint32_t i = 0; i < m_bar_width; ++i)
	{
		if (i * kBlockCount <= m_steps)
			bar += kBlocks[kBlockCount];
		else if (i * kBlockCount > m_steps + kBlockCount)
			bar += kBlocks[0];
		else
			bar += kBlocks[1 + m_steps % kBlockCount];
	}

	// make the bar more colorfull
	struct color_type
	{
		uint8_t r, g, b;
	} fg{ 0, 3, 5 }, bg{ 0, 1, 2 };

	auto esc_1 = std::format("\033[38;5;{}m\033[48;5;{}m",
		16 + (fg.r * 36) + (fg.g * 6) + fg.b,
		16 + (bg.r * 36) + (bg.g * 6) + bg.b);
	std::string esc_2("\033[0m");

	bar = esc_1 + bar + esc_2;

	std::string msg = m_message.length() <= msg_width
	                      ? m_message
	                      : m_message.substr(0, msg_width - 3) + "...";

	std::cout << std::format("\0337\033[?25l\033[{};{}f{:{}} {} {:3d}%\033[?25h\0338", m_height, 1,
					 msg, msg_width,
					 bar,
					 static_cast<int>(std::ceil(m_progress * 100)))
			  << std::flush;
}

void fancy_progress_bar_impl::print_done()
{
	// wipe out progress bar first
	std::tie(m_width, m_height) = get_terminal_width_and_height();
	std::cout << std::format("\0337\033[{};{}H{}\033[{};{}r\0338", m_height, 0,
					 std::string(m_width, ' '), 0, m_height)
			  << std::flush;
	;
	progress_bar_impl::print_done();
}

// --------------------------------------------------------------------

progress_bar::progress_bar(int64_t max_value, const std::string &message)
{
	if (VERBOSE >= 0)
	{
		if (isatty(STDOUT_FILENO) and get_terminal_width() > kMinBarWidth)
			m_impl = new fancy_progress_bar_impl(max_value, message);
		else
			m_impl = new simple_progress_bar_impl(max_value, message);
	}
}

progress_bar::~progress_bar()
{
	flush();
}

void progress_bar::consumed(int64_t inConsumed)
{
	if (m_impl != nullptr)
		m_impl->consumed(inConsumed);
}

void progress_bar::progress(int64_t value)
{
	if (m_impl != nullptr)
		m_impl->progress(value);
}

void progress_bar::message(const std::string &message)
{
	if (m_impl != nullptr)
		m_impl->message(message);
}

void progress_bar::flush()
{
	if (m_impl)
	{
		delete m_impl;
		m_impl = nullptr;
	}
}

} // namespace cif

// --------------------------------------------------------------------
//
// Try to find a named resource. Can be either a local file with this name,
// a file located in our cache directory or a resource linked in with mrc.
//
// We have a special, private version of mrsrc here. To be able to create
// shared libraries and still be able to link when there's no mrc used.

namespace mrsrc
{
/// \brief Internal data structure as generated by mrc
struct rsrc_imp
{
	unsigned int m_next;
	unsigned int m_child;
	unsigned int m_name;
	unsigned int m_size;
	unsigned int m_data;
};
} // namespace mrsrc

#if _WIN32

# if __MINGW32__

extern "C" __attribute__((weak, alias("gResourceIndexDefault"))) const mrsrc::rsrc_imp gResourceIndex[];
extern "C" __attribute__((weak, alias("gResourceDataDefault"))) const char gResourceData[];
extern "C" __attribute__((weak, alias("gResourceNameDefault"))) const char gResourceName[];

# else

extern "C" const mrsrc::rsrc_imp *gResourceIndexDefault[1] = {};
extern "C" const char *gResourceDataDefault[1] = {};
extern "C" const char *gResourceNameDefault[1] = {};

extern "C" const mrsrc::rsrc_imp gResourceIndex[];
extern "C" const char gResourceData[];
extern "C" const char gResourceName[];

#  pragma comment(linker, "/alternatename:gResourceIndex=gResourceIndexDefault")
#  pragma comment(linker, "/alternatename:gResourceData=gResourceDataDefault")
#  pragma comment(linker, "/alternatename:gResourceName=gResourceNameDefault")

# endif

#else
extern const __attribute__((weak)) mrsrc::rsrc_imp gResourceIndex[];
extern const __attribute__((weak)) char gResourceData[];
extern const __attribute__((weak)) char gResourceName[];

const mrsrc::rsrc_imp gResourceIndex[1] = {};
const char gResourceData[1] = {};
const char gResourceName[1] = {};

#endif

namespace mrsrc
{
class rsrc_data
{
  public:
	static rsrc_data &instance()
	{
		static rsrc_data s_instance;
		return s_instance;
	}

	[[nodiscard]] const rsrc_imp *index() const { return m_index; }

	[[nodiscard]] const char *data(unsigned int offset) const
	{
		return m_data + offset;
	}

	[[nodiscard]] const char *name(unsigned int offset) const
	{
		return m_name + offset;
	}

  private:
	rsrc_data()
	{
		// if (gResourceIndex and (gResourceIndex[0].m_child > 0 or gResourceIndex[0].m_size > 0) and gResourceIndex and gResourceName)
		if (gResourceIndex[0].m_child > 0 or gResourceIndex[0].m_size > 0)
		{
			m_index = gResourceIndex;
			m_data = gResourceData;
			m_name = gResourceName;
		}
	}

	rsrc_imp m_dummy = {};
	const rsrc_imp *m_index = &m_dummy;
	const char *m_data = "";
	const char *m_name = "";
};

/// \brief Class mrsrc::rsrc contains a pointer to the data in the
/// resource, as well as offering an iterator interface to its
/// children.

class rsrc
{
  public:
	rsrc()
		: m_impl(rsrc_data::instance().index())
	{
	}

	rsrc(const rsrc &other) = default;
	rsrc &operator=(const rsrc &other) = default;

	rsrc(std::filesystem::path path);

	[[nodiscard]] std::string name() const { return rsrc_data::instance().name(m_impl->m_name); }

	[[nodiscard]] const char *data() const { return rsrc_data::instance().data(m_impl->m_data); }

	[[nodiscard]] unsigned long size() const { return m_impl->m_size; }

	explicit operator bool() const { return m_impl != nullptr and m_impl->m_size > 0; }

	template <typename RSRC>
	class iterator_t
	{
	  public:
		using iterator_category = std::input_iterator_tag;
		using value_type = RSRC;
		using difference_type = std::ptrdiff_t;
		using pointer = value_type *;
		using reference = value_type &;

		iterator_t(const rsrc_imp *cur)
			: m_cur(cur)
		{
		}

		iterator_t(const iterator_t &i) = default;
		iterator_t &operator=(const iterator_t &i) = default;

		reference operator*() { return m_cur; }
		pointer operator->() { return &m_cur; }

		iterator_t &operator++()
		{
			if (m_cur.m_impl->m_next)
				m_cur.m_impl = rsrc_data::instance().index() + m_cur.m_impl->m_next;
			else
				m_cur.m_impl = nullptr;
			return *this;
		}

		iterator_t operator++(int)
		{
			auto tmp(*this);
			this->operator++();
			return tmp;
		}

		bool operator==(const iterator_t &rhs) const { return m_cur.m_impl == rhs.m_cur.m_impl; }
		bool operator!=(const iterator_t &rhs) const { return m_cur.m_impl != rhs.m_cur.m_impl; }

	  private:
		value_type m_cur;
	};

	using iterator = iterator_t<rsrc>;

	[[nodiscard]] iterator begin() const
	{
		const rsrc_imp *impl = nullptr;
		if (m_impl and m_impl->m_child)
			impl = rsrc_data::instance().index() + m_impl->m_child;
		return { impl };
	}

	[[nodiscard]] iterator end() const
	{
		return { nullptr };
	}

  private:
	rsrc(const rsrc_imp *imp)
		: m_impl(imp)
	{
	}

	const rsrc_imp *m_impl;
};

inline rsrc::rsrc(std::filesystem::path p)
{
	m_impl = rsrc_data::instance().index();

	// using std::filesytem::path would have been natural here of course...

	auto pb = p.begin();
	auto pe = p.end();

	while (m_impl != nullptr and pb != pe)
	{
		auto name = *pb++;

		const rsrc_imp *impl = nullptr;
		for (rsrc child : *this)
		{
			if (child.name() == name)
			{
				impl = child.m_impl;
				break;
			}
		}

		m_impl = impl;
	}

	if (pb != pe) // not found
		m_impl = nullptr;
}

// --------------------------------------------------------------------

template <typename CharT, typename Traits>
class basic_streambuf : public std::basic_streambuf<CharT, Traits>
{
  public:
	using char_type = CharT;
	using traits_type = Traits;
	using int_type = typename traits_type::int_type;
	using pos_type = typename traits_type::pos_type;
	using off_type = typename traits_type::off_type;

	/// \brief constructor taking a \a path to the resource in memory
	basic_streambuf(const std::string &path)
		: m_rsrc(path)
	{
		init();
	}

	/// \brief constructor taking a \a rsrc
	basic_streambuf(const rsrc &rsrc)
		: m_rsrc(rsrc)
	{
		init();
	}

	basic_streambuf(const basic_streambuf &) = delete;

	basic_streambuf(basic_streambuf &&rhs) noexcept
		: basic_streambuf(rhs.m_rsrc)
	{
	}

	basic_streambuf &operator=(const basic_streambuf &) = delete;

	basic_streambuf &operator=(basic_streambuf &&rhs) noexcept
	{
		swap(rhs);
		return *this;
	}

	void swap(basic_streambuf &rhs) noexcept
	{
		std::swap(m_begin, rhs.m_begin);
		std::swap(m_end, rhs.m_end);
		std::swap(m_current, rhs.m_current);
	}

  private:
	void init()
	{
		m_begin = reinterpret_cast<const char_type *>(m_rsrc.data());
		m_end = reinterpret_cast<const char_type *>(m_rsrc.data() + m_rsrc.size());
		m_current = m_begin;
	}

	int_type underflow() override
	{
		if (m_current == m_end)
			return traits_type::eof();

		return traits_type::to_int_type(*m_current);
	}

	int_type uflow() override
	{
		if (m_current == m_end)
			return traits_type::eof();

		return traits_type::to_int_type(*m_current++);
	}

	int_type pbackfail(int_type ch) override
	{
		if (m_current == m_begin or (ch != traits_type::eof() and ch != m_current[-1]))
			return traits_type::eof();

		return traits_type::to_int_type(*--m_current);
	}

	std::streamsize showmanyc() override
	{
		assert(std::less_equal<const char *>()(m_current, m_end));
		return m_end - m_current;
	}

	pos_type seekoff(off_type off, std::ios_base::seekdir dir, std::ios_base::openmode which) override
	{
		switch (dir)
		{
			case std::ios_base::beg:
				m_current = m_begin + off;
				break;

			case std::ios_base::end:
				m_current = m_end + off;
				break;

			case std::ios_base::cur:
				m_current += off;
				break;

			default:
				break;
		}

		if (m_current < m_begin)
			m_current = m_begin;

		if (m_current > m_end)
			m_current = m_end;

		return m_current - m_begin;
	}

	pos_type seekpos(pos_type pos, std::ios_base::openmode which) override
	{
		m_current = m_begin + pos;

		if (m_current < m_begin)
			m_current = m_begin;

		if (m_current > m_end)
			m_current = m_end;

		return m_current - m_begin;
	}

  private:
	rsrc m_rsrc;
	const char_type *m_begin;
	const char_type *m_end;
	const char_type *m_current;
};

using streambuf = basic_streambuf<char, std::char_traits<char>>;

// --------------------------------------------------------------------
// class mrsrc::istream

template <typename CharT, typename Traits>
class basic_istream : public std::basic_istream<CharT, Traits>
{
  public:
	using char_type = CharT;
	using traits_type = Traits;
	using int_type = typename traits_type::int_type;
	using pos_type = typename traits_type::pos_type;
	using off_type = typename traits_type::off_type;

  private:
	using _streambuf_type = basic_streambuf<CharT, Traits>;
	using _istream_type = std::basic_istream<CharT, Traits>;

	_streambuf_type m_buffer;

  public:
	basic_istream(const std::string &path)
		: _istream_type(&m_buffer)
		, m_buffer(path)
	{
		this->init(&m_buffer);
	}

	basic_istream(rsrc &resource)
		: _istream_type(&m_buffer)
		, m_buffer(resource)
	{
		this->init(&m_buffer);
	}

	basic_istream(const basic_istream &) = delete;

	basic_istream(basic_istream &&rhs)
		: _istream_type(std::move(rhs))
		, m_buffer(std::move(rhs.m_buffer))
	{
		_istream_type::set_rdbuf(&m_buffer);
	}

	basic_istream &operator=(const basic_istream &) = delete;

	basic_istream &operator=(basic_istream &&rhs)
	{
		_istream_type::operator=(std::move(rhs));
		m_buffer = std::move(rhs.m_buffer);
		return *this;
	}

	void swap(basic_istream &rhs)
	{
		_istream_type::swap(rhs);
		m_buffer.swap(rhs.m_buffer);
	}

	_streambuf_type *rdbuf() const
	{
		return const_cast<_streambuf_type *>(&m_buffer);
	}
};

using istream = basic_istream<char, std::char_traits<char>>;
} // namespace mrsrc

// --------------------------------------------------------------------

namespace cif
{

// --------------------------------------------------------------------

class resource_pool
{
  public:
	static resource_pool &instance()
	{
		static std::unique_ptr<resource_pool> s_instance(new resource_pool);
		return *s_instance;
	}

	void pushDir(fs::path dir)
	{
		std::error_code ec;

		if (fs::exists(dir, ec) and not ec)
			mDirs.push_front(dir);
	}

	void pushDir(const char *path)
	{
		if (path != nullptr)
			pushDir(fs::path(path));
	}

	void pushAlias(const std::string &name, std::filesystem::path dataFile)
	{
		std::error_code ec;
		if (not fs::exists(dataFile, ec) or ec)
			throw std::runtime_error("Attempt to add a file resource for " + name + " that cannot be used (" + dataFile.string() + ") :" + ec.message());

		mLocalResources[name] = dataFile;
	}

	std::unique_ptr<std::istream> load(fs::path name);

	const auto &data_directories() { return mDirs; }
	const auto &file_resources() { return mLocalResources; }

  private:
	resource_pool();

	std::unique_ptr<std::ifstream> open(fs::path &p)
	{
		std::unique_ptr<std::ifstream> result;

		std::error_code ec;
		if (fs::exists(p, ec))
		{
			std::unique_ptr<std::ifstream> file = std::make_unique<std::ifstream>(p, std::ios::binary);
			if (file->is_open())
				result.reset(file.release());
		}
		if (ec != std::errc{} or result == nullptr)
			std::cerr << "Error opening resource file " << std::quoted(p.string()) << '\n';

		return result;
	}

	std::map<std::string, std::filesystem::path> mLocalResources;
	std::deque<fs::path> mDirs;
};

resource_pool::resource_pool()
{
	// directories are searched in reverse order

	// As a last resort, try the location that might have been
	// used during installation, works only when running on an
	// OS with a proc file system.

	std::error_code ec;
	if (auto exefile = fs::read_symlink("/proc/self/exe", ec); not ec and exefile.parent_path().filename() == "bin")
	{
		auto install_prefix = exefile.parent_path().parent_path();
		auto data_dir = install_prefix / "share" / "libcifpp";
		if (fs::exists(data_dir, ec))
			pushDir(data_dir);
	}

#if defined(DATA_DIR)
	pushDir(DATA_DIR);
#endif

	pushDir(getenv("LIBCIFPP_DATA_DIR"));

	auto ccp4 = getenv("CCP4");
	if (ccp4 != nullptr)
		pushDir(fs::path(ccp4) / "share" / "libcifpp");

#if defined(CACHE_DIR)
	pushDir(CACHE_DIR);
#endif
}

std::unique_ptr<std::istream> resource_pool::load(fs::path name)
{
	std::unique_ptr<std::istream> result;
	std::error_code ec;

	fs::path p = name;

	if (mLocalResources.count(name.string()))
		result = open(mLocalResources[name.string()]);

	if (fs::exists(p, ec) and not ec)
		result = open(p);

	for (auto di = mDirs.begin(); not result and di != mDirs.end(); ++di)
	{
		auto p2 = *di / p;
		if (fs::exists(p2, ec) and not ec)
			result = open(p2);
	}

	// if (not result and gResourceData)
	if (not result and (gResourceIndex[0].m_child > 0 or gResourceIndex[0].m_size > 0))
	{
		mrsrc::rsrc rsrc(name);
		if (rsrc)
			result = std::make_unique<mrsrc::istream>(rsrc);
	}

	return result;
}

// --------------------------------------------------------------------

void add_data_directory(std::filesystem::path dataDir)
{
	resource_pool::instance().pushDir(dataDir);
}

void add_file_resource(const std::string &name, std::filesystem::path dataFile)
{
	resource_pool::instance().pushAlias(name, dataFile);
}

std::unique_ptr<std::istream> load_resource(std::filesystem::path name)
{
	return resource_pool::instance().load(name);
}

void list_file_resources(std::ostream &os)
{
	auto &file_resources = resource_pool::instance().file_resources();

	if (not file_resources.empty())
	{
		os << "\nThe following named resources were loaded:\n";
		for (const auto &[name, path] : file_resources)
			os << name << " -> " << std::quoted(path.string()) << '\n';
	}
}

void list_data_directories(std::ostream &os)
{
	for (auto &p : resource_pool::instance().data_directories())
		os << p << '\n';
}

} // namespace cif
