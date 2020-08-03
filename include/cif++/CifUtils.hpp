// cif parsing library

#pragma once

#include <vector>
#include <set>
#include <cassert>

#include "cif++/Config.hpp"

struct rsrc_imp
{
	unsigned int	m_next;
	unsigned int	m_child;
	unsigned int	m_name;
	unsigned int	m_size;
	unsigned int	m_data;
};

#if USE_RSRC
extern const rsrc_imp gResourceIndex[];
extern const char gResourceData[];
extern const char gResourceName[];
#endif

namespace cif
{

// some basic utilities: Since we're using ASCII input only, we define for optimisation
// our own case conversion routines.

bool iequals(const std::string& a, const std::string& b);
int icompare(const std::string& a, const std::string& b);

bool iequals(const char* a, const char* b);
int icompare(const char* a, const char* b);

void toLower(std::string& s);
std::string toLowerCopy(const std::string& s);

// To make life easier, we also define iless and iset using iequals

struct iless
{
	bool operator()(const std::string& a, const std::string& b) const
	{
		return icompare(a, b) < 0;
	}
};

typedef std::set<std::string, iless>	iset;

// --------------------------------------------------------------------
// This really makes a difference, having our own tolower routines

extern const uint8_t kCharToLowerMap[256];

inline char tolower(char ch)
{
	return static_cast<char>(kCharToLowerMap[static_cast<uint8_t>(ch)]);
}

// --------------------------------------------------------------------

std::tuple<std::string,std::string> splitTagName(const std::string& tag);

// --------------------------------------------------------------------
//	custom wordwrapping routine

std::vector<std::string> wordWrap(const std::string& text, unsigned int width);

// --------------------------------------------------------------------
//	Code helping with terminal i/o

uint32_t get_terminal_width();

// --------------------------------------------------------------------
//	Path of the current executable

std::string get_executable_path();

// --------------------------------------------------------------------
//	some manipulators to write coloured text to terminals

enum StringColour {
	scBLACK = 0, scRED, scGREEN, scYELLOW, scBLUE, scMAGENTA, scCYAN, scWHITE, scNONE = 9 };

template<typename String, typename CharT>
struct ColouredString
{
	static_assert(std::is_reference<String>::value or std::is_pointer<String>::value, "String type must be pointer or reference");
	
	ColouredString(String s, StringColour fore, StringColour back, bool bold = true)
		: m_s(s), m_fore(fore), m_back(back), m_bold(bold) {}
	
	ColouredString& operator=(const ColouredString&) = delete;
	
	String m_s;
	StringColour m_fore, m_back;
	bool m_bold;
};

template<typename CharT, typename Traits>
std::basic_ostream<CharT, Traits>& operator<<(std::basic_ostream<CharT, Traits>& os, const ColouredString<const CharT*,CharT>& s)
{
	if (isatty(STDOUT_FILENO))
	{
		std::basic_ostringstream<CharT, Traits> ostr;
		ostr << "\033[" << (30 + s.m_fore) << ';' << (s.m_bold ? "1" : "22") << ';' << (40 + s.m_back) << 'm'
			 << s.m_s
			 << "\033[0m";
		
		return os << ostr.str();
	}
	else
		return os << s.m_s;
}

template<typename CharT, typename Traits, typename String>
std::basic_ostream<CharT, Traits>& operator<<(std::basic_ostream<CharT, Traits>& os, const ColouredString<String,CharT>& s)
{
	if (isatty(STDOUT_FILENO))
	{
		std::basic_ostringstream<CharT, Traits> ostr;
		ostr << "\033[" << (30 + s.m_fore) << ';' << (s.m_bold ? "1" : "22") << ';' << (40 + s.m_back) << 'm'
			 << s.m_s
			 << "\033[0m";
		
		return os << ostr.str();
	}
	else
		return os << s.m_s;
}

template<typename CharT>
inline auto coloured(const CharT* s, StringColour fore = scWHITE, StringColour back = scRED, bool bold = true)
{
	return ColouredString<const CharT*, CharT>(s, fore, back, bold);
}

template<typename CharT, typename Traits, typename Alloc>
inline auto coloured(const std::basic_string<CharT, Traits, Alloc>& s, StringColour fore = scWHITE, StringColour back = scRED, bool bold = true)
{
	return ColouredString<const std::basic_string<CharT, Traits, Alloc>, CharT>(s, fore, back, bold);
}

template<typename CharT, typename Traits, typename Alloc>
inline auto coloured(std::basic_string<CharT, Traits, Alloc>& s, StringColour fore = scWHITE, StringColour back = scRED, bool bold = true)
{
	return ColouredString<std::basic_string<CharT, Traits, Alloc>, CharT>(s, fore, back, bold);
}

// --------------------------------------------------------------------
//	A progress bar

class Progress
{
  public:
				Progress(int64_t inMax, const std::string& inAction);
	virtual		~Progress();
	
	void		consumed(int64_t inConsumed);	// consumed is relative
	void		progress(int64_t inProgress);	// progress is absolute

	void		message(const std::string& inMessage);

  private:
				Progress(const Progress&);
	Progress&	operator=(const Progress&);

	struct ProgressImpl*	mImpl;
};

// --------------------------------------------------------------------
// The new default is to load 'resource' files from the file system
// since not everyone likes mrc. But if you do want to use mrc to load
// resources, specify the compile time flag USE_RSRC.

/// \brief Simple class containing the data for a resource
class rsrc
{
  public:
	rsrc()
		: m_data(nullptr), m_size(0) {}
	rsrc(const char* data, size_t size)
		: m_data(data), m_size(size) {}
	rsrc(const rsrc& rhs)
		: m_data(rhs.m_data), m_size(rhs.m_size) {}

	rsrc& operator=(const rsrc& rhs)
	{
		m_data = rhs.m_data;
		m_size = rhs.m_size;
		return *this;
	}

	explicit operator bool()
	{
		return m_data != nullptr and m_size > 0;
	}

	const char* data() const			{ return m_data; }
	size_t size() const					{ return m_size; }

  private:
	const char* m_data;
	size_t m_size;
};

/// \brief loader types
enum class rsrc_loader_type { mrsrc, file };

/// \brief loader info
struct rsrc_loader_info
{
	rsrc_loader_type type;
	std::string info;
	const void* ptrs[3];
};

class rsrc_loader
{
  public:

	static void init(std::initializer_list<rsrc_loader_info> loaders =
	{
#if USE_RSRC
		{ rsrc_loader_type::mrsrc, "", { gResourceIndex, gResourceData, gResourceName } },
#endif
		{ rsrc_loader_type::file, "." }
	})
	{
		assert(not s_instance);
		s_instance.reset(new rsrc_loader(loaders));
	}

	static rsrc load(const std::string& name)
	{
		assert(s_instance);
		if (not s_instance)
			init();

		return s_instance->do_load(name);
	}

	~rsrc_loader();

  private:

	rsrc_loader(const std::initializer_list<rsrc_loader_info>& loaders);

	rsrc do_load(const std::string& name);

	static std::unique_ptr<rsrc_loader> s_instance;
	
	std::list<struct rsrc_loader_impl*> m_loaders;
};

}


