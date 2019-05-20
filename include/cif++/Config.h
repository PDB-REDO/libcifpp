// Lib for working with structures as contained in mmCIF and PDB files

#pragma once

#include <string>
#include <cstdint>

#define HAVE_CPP0X_TEMPLATE_ALIASES 1
#define HAVE_CPP0X_VARIADIC_TEMPLATES 1
#define HAVE_CPP0X_INITIALIZER_LISTS 1

#if defined(_MSC_VER)

// These are Microsoft Visual C++ special settings
// the iso646 file contains the C++ keywords that are
// otherwise not recognized.
#include <ciso646>
#define snprintf _snprintf

// Disable some warnings
#pragma warning (disable : 4996)
#pragma warning (disable : 4355)
#endif
