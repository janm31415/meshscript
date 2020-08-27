#pragma once

#ifdef _WIN32
#if defined(libpoisson_EXPORTS)
#  define LIBPOISSON_API __declspec(dllexport)
#else
#  define LIBPOISSON_API __declspec(dllimport)
#endif
#else
#define LIBPOISSON_API
#endif