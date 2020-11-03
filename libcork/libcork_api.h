#pragma once

#ifdef _WIN32
#if defined(libcork_EXPORTS)
#  define LIBCORK_API __declspec(dllexport)
#else
#  define LIBCORK_API __declspec(dllimport)
#endif
#else
#define LIBCORK_API
#endif