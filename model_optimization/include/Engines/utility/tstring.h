#pragma once
#include <string>
#include <tchar.h>

#ifdef _UNICODE
#define CHAR_SIZE  2
#define LL         L
#define ToLongNumber wcstol
typedef std::wstring _tString;
typedef wchar_t    tChar;
#else
#define CHAR_SIZE  1
#define LL
#define ToLongNumber strtol
typedef std::string _tString;
typedef char       tChar;
#endif