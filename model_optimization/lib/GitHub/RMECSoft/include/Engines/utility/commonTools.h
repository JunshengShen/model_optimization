#pragma once
/**
Common utility functions.
*/
#include <strstream>
#include <cstdarg>

namespace SU{

	void Output(const _tString &message);
	std::wstring s2ws(const std::string& s);
	std::string ws2s(const std::wstring& s);
	int WstringToInt( const std::wstring ws );
	double WstringToDouble(const std::wstring ws);
	int StringToInt(const std::string s);
	/**
	* Format a string object
	*/
	_tString formatString(const tChar* format, ...);

	/**
	* Converts a type to a _tString.
	* Convert a type to a string by streaming it. Requires that there's an ostream
	* inserter available for type T.
	*/
	template <class T>
	_tString ToString(T num)
	{
		_tString strNum = _T("");

		{
			std::strstream buf;

			buf << num << std::ends;

#ifdef _UNICODE
			std::string temp = buf.str();

			USES_CONVERSION;

			strNum = A2W(temp.c_str());
#else
			strNum = buf.str();
#endif
			buf.freeze(false);
		}

		return strNum;
	}

	template <class T>
	bool ToBool(const T &value)
	{
		return (0 != value);
	}

	inline bool BOOL_to_bool(const BOOL bResult)
	{
		// Convert a make believe BOOL into a real bool.
		// Removes warning C4800...

		return (TRUE == bResult);
	}

	_tString HexToString(   const BYTE *pBuffer, size_t iBytes);

	void StringToHex(   const _tString &str,    BYTE *pBuffer,    size_t nBytes);

	_tString GetLastErrorMessage(   DWORD last_error);

	HMODULE  GetCurrentModule();                              //!<Get handle of current EXE/DLL 
	_tString GetModuleFileName(HINSTANCE hModule = NULL);
	_tString GetModuleDir(HMODULE hModule = NULL);            //!<��ȡִ�г���Ŀ¼
	_tString GetCurrentDirectory();      //!<��ȡ����ʱĿ¼
	
	_tString GetDateStamp();
	_tString GetCurTime();
	_tString GetTimeStamp();             //!<����YYYYMMDDHHMMSS


	_tString ToHex(BYTE c);

	_tString GetComputerName();        	
	_tString GetUserName();
	
	_tString GetFileVersion();

	_tString StripLeading(   const _tString &source,    const char toStrip);

	_tString StripTrailing(   const _tString &source,    const char toStrip);

	std::string GBK2UTF(std::string gbkStr);

	

}//end name space SU