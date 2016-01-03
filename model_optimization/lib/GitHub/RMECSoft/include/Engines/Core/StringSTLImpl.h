/*
String.h
Written by Matthew Fisher

String class and related conversions.
*/

#pragma once

class String
{
public:
    struct LexicographicComparison
    {
        bool operator()(const String &L, const String &R) const
        {
            return strcmp(L.CString(), R.CString()) < 0;
        }
    };

    String()
    {
        FreeMemory();
    }

    ~String()
    {
        
    }

    //
    // All constructors are made explicit to avoid undesired auto-casting.
    //
    String(const String &S)
    {
        FreeMemory();
        *this = S;
    }
    explicit String(UINT Value)
    {
        FreeMemory();
        *this = Value;
    }
    explicit String(UINT64 Value)
    {
        FreeMemory();
        *this = Value;
    }
    explicit String(DWORD Value)
    {
        FreeMemory();
        *this = int(Value);
    }
    explicit String(int Value)
    {
        FreeMemory();
        *this = Value;
    }
    explicit String(float Value)
    {
        FreeMemory();
        *this = Value;
    }
    explicit String(double Value)
    {
        FreeMemory();
        *this = float(Value);
    }
    String(const char *Text)
    {
        FreeMemory();
        *this = Text;
    }
    explicit String(char Character)
    {
        FreeMemory();
        *this = Character;
    }
    explicit String(bool B)
    {
        FreeMemory();
        if(B)
        {
            *this = "true";
        }
        else
        {
            *this = "false";
        }
    }

    //
    // Memory
    //
    __forceinline void FreeMemory()
    {
        _Data.Allocate(1);
        _Data[0] = '\0';
    }

    __forceinline void ReSize(UINT NewLength)
    {
        _Data.ReSize(NewLength + 1);
        _Data.Last() = '\0';
    }

    //
    // Accessors
    //
    __forceinline char* CString()
    {
        return _Data.CArray();
    }

    __forceinline const char* CString() const
    {
        return _Data.CArray();
    }

    __forceinline UINT Length() const
    {
        return _Data.Length() - 1;
    }

    __forceinline char Last() const
    {
        return _Data[_Data.Length() - 2];
    }

    __forceinline char& operator [] (UINT k)
    {
        return _Data[k];
    }

    __forceinline char& operator [] (int k) 
    {
        return _Data[k];
    }

    __forceinline const char& operator [] (UINT k) const
    {
        return _Data[k];
    }

    __forceinline const char& operator [] (int k) const
    {
        return _Data[k];
    }

    //
    // Conversions
    //
    int ConvertToInteger() const
    {
        return atoi(CString());
    }

    ULONG ConvertToUnsignedInteger() const
    {
        return ULONG(_atoi64(CString()));
    }

    LONGLONG ConvertToLongInteger() const
    {
        return _atoi64(CString());
    }

    double ConvertToDouble() const
    {
        return atof(CString());
    }

    bool ConvertToBoolean() const;

    float ConvertToFloat() const
    {
        return float(atof(CString()));
    }

    //
    // Query
    //
    bool ExactMatchAtOffset(const String &Find, UINT Offset) const;
    bool Contains(const String &Find) const;
    bool IsNumeric() const;
    bool IsSuffix(const String &EndCanidate) const;
    bool IsPrefix(const String &StartCanidate) const;
    void Partition(char Seperator, Vector<String> &Output) const;
    Vector<String> Partition(char Seperator) const;
    void PartitionAboutIndex(UINT Index, String &Left, String &Right) const;
    void Partition(const String &Seperator, Vector<String> &Output) const;
    int FindFirstIndex(char Seperator) const;
    int FindLastIndex(char Seperator) const;
    String FindAndReplace(const String &Find, const String &Replace) const;
    String MakeLowercase() const;
    String RemoveSuffix(const String &EndCandidate) const;
    UINT Hash() const;
    
    //
    // Assignment
    //
    String& operator = (const char *String);
    String& operator = (char Character);
    String& operator = (float Value);
    String& operator = (int Value);
    String& operator = (UINT Value);
    String& operator = (UINT64 Value);
    String& operator = (const String &S);

    //
    // Modifiers
    //
    String& operator += (const String &S);
    __forceinline String& operator += (char C)
    {
        PushEnd(C);
        return *this;
    }

    __forceinline void PushEnd(char C)
    {
        _Data.Last() = C;
        _Data.PushEnd('\0');
    }

    __forceinline void PopEnd()
    {
        _Data.PopEnd();
        _Data.Last() = '\0';
    }

    void PopFront()
    {
        _Data.PopFront();
    }

    //
    // Formatting
    //
    static String UnsignedIntAsHex(UINT Value);
    static String ZeroPad(const String &S, UINT ZeroPadding);

private:
    void ResizeToCStringLength();

    Vector<char> _Data;
};

//
// String Comparison
//
bool operator == (const String &L, const String &R);
bool operator == (const char *L, const String &R);
bool operator == (const String &R, const char *L);
__forceinline bool operator != (const String &L, const String &R) {return !(L == R);}
__forceinline bool operator != (const char *L, const String &R) {return !(L == R);}
__forceinline bool operator != (const String &R, const char *L) {return !(L == R);}

//
// String Operations
//
String operator + (const String &L, const String &R);
String operator + (const char *L, const String &R);
String operator + (const String &R, const char *L);
ostream& operator << (ostream &os, const String &S);

