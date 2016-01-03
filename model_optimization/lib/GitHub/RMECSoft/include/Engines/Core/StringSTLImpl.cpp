/*
String.cpp
Written by Matthew Fisher

Implementation of String class.
*/

String String::UnsignedIntAsHex(UINT Value)
{
    String S;
    S.ReSize(20);
    sprintf(S.CString(), "%x", Value);
    S.ResizeToCStringLength();
    return S;
}

String String::ZeroPad(const String &S, UINT ZeroPadding)
{
    String Prefix;
    for(UINT i = S.Length(); i < ZeroPadding; i++)
    {
        Prefix.PushEnd('0');
    }
    return Prefix + S;
}

bool String::ConvertToBoolean() const
{
    if(*this == "true")
    {
        return true;
    }
    if(*this == "false")
    {
        return false;
    }
    SignalError("Invalid string -> boolean conversion");
    return false;
}

void String::ResizeToCStringLength()
{
    for(UINT i = 0; i < _Data.Length(); i++)
    {
        if(_Data[i] == '\0')
        {
            _Data.ReSize(i + 1);
            return;
        }
    }
    SignalError("No null terminator found in String::ResizeToCStringLength");
}

String& String::operator = (float Value)
{
    ostringstream S(ostringstream::out);
    S << Value;
    *this = S.str().c_str();
    return *this;
}

String& String::operator = (char Character)
{
    _Data.ReSize(2);
    _Data[0] = Character;
    _Data[1] = '\0';
    return *this;
}

String& String::operator = (int Value)
{
    ostringstream S(ostringstream::out);
    S << Value;
    *this = S.str().c_str();
    return *this;
}

String& String::operator = (UINT Value)
{
    ostringstream S(ostringstream::out);
    S << Value;
    *this = S.str().c_str();
    return *this;
}

String& String::operator = (UINT64 Value)
{
    ostringstream S(ostringstream::out);
    S << Value;
    *this = S.str().c_str();
    return *this;
}

String& String::operator = (const char *String)
{
    UINT NewLength = UINT(strlen(String)) + 1;
    _Data.ReSize(NewLength);
    memcpy(CString(), String, NewLength);
    return *this;
}

String& String::operator = (const String &S)
{
    _Data = S._Data;
    return *this;
}

String& String::operator += (const String &S)
{
    UINT BaseLength = Length();
    _Data.ReSize(BaseLength + S.Length() + 1);
    for(UINT i = 0; i < S.Length(); i++)
    {
        _Data[i + BaseLength] = S[i];
    }
    _Data.Last() = '\0';
    return *this;
}

String String::MakeLowercase() const
{
    String Result = *this;
    for(UINT i = 0; i < Result.Length(); i++)
    {
        Result[i] = tolower(Result[i]);
    }
    return Result;
}

void String::Partition(char Seperator, Vector<String> &Output) const
{
    Output.FreeMemory();
    String CurEntry;
    for(UINT CurIndex = 0; CurIndex < Length(); CurIndex++)
    {
        if((*this)[CurIndex] == Seperator)
        {
            if(CurEntry.Length() > 0)
            {
                Output.PushEnd(CurEntry);
                CurEntry.FreeMemory();
            }
        }
        else
        {
            CurEntry.PushEnd((*this)[CurIndex]);
        }
    }
    if(CurEntry.Length() > 0)
    {
        Output.PushEnd(CurEntry);
    }
}

Vector<String> String::Partition(char Seperator) const
{
    Vector<String> Result;
    Partition(Seperator, Result);
    return Result;
}

int String::FindFirstIndex(char Seperator) const
{
    UINT L = Length();
    for(UINT CurIndex = 0; CurIndex < L; CurIndex++)
    {
        if((*this)[CurIndex] == Seperator)
        {
            return CurIndex;
        }
    }
    return -1;
}

int String::FindLastIndex(char Seperator) const
{
    UINT L = Length();
    for(int CurIndex = L - 1; CurIndex >= 0; CurIndex--)
    {
        if((*this)[CurIndex] == Seperator)
        {
            return CurIndex;
        }
    }
    return -1;
}

void String::PartitionAboutIndex(UINT Index, String &Left, String &Right) const
{
    UINT L = Length();
    Assert(Index < L, "Index invalid in String::PartitionAboutIndex");
    Left.ReSize(Index);
    Right.ReSize(L - Index - 1);
    for(UINT LeftIndex = 0; LeftIndex < Index; LeftIndex++)
    {
        Left[LeftIndex] = (*this)[LeftIndex];
    }
    for(UINT RightIndex = 0; RightIndex < Right.Length(); RightIndex++)
    {
        Right[RightIndex] = (*this)[Index + 1 + RightIndex];
    }
}

void String::Partition(const String &Seperator, Vector<String> &Output) const
{
    Assert(Seperator.Length() >= 1, "empty seperator");
    Output.FreeMemory();
    String CurEntry;
    for(UINT CurIndex = 0; CurIndex < Length(); CurIndex++)
    {
        bool IsSeperator = true;
        for(UINT TestIndex = 0; TestIndex < Seperator.Length() && CurIndex + TestIndex < Length() && IsSeperator; TestIndex++)
        {
            if((*this)[CurIndex + TestIndex] != Seperator[TestIndex])
            {
                IsSeperator = false;
            }
        }

        if(IsSeperator)
        {
            if(CurEntry.Length() > 0)
            {
                Output.PushEnd(CurEntry);
                CurEntry.FreeMemory();
            }
            CurIndex += Seperator.Length() - 1;
        }
        else
        {
            CurEntry.PushEnd((*this)[CurIndex]);
        }
    }
    if(CurEntry.Length() > 0)
    {
        Output.PushEnd(CurEntry);
    }
}

bool String::ExactMatchAtOffset(const String &Find, UINT Offset) const
{
    UINT MatchLength = 0;
    for(UINT Index = 0; Index + Offset < Length() && Index < Find.Length(); Index++)
    {
        if((*this)[Index + Offset] == Find[Index])
        {
            MatchLength++;
            if(MatchLength == Find.Length())
            {
                return true;
            }
        }
    }
    return false;
}

bool String::Contains(const String &Find) const
{
    for(UINT Index = 0; Index < Length(); Index++)
    {
        if(ExactMatchAtOffset(Find, Index))
        {
            return true;
        }
    }
    return false;
}

String String::FindAndReplace(const String &Find, const String &Replace) const
{
    String Result;
    for(UINT Index = 0; Index < Length(); Index++)
    {
        if(ExactMatchAtOffset(Find, Index))
        {
            Result += Replace;
            Index += Find.Length() - 1;
        }
        else
        {
            Result += (*this)[Index];
        }
    }
    return Result;
}

bool String::IsNumeric() const
{
    return (Length() > 0 &&
           (ConvertToDouble() != 0.0 ||
            _Data[0] == '0' ||
            _Data[0] == '.' ||
            _Data[0] == '-'));
}

bool String::IsSuffix(const String &EndCanidate) const
{
    if(Length() < EndCanidate.Length())
    {
        return false;
    }
    for(UINT i = 0; i < EndCanidate.Length(); i++)
    {
        if((*this)[Length() - EndCanidate.Length() + i] != EndCanidate[i])
        {
            return false;
        }
    }
    return true;
}

String String::RemoveSuffix(const String &EndCandidate) const
{
    Assert(IsSuffix(EndCandidate), "Removing invalid suffix");
    String Result = *this;
    for(UINT i = 0; i < EndCandidate.Length(); i++)
    {
        Result.PopEnd();
    }
    return Result;
}

bool String::IsPrefix(const String &StartCanidate) const
{
    if(Length() < StartCanidate.Length())
    {
        return false;
    }
    for(UINT i = 0; i < StartCanidate.Length(); i++)
    {
        if((*this)[i] != StartCanidate[i])
        {
            return false;
        }
    }
    return true;
}

UINT String::Hash() const
{
    return Utility::Hash32((const BYTE *)_Data.CArray(), _Data.Length());
}

//
// Overloaded operators
//

String operator + (const String &L, const String &R)
{
    String Ret = L;
    for(UINT i = 0; i < R.Length(); i++)
    {
        Ret.PushEnd(R[i]);
    }
    return Ret;
}

String operator + (const char *L, const String &R)
{
    String Ret = L;
    for(UINT i = 0; i < R.Length(); i++)
    {
        Ret.PushEnd(R[i]);
    }
    return Ret;
}

String operator + (const String &L, const char *R)
{
    UINT RLen = UINT(strlen(R));
    String Ret = L;
    for(UINT i = 0; i < RLen; i++)
    {
        Ret.PushEnd(R[i]);
    }
    return Ret;
}

ostream& operator << (ostream &os, const String &S)
{
    if(S.Length() > 0)
    {
        os.write(S.CString(), S.Length());
    }
    return os;
}

bool operator == (const String &L, const String &R)
{
    UINT LLen = L.Length();
    UINT RLen = R.Length();
    if(LLen != RLen)
    {
        return false;
    }
    for(UINT i = 0; i < LLen; i++)
    {
        if(L[i] != R[i])
        {
            return false;
        }
    }
    return true;
}

bool operator == (const char *L, const String &R)
{
    UINT LLen = UINT(strlen(L));
    UINT RLen = R.Length();
    if(LLen != RLen)
    {
        return false;
    }
    for(UINT i = 0; i < LLen; i++)
    {
        if(L[i] != R[i])
        {
            return false;
        }
    }
    return true;
}

bool operator == (const String &R, const char *L)
{
    return (L == R);
}
