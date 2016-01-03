/*
Console.cpp
Written by Matthew Fisher

Functions for interacting with the console.
*/

namespace Console
{
    ofstream ConsoleFile("Console.txt");

    ofstream& File()
    {
        return ConsoleFile;
    }

    void AdvanceLine()
    {
        WriteLine("");
    }

    void WriteLine(const String &S)
    {
        HANDLE StdHandle = GetStdHandle(STD_OUTPUT_HANDLE);
        DWORD BytesWritten;
        String FinalString = S + String("\n");
        WriteConsole(StdHandle, FinalString.CString(), FinalString.Length(), &BytesWritten, NULL);
        ConsoleFile << S << endl;
        ConsoleFile.flush();
    }

    void OverwriteLine(const String &S)
    {
        const UINT ConsoleWidth = 79;
        HANDLE StdHandle = GetStdHandle(STD_OUTPUT_HANDLE);
        CONSOLE_SCREEN_BUFFER_INFO CursorInfo;
        GetConsoleScreenBufferInfo(StdHandle, &CursorInfo);
        CursorInfo.dwCursorPosition.X = 0;
        CursorInfo.dwCursorPosition.Y = Math::Max(CursorInfo.dwCursorPosition.Y - 1, 0);
        SetConsoleCursorPosition(StdHandle, CursorInfo.dwCursorPosition);
        DWORD BytesWritten;
        String FinalString = S;
        while(FinalString.Length() < ConsoleWidth)
        {
            FinalString.PushEnd(' ');
        }
        FinalString.PushEnd('\n');
        WriteConsole(StdHandle, FinalString.CString(), FinalString.Length(), &BytesWritten, NULL);
        ConsoleFile << S << endl;
        ConsoleFile.flush();
    }

    void WriteString(const String &S)
    {
        HANDLE StdHandle = GetStdHandle(STD_OUTPUT_HANDLE);
        DWORD BytesWritten;
        WriteConsole(StdHandle, S.CString(), S.Length(), &BytesWritten, NULL);
        ConsoleFile << S;
        ConsoleFile.flush();
    }
}