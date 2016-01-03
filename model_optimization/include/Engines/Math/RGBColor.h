/*
RGBColor.h
Written by Matthew Fisher

Generic 32-bit RGBA triplet structure.  Rather self-explanitory
*/

struct RGBColor
{
    RGBColor();
    RGBColor(BYTE _r, BYTE _g, BYTE _b);
    RGBColor(BYTE _r, BYTE _g, BYTE _b, BYTE _a);
    explicit RGBColor(const Vec3f &V);
    explicit RGBColor(const Vec4f &V);

    String CommaSeperatedString() const;
    
#ifdef USE_D3D
    operator D3DCOLOR();
#endif

	friend bool operator < (const struct RGBColor &ls, const struct RGBColor &rs)
	{
		return ls.r < rs.r || (ls.r == rs.r && ls.g < rs.g) || (ls.r == rs.r && ls.g == rs.g && ls.b <rs.b) ;
	}

    __forceinline RGBColor& operator += (RGBColor c)
    {
        r = Utility::BoundToByte(int(r) + int(c.r));
        g = Utility::BoundToByte(int(g) + int(c.g));
        b = Utility::BoundToByte(int(b) + int(c.b));
        a = Utility::BoundToByte(int(a) + int(c.a));
        return *this;
    }

    __forceinline static UINT DistL1(RGBColor a, RGBColor b)
    {
        return Math::Abs(int(a.r) - int(b.r)) +
               Math::Abs(int(a.g) - int(b.g)) +
               Math::Abs(int(a.b) - int(b.b));
    }
    __forceinline static UINT DistL2(RGBColor a, RGBColor b)
    {
        int DiffR = int(a.r) - int(b.r);
        int DiffG = int(a.g) - int(b.g);
        int DiffB = int(a.b) - int(b.b);
        return DiffR * DiffR + DiffG * DiffG + DiffB * DiffB;
    }

    __forceinline static RGBColor RandomColor()
    {
        return RGBColor(rand() & 255, rand() & 255, rand() & 255);
    }

    __forceinline RGBColor FlipBlueAndRed()
    {
        return RGBColor(b, g, r, a);
    }

    __forceinline RGBColor Inverse()
    {
        return RGBColor(255 - r, 255 - g, 255 - b, 255 - a);
    }

    RGBColor RGBToHSV() const;

    static RGBColor Interpolate(RGBColor LowColor, RGBColor HighColor, float s);
    static const RGBColor White;
    static const RGBColor Red;
    static const RGBColor Green;
    static const RGBColor Gray;
    static const RGBColor Blue;
    static const RGBColor Yellow;
    static const RGBColor Orange;
    static const RGBColor Magenta;
    static const RGBColor Black;
    static const RGBColor Cyan;
    static const RGBColor Purple;
	static const RGBColor DarkRed;  //1
	static const RGBColor LightViolet;
	static const RGBColor LightYellowOrange;
	static const RGBColor DarkerRedOrange;
	static const RGBColor PaleWarmBrown;
	static const RGBColor CMYKGreen;
	static const RGBColor LightCyan;

    BYTE r, g, b, a;
};

//
// Exact comparison functions.  Does not match the alpha channel.
//
__forceinline bool operator == (RGBColor left, RGBColor right)
{
    return ((left.r == right.r) && (left.g == right.g) && (left.b == right.b));
}

__forceinline bool operator != (RGBColor left, RGBColor right)
{
    return ((left.r != right.r) || (left.g != right.g) || (left.b != right.b));
}

#ifdef USE_D3D
//
// Helper function since D3DCOLORVALUE lacks a constructor
//
D3DCOLORVALUE MakeColor(float r, float g, float b, float a);
#endif
