#pragma once

#include <vector>
#include <assert.h>

class Vector;
class Point;
class Normal;

inline float Clamp(float val, float low, float high)
{
	if (val < low)
	{
		return low;
	}
	else if (val > high)
	{
		return high;
	}
	else
	{
		return val;
	}
}


inline int Clamp(int val, int low, int high)
{
	if (val < low)
	{
		return low;
	}
	else if (val > high)
	{
		return high;
	}
	else
	{
		return val;
	}
}

inline bool Quadratic(float A, float B, float C, float* t0, float* t1)
{
	// Find quadratic discriminant
	float discrim = B * B - 4.f * A * C;
	if (discrim <= 0.)
	{
		return false;
	}
	float rootDiscrim = sqrt(discrim);
	// Compute quadratic _t_ values
	float q;
	if (B < 0)
	{
		q = -.5f * (B - rootDiscrim);
	}
	else
	{
		q = -.5f * (B + rootDiscrim);
	}
	*t0 = q / A;
	*t1 = C / q;
	//*t0 = (-B - rootDiscrim) / 2 / A;
	//*t1 = (-B + rootDiscrim) / 2 / A;
	if (*t0 > *t1)
	{
		std::swap(*t0, *t1);
	}
	return true;
}

// Geometry Declarations
class Vector
{
	public:
		// Vector Public Methods
		Vector()
		{
			x = y = z = 0.f;
		}
		Vector(float xx, float yy, float zz)
			: x(xx), y(yy), z(zz)
		{
			assert(!HasNaNs());
		}
		bool HasNaNs() const
		{
			return isnan(x) || isnan(y) || isnan(z);
		}
		explicit Vector(const Point &p);
#ifndef NDEBUG
		// The default versions of these are fine for release builds; for debug
		// we define them so that we can add the assert checks.
		Vector(const Vector &v)
		{
			assert(!v.HasNaNs());
			x = v.x;
			y = v.y;
			z = v.z;
		}

		Vector &operator=(const Vector &v)
		{
			assert(!v.HasNaNs());
			x = v.x;
			y = v.y;
			z = v.z;
			return *this;
		}
#endif // !NDEBUG
		Vector operator+(const Vector &v) const
		{
			assert(!v.HasNaNs());
			return Vector(x + v.x, y + v.y, z + v.z);
		}

		Vector& operator+=(const Vector &v)
		{
			assert(!v.HasNaNs());
			x += v.x;
			y += v.y;
			z += v.z;
			return *this;
		}
		Vector operator-(const Vector &v) const
		{
			assert(!v.HasNaNs());
			return Vector(x - v.x, y - v.y, z - v.z);
		}

		Vector& operator-=(const Vector &v)
		{
			assert(!v.HasNaNs());
			x -= v.x;
			y -= v.y;
			z -= v.z;
			return *this;
		}
		Vector operator*(float f) const
		{
			return Vector(f*x, f*y, f*z);
		}

		Vector &operator*=(float f)
		{
			assert(!isnan(f));
			x *= f;
			y *= f;
			z *= f;
			return *this;
		}
		Vector operator/(float f) const
		{
			assert(f != 0);
			float inv = 1.f / f;
			return Vector(x * inv, y * inv, z * inv);
		}

		Vector &operator/=(float f)
		{
			assert(f != 0);
			float inv = 1.f / f;
			x *= inv;
			y *= inv;
			z *= inv;
			return *this;
		}
		Vector operator-() const
		{
			return Vector(-x, -y, -z);
		}
		float operator[](int i) const
		{
			assert(i >= 0 && i <= 2);
			return (&x)[i];
		}

		float &operator[](int i)
		{
			assert(i >= 0 && i <= 2);
			return (&x)[i];
		}
		float LengthSquared() const
		{
			return x*x + y*y + z*z;
		}
		float Length() const
		{
			return sqrtf(LengthSquared());
		}
		explicit Vector(const Normal &n);

		bool operator==(const Vector &v) const
		{
			return x == v.x && y == v.y && z == v.z;
		}
		bool operator!=(const Vector &v) const
		{
			return x != v.x || y != v.y || z != v.z;
		}

		// Vector Public Data
		float x, y, z;
};


class Point
{
	public:
		// Point Public Methods
		Point()
		{
			x = y = z = 0.f;
		}
		Point(float xx, float yy, float zz)
			: x(xx), y(yy), z(zz)
		{
			assert(!HasNaNs());
		}
#ifndef NDEBUG
		Point(const Point &p)
		{
			assert(!p.HasNaNs());
			x = p.x;
			y = p.y;
			z = p.z;
		}

		Point &operator=(const Point &p)
		{
			assert(!p.HasNaNs());
			x = p.x;
			y = p.y;
			z = p.z;
			return *this;
		}
#endif // !NDEBUG
		Point operator+(const Vector &v) const
		{
			assert(!v.HasNaNs());
			return Point(x + v.x, y + v.y, z + v.z);
		}

		Point &operator+=(const Vector &v)
		{
			assert(!v.HasNaNs());
			x += v.x;
			y += v.y;
			z += v.z;
			return *this;
		}
		Vector operator-(const Point &p) const
		{
			assert(!p.HasNaNs());
			return Vector(x - p.x, y - p.y, z - p.z);
		}

		Point operator-(const Vector &v) const
		{
			assert(!v.HasNaNs());
			return Point(x - v.x, y - v.y, z - v.z);
		}

		Point &operator-=(const Vector &v)
		{
			assert(!v.HasNaNs());
			x -= v.x;
			y -= v.y;
			z -= v.z;
			return *this;
		}
		Point &operator+=(const Point &p)
		{
			assert(!p.HasNaNs());
			x += p.x;
			y += p.y;
			z += p.z;
			return *this;
		}
		Point operator+(const Point &p) const
		{
			assert(!p.HasNaNs());
			return Point(x + p.x, y + p.y, z + p.z);
		}
		Point operator* (float f) const
		{
			return Point(f*x, f*y, f*z);
		}
		Point &operator*=(float f)
		{
			x *= f;
			y *= f;
			z *= f;
			return *this;
		}
		Point operator/ (float f) const
		{
			float inv = 1.f/f;
			return Point(inv*x, inv*y, inv*z);
		}
		Point &operator/=(float f)
		{
			float inv = 1.f/f;
			x *= inv;
			y *= inv;
			z *= inv;
			return *this;
		}
		float operator[](int i) const
		{
			assert(i >= 0 && i <= 2);
			return (&x)[i];
		}

		float &operator[](int i)
		{
			assert(i >= 0 && i <= 2);
			return (&x)[i];
		}
		bool HasNaNs() const
		{
			return isnan(x) || isnan(y) || isnan(z);
		}

		bool operator==(const Point &p) const
		{
			return x == p.x && y == p.y && z == p.z;
		}
		bool operator!=(const Point &p) const
		{
			return x != p.x || y != p.y || z != p.z;
		}

		// Point Public Data
		float x, y, z;
};


class Normal
{
	public:
		// Normal Public Methods
		Normal()
		{
			x = y = z = 0.f;
		}
		Normal(float xx, float yy, float zz)
			: x(xx), y(yy), z(zz)
		{
			assert(!HasNaNs());
		}
		Normal operator-() const
		{
			return Normal(-x, -y, -z);
		}
		Normal operator+ (const Normal &n) const
		{
			assert(!n.HasNaNs());
			return Normal(x + n.x, y + n.y, z + n.z);
		}

		Normal& operator+=(const Normal &n)
		{
			assert(!n.HasNaNs());
			x += n.x;
			y += n.y;
			z += n.z;
			return *this;
		}
		Normal operator- (const Normal &n) const
		{
			assert(!n.HasNaNs());
			return Normal(x - n.x, y - n.y, z - n.z);
		}

		Normal& operator-=(const Normal &n)
		{
			assert(!n.HasNaNs());
			x -= n.x;
			y -= n.y;
			z -= n.z;
			return *this;
		}
		bool HasNaNs() const
		{
			return isnan(x) || isnan(y) || isnan(z);
		}
		Normal operator*(float f) const
		{
			return Normal(f*x, f*y, f*z);
		}

		Normal &operator*=(float f)
		{
			x *= f;
			y *= f;
			z *= f;
			return *this;
		}
		Normal operator/(float f) const
		{
			assert(f != 0);
			float inv = 1.f/f;
			return Normal(x * inv, y * inv, z * inv);
		}

		Normal &operator/=(float f)
		{
			assert(f != 0);
			float inv = 1.f/f;
			x *= inv;
			y *= inv;
			z *= inv;
			return *this;
		}
		float LengthSquared() const
		{
			return x*x + y*y + z*z;
		}
		float Length() const
		{
			return sqrtf(LengthSquared());
		}

#ifndef NDEBUG
		Normal(const Normal &n)
		{
			assert(!n.HasNaNs());
			x = n.x;
			y = n.y;
			z = n.z;
		}

		Normal &operator=(const Normal &n)
		{
			assert(!n.HasNaNs());
			x = n.x;
			y = n.y;
			z = n.z;
			return *this;
		}
#endif // !NDEBUG
		explicit Normal(const Vector &v)
			: x(v.x), y(v.y), z(v.z)
		{
			assert(!v.HasNaNs());
		}
		float operator[](int i) const
		{
			assert(i >= 0 && i <= 2);
			return (&x)[i];
		}

		float &operator[](int i)
		{
			assert(i >= 0 && i <= 2);
			return (&x)[i];
		}

		bool operator==(const Normal &n) const
		{
			return x == n.x && y == n.y && z == n.z;
		}
		bool operator!=(const Normal &n) const
		{
			return x != n.x || y != n.y || z != n.z;
		}

		// Normal Public Data
		float x, y, z;
};


class Ray
{
	public:
		// Ray Public Methods
		Ray() : mint(0.f), maxt(INFINITY), time(0.f), depth(0) { }
		Ray(const Point &origin, const Vector &direction)
			: o(origin), d(direction), mint(0), maxt(INFINITY), time(0.f), depth(0) { }
		Ray(const Point &origin, const Vector &direction,
		    float start, float end = INFINITY, float t = 0.f, int d = 0)
			: o(origin), d(direction), mint(start), maxt(end), time(t), depth(d) { }
		Ray(const Point &origin, const Vector &direction, const Ray &parent,
		    float start, float end = INFINITY)
			: o(origin), d(direction), mint(start), maxt(end),
			  time(parent.time), depth(parent.depth+1) { }
		Point operator()(float t) const
		{
			return o + d * t;
		}
		bool HasNaNs() const
		{
			return (o.HasNaNs() || d.HasNaNs() ||
			        isnan(mint) || isnan(maxt));
		}

		// Ray Public Data
		Point o;
		Vector d;
		mutable float mint, maxt;
		float time;
		int depth;
};



// Geometry Inline Functions
inline Vector::Vector(const Point &p)
	: x(p.x), y(p.y), z(p.z)
{
	assert(!HasNaNs());
}


inline Vector operator*(float f, const Vector &v)
{
	return v*f;
}
inline float Dot(const Vector &v1, const Vector &v2)
{
	assert(!v1.HasNaNs() && !v2.HasNaNs());
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}


inline float AbsDot(const Vector &v1, const Vector &v2)
{
	assert(!v1.HasNaNs() && !v2.HasNaNs());
	return fabsf(Dot(v1, v2));
}


inline Vector Cross(const Vector &v1, const Vector &v2)
{
	assert(!v1.HasNaNs() && !v2.HasNaNs());
	float v1x = v1.x, v1y = v1.y, v1z = v1.z;
	float v2x = v2.x, v2y = v2.y, v2z = v2.z;
	return Vector((v1y * v2z) - (v1z * v2y),
	              (v1z * v2x) - (v1x * v2z),
	              (v1x * v2y) - (v1y * v2x));
}


inline Vector Cross(const Vector &v1, const Normal &v2)
{
	assert(!v1.HasNaNs() && !v2.HasNaNs());
	float v1x = v1.x, v1y = v1.y, v1z = v1.z;
	float v2x = v2.x, v2y = v2.y, v2z = v2.z;
	return Vector((v1y * v2z) - (v1z * v2y),
	              (v1z * v2x) - (v1x * v2z),
	              (v1x * v2y) - (v1y * v2x));
}


inline Vector Cross(const Normal &v1, const Vector &v2)
{
	assert(!v1.HasNaNs() && !v2.HasNaNs());
	float v1x = v1.x, v1y = v1.y, v1z = v1.z;
	float v2x = v2.x, v2y = v2.y, v2z = v2.z;
	return Vector((v1y * v2z) - (v1z * v2y),
	              (v1z * v2x) - (v1x * v2z),
	              (v1x * v2y) - (v1y * v2x));
}


inline Vector Normalize(const Vector &v)
{
	return v / v.Length();
}
inline void CoordinateSystem(const Vector &v1, Vector *v2, Vector *v3)
{
	if (fabsf(v1.x) > fabsf(v1.y))
	{
		float invLen = 1.f / sqrtf(v1.x*v1.x + v1.z*v1.z);
		*v2 = Vector(-v1.z * invLen, 0.f, v1.x * invLen);
	}
	else
	{
		float invLen = 1.f / sqrtf(v1.y*v1.y + v1.z*v1.z);
		*v2 = Vector(0.f, v1.z * invLen, -v1.y * invLen);
	}
	*v3 = Cross(v1, *v2);
}


inline float Distance(const Point &p1, const Point &p2)
{
	return (p1 - p2).Length();
}


inline float DistanceSquared(const Point &p1, const Point &p2)
{
	return (p1 - p2).LengthSquared();
}


inline Point operator*(float f, const Point &p)
{
	assert(!p.HasNaNs());
	return p*f;
}


inline Normal operator*(float f, const Normal &n)
{
	return Normal(f*n.x, f*n.y, f*n.z);
}


inline Normal Normalize(const Normal &n)
{
	return n / n.Length();
}


inline Vector::Vector(const Normal &n)
	: x(n.x), y(n.y), z(n.z)
{
	assert(!n.HasNaNs());
}


inline float Dot(const Normal &n1, const Vector &v2)
{
	assert(!n1.HasNaNs() && !v2.HasNaNs());
	return n1.x * v2.x + n1.y * v2.y + n1.z * v2.z;
}


inline float Dot(const Vector &v1, const Normal &n2)
{
	assert(!v1.HasNaNs() && !n2.HasNaNs());
	return v1.x * n2.x + v1.y * n2.y + v1.z * n2.z;
}


inline float Dot(const Normal &n1, const Normal &n2)
{
	assert(!n1.HasNaNs() && !n2.HasNaNs());
	return n1.x * n2.x + n1.y * n2.y + n1.z * n2.z;
}


inline float AbsDot(const Normal &n1, const Vector &v2)
{
	assert(!n1.HasNaNs() && !v2.HasNaNs());
	return fabsf(n1.x * v2.x + n1.y * v2.y + n1.z * v2.z);
}


inline float AbsDot(const Vector &v1, const Normal &n2)
{
	assert(!v1.HasNaNs() && !n2.HasNaNs());
	return fabsf(v1.x * n2.x + v1.y * n2.y + v1.z * n2.z);
}


inline float AbsDot(const Normal &n1, const Normal &n2)
{
	assert(!n1.HasNaNs() && !n2.HasNaNs());
	return fabsf(n1.x * n2.x + n1.y * n2.y + n1.z * n2.z);
}


inline Normal Faceforward(const Normal &n, const Vector &v)
{
	return (Dot(n, v) < 0.f) ? -n : n;
}


inline Normal Faceforward(const Normal &n, const Normal &n2)
{
	return (Dot(n, n2) < 0.f) ? -n : n;
}



inline Vector Faceforward(const Vector &v, const Vector &v2)
{
	return (Dot(v, v2) < 0.f) ? -v : v;
}



inline Vector Faceforward(const Vector &v, const Normal &n2)
{
	return (Dot(v, n2) < 0.f) ? -v : v;
}


inline Vector SphericalDirection(float sintheta,
                                 float costheta, float phi)
{
	return Vector(sintheta * cosf(phi),
	              sintheta * sinf(phi),
	              costheta);
}


inline Vector SphericalDirection(float sintheta, float costheta,
                                 float phi, const Vector &x,
                                 const Vector &y, const Vector &z)
{
	return sintheta * cosf(phi) * x +
	       sintheta * sinf(phi) * y + costheta * z;
}


inline float SphericalTheta(const Vector &v)
{
	return acosf(Clamp(v.z, -1.f, 1.f));
}


inline float SphericalPhi(const Vector &v)
{
	float p = atan2f(v.y, v.x);
	return (p < 0.f) ? p + 2.f*M_PI : p;
}

