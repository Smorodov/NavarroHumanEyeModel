#if defined(_MSC_VER)
#pragma once
#endif

#ifndef HUMAN_EYE_H
#define HUMAN_EYE_H
#include "gdt.h"
#include "vec.h"

#include "opencv2/opencv.hpp"
#include <string>



class Ray
{
public:
	// Ray Public Methods
	Ray() : mint(0.f), maxt(INFINITY), time(0.f), depth(0) { }
	Ray(const gdt::vec3f& origin, const gdt::vec3f& direction)
		: o(origin), d(direction), mint(0), maxt(INFINITY), time(0.f), depth(0) { }
	Ray(const gdt::vec3f& origin, const gdt::vec3f& direction,
		float start, float end = INFINITY, float t = 0.f, int d = 0)
		: o(origin), d(direction), mint(start), maxt(end), time(t), depth(d) { }
	Ray(const gdt::vec3f& origin, const gdt::vec3f& direction, const Ray& parent,
		float start, float end = INFINITY)
		: o(origin), d(direction), mint(start), maxt(end),
		time(parent.time), depth(parent.depth + 1) { }
	gdt::vec3f operator()(float t) const
	{
		return o + d * t;
	}


	// Ray Public Data
	gdt::vec3f o;
	gdt::vec3f d;
	mutable float mint, maxt;
	float time;
	int depth;
};

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


inline bool HasNaNs(gdt::vec3f& v)
{
	return isnan(v.x) || isnan(v.y) || isnan(v.z);
}

inline bool HasNaNs(Ray& r)
{
	return (HasNaNs(r.o) || HasNaNs(r.d) ||
		isnan(r.mint) || isnan(r.maxt));
}

class Lens
{
	public:

		inline void print(void)
		{
			std::cout << "-----------------" << std::endl;
			std::cout << "radius=" << radius << std::endl <<
				"asphericity=" << asphericity << std::endl <<
				"zPos=" << zPos << std::endl <<
				"refractionRatio=" << refractionRatio << std::endl <<
				"aperture=" << aperture << std::endl;
			std::cout << "-----------------" << std::endl;
		}

		Lens(float rad, float asph, float zpos, float refr, float aper);
		bool RefractRay(const Ray &inRay, Ray &outRay,float kn) const;
		gdt::vec3f GetNormal(const Ray * ray, const gdt::vec3f & p) const;
		bool OnLens(gdt::vec3f p) const;
		float yToZ(float);

		float radius;
		float asphericity;
		float zPos;
		float refractionRatio;
		float aperture;
};

class HumanEye
{
	public:
		HumanEye();
		~HumanEye();
		float GenerateRay(float focus, float x, float y, float z, float pupil, std::vector<gdt::vec3f>& hitPoints, float kn) const;
		void ParseSpecfile(const std::string & specfile);
		bool TraceLenses(const Ray &inRay, Ray &outRay, int start, int end,float kn) const;
		bool TraceLenses(const Ray &inRay, Ray &outRay, int start, int end, std::vector<gdt::vec3f>& points,float kn) const;
		
		gdt::vec3f getRandomPointOnLens() const;
		bool SetFocalLength(float focus);
		void SetDP(float);
		void SetPupilSize(float pupilSize)
		{
			lenses[3].aperture = pupilSize;
		}

	//private:
		// film plane's z position
		float filmPlane;
		float filmDistance;
		// lens disc z position
		float discZ;
		float apertureDiameter;
		float filmSize;
		std::vector<Lens> lenses;
		float dp;
		gdt::vec3f testPoint;
		float phi;
};

#endif
