#if defined(_MSC_VER)
#pragma once
#endif

#ifndef HUMAN_EYE_H
#define HUMAN_EYE_H
#include "geometry.h"
#include "opencv2/opencv.hpp"
#include <string>
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
		Vector GetNormal(const Ray * ray, const Point & p) const;
		bool OnLens(Point p) const;
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
		float GenerateRay(float focus, float x, float y, float z, float pupil, std::vector<Point>& hitPoints, float kn) const;
		void ParseSpecfile(const std::string & specfile);
		bool TraceLenses(const Ray &inRay, Ray &outRay, int start, int end,float kn) const;
		bool TraceLenses(const Ray &inRay, Ray &outRay, int start, int end, std::vector<Point>& points,float kn) const;
		Point getRandomPointOnLens() const;
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
		Point testPoint;
		float phi;
};

#endif
