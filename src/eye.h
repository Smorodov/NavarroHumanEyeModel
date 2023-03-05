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
		Lens(float rad, float asph, float zpos, float refr, float aper);
		bool RefractRay(const Ray &inRay, Ray &outRay) const;
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
		float GenerateRay(float focus, float x, float y, float z, float pupil, cv::Mat& screen,cv::Mat& screenHitCount,cv::Vec3f c) const;
		void ParseSpecfile(const std::string & specfile);
		bool TraceLenses(const Ray &inRay, Ray &outRay, int start, int end) const;
		bool TraceLenses(const Ray &inRay, Ray &outRay, int start, int end, std::vector<Point>& points) const;
		Point getRandomPointOnLens() const;
		bool SetFocalLength(float focus);
		void SetDP(float);
		void SetPupilSize(float pupilSize)
		{
			lenses[3].aperture = pupilSize;
		}

	private:
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