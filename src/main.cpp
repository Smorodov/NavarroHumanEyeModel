#include "eye.h"
#include <direct.h>
#include "opencv2/opencv.hpp"

// https://peteroupc.github.io/colorgen.html
// waveLength = 650 - 250 / 270 * H

int main()
{
	float focus = 2120;
	float srcZ = 120;
	float pupil = 3;
	cv::Mat src = cv::imread("F:/ImagesForTest/lena.jpg");
	cv::Mat screen = cv::Mat::zeros(1000, 1000, CV_32FC3);
	cv::Mat screenHitCount = cv::Mat::zeros(screen.size(), CV_32FC1);
	HumanEye camera;
	camera.SetFocalLength(focus);
	camera.SetPupilSize(pupil);
	float density = 1;
	cv::resize(src, src, cv::Size(), density, density);
	float xs = -float(src.cols ) / density;
	float ys = -float(src.rows) / density;
	float xe = float(src.cols ) / density;
	float ye = float(src.rows) / density;
	float dx = xe - xs;
	float dy = ye - ys;
	float stepx = dx / float(src.cols);
	float stepy = dy / float(src.rows);
	float cx = 0;
	float cy = 0;
	float rad = 28;
	cx *= screen.cols / 2 / 30.0;
	cy *= screen.rows / 2 / 30.0;
	cx += screen.cols / 2;
	cy = screen.rows / 2 - cy;
	rad *= screen.cols / 2 / 30.0;
	// just sample 1/4 circle is enough (symmetry)
	int row = 0;
	for (float srcY = ys; srcY < ye; srcY+=stepy)
	{
		int col = 0;
		for (float srcX = xs; srcX < xe; srcX+=stepx)
		{
			cv::Vec3b color(0, 0, 0);
			if (row >= 0 &&
			        row < src.rows &&
			        col >= 0 &&
			        col < src.cols
			   )
			{
				color = src.at<cv::Vec3b>(row, col);
			}
			cv::Vec3f colorf( (float)color[0]/255.0, (float)color[1]/255.0, (float)color[2]/255.0);
			camera.GenerateRay(focus, srcX, srcY, srcZ, pupil, screen, screenHitCount, colorf);
			++col;
		}
		++row;
	}
	cv::circle(screen, cv::Point(cx, cy), rad, cv::Scalar(255, 255, 255), 2);
	//cv::normalize(screen,screen,0,1,cv::NORM_MINMAX);
	std::vector<cv::Mat> ch;
	cv::split(screen, ch);
for (auto& c : ch)
	{
		cv::divide(c, screenHitCount, c);
	}
	cv::merge(ch, screen);
	cv::imshow("screen", screen);
	cv::imwrite("res1.png",screen*255);
	cv::waitKey();
	return 1;
}