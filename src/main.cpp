#include "eye.h"

// https://peteroupc.github.io/colorgen.html
// waveLength = 650 - 250 / 270 * H


int main()
{
	float focus = 120;
	float srcZ = 120;
	float pupil = 3;
	cv::Mat src = cv::imread("F:/ImagesForTest/lena.jpg");
	cv::Mat screen = cv::Mat::zeros(1000, 1000, CV_32FC3);
	cv::Mat screenHitCount = cv::Mat::zeros(screen.size(), CV_32FC1);
	HumanEye camera;
	camera.SetFocalLength(focus);
	camera.SetPupilSize(pupil);
	float density = 4;
	cv::resize(src, src, cv::Size(), density, density);
	float xs = -float(src.cols) / density;
	float ys = -float(src.rows) / density;
	float xe = float(src.cols) / density;
	float ye = float(src.rows) / density;
	float dx = xe - xs;
	float dy = ye - ys;
	float stepx = dx / float(src.cols);
	float stepy = dy / float(src.rows);
	float cx = 0;
	float cy = 0;
	float rad = 29;
	cx *= screen.cols / 2.0 / 30.0;
	cy *= screen.rows / 2.0 / 30.0;
	cx += screen.cols / 2.0;
	cy = screen.rows / 2.0 - cy;
	rad *= screen.cols / 2.0 / 30.0;
	
	int row = 0;
	for (float srcY = ys; srcY < ye; srcY += stepy)
	{
		int col = 0;
		for (float srcX = xs; srcX < xe; srcX += stepx)
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
			cv::Vec3f colorf((float)color[0] / 255.0, (float)color[1] / 255.0, (float)color[2] / 255.0);
			camera.GenerateRay(focus, srcX, srcY, srcZ, pupil, screen, screenHitCount, colorf);
			++col;
		}
		++row;
	}
	std::vector<cv::Mat> ch;
	cv::split(screen, ch);
	for (auto& c : ch)
	{
		cv::divide(c, screenHitCount, c);
	}
	cv::merge(ch, screen);

	cv::circle(screen, cv::Point(cx, cy), rad, cv::Scalar(255, 255, 255), 2);
	cv::imshow("screen", screen);
	cv::imwrite("res1.png", screen * 255);
	cv::waitKey();
	return 1;
}


/*
int main()
{
	float focus = 2120;
	float srcZ = 120;
	float pupil = 3;
	cv::Mat screen = cv::Mat::zeros(1000, 1000, CV_32FC3);
	cv::Mat screenHitCount = cv::Mat::zeros(screen.size(), CV_32FC1);
	HumanEye camera;
	camera.SetFocalLength(focus);
	camera.SetPupilSize(pupil);

	float cx = 0;
	float cy = 0;
	float rad = 29;
	cx *= screen.cols / 2 / 30.0;
	cy *= screen.rows / 2 / 30.0;
	cx += screen.cols / 2;
	cy = screen.rows / 2 - cy;
	rad *= screen.cols / 2 / 30.0;

			cv::Vec3b color(255, 255, 255);
			cv::Vec3f colorf((float)color[0] / 255.0, (float)color[1] / 255.0, (float)color[2] / 255.0);
			camera.GenerateRay(focus, 0, 0, srcZ, pupil, screen, screenHitCount, colorf);

	cv::circle(screen, cv::Point(cx, cy), rad, cv::Scalar(255, 255, 255), 2);

	std::vector<cv::Mat> ch;
	cv::split(screen, ch);
	for (auto& c : ch)
	{
		cv::divide(c, screenHitCount, c);
	}
	cv::merge(ch, screen);
	cv::imshow("screen", screen);
	cv::imwrite("res1.png", screen * 255);
	cv::waitKey();
	return 1;
}*/