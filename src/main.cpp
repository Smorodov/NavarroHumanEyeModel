#include "eye.h"

// https://peteroupc.github.io/colorgen.html
// waveLength = 650 - 250 / 270 * H


int main()
{
	float density = 2;

	float focus = 120;
	float srcZ = 520;
	float pupil = 6;
	cv::Mat src = cv::imread("F:/ImagesForTest/lena.jpg");
	cv::Mat screen = cv::Mat::zeros(1000, 1000, CV_32FC3);
	cv::Mat screenHitCount = cv::Mat::zeros(screen.size(), CV_32FC1);
	HumanEye camera;
	camera.SetFocalLength(focus);
	camera.SetPupilSize(pupil);
	
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
	float rad = 24;
	cx *= screen.cols / 2.0 / 24.0;
	cy *= screen.rows / 2.0 / 24.0;
	cx += screen.cols / 2.0;
	cy = screen.rows / 2.0 - cy;
	rad *= screen.cols / 2.0 / 24.0;
	
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
	cv::imwrite("res2.png", screen * 255);
	cv::waitKey();
	return 1;
}


/*
int main()
{
	float focus = 220;
	float srcZ = 120;
	float pupil = 6;
	cv::Mat screen = cv::Mat::zeros(1000, 1000, CV_32FC3);
	cv::Mat screenHitCount = cv::Mat::zeros(screen.size(), CV_32FC1);
	HumanEye camera;
	camera.SetFocalLength(focus);
	camera.SetPupilSize(pupil);
	
	for (auto l : camera.lenses)
	{
		l.print();
	}

	float cx = 0;
	float cy = 0;
	float rad = 24;
	cx += screen.cols / 2;
	cy += screen.rows / 2;
	rad *= screen.cols / 60.0;

			cv::Vec3b color(255, 255, 255);
			cv::Vec3f colorf((float)color[0] / 255.0, (float)color[1] / 255.0, (float)color[2] / 255.0);
			float R = 700;
			
			for (float r = 0; r < R; r += R/15)
			{
				for (float a = 0; a < 2 * M_PI; a += M_PI / 60)
				{
					float x = r * cos(a);
					float y = r * sin(a);
					camera.GenerateRay(focus, x, y, srcZ, pupil, screen, screenHitCount, colorf);
				}
			}
			
			float x = R * cos(0);
			float y = R * sin(0);
			camera.GenerateRay(focus, x, y, srcZ, pupil, screen, screenHitCount, colorf);
	

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
*/