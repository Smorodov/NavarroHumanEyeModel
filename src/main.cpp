#include "eye.h"

// https://peteroupc.github.io/colorgen.html
// waveLength = 650 - 250 / 270 * H

#if 1
int main()
{
	float density = 2;

	float focus = 320;
	float srcZ = 320;
	float pupil = 6;
	cv::Mat src = cv::imread("F:/ImagesForTest/lena.jpg");
	cv::Mat screen = cv::Mat::zeros(1000, 1000, CV_32FC3);
	cv::Mat screenHitCount = cv::Mat::zeros(screen.size(), CV_32FC3);
	HumanEye camera;
	camera.SetFocalLength(focus);
	camera.SetPupilSize(pupil);

	cv::resize(src, src, cv::Size(), density, density,cv::INTER_NEAREST);
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

	float kns[3] = { 0.99,1,1.01 };

	for (int channel = 0; channel < 3; channel++)
	{

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
					color[channel] = src.at<cv::Vec3b>(row, col)[channel];
				}
				cv::Vec3f colorf((float)color[0] / 255.0, (float)color[1] / 255.0, (float)color[2] / 255.0);

				std::vector<Point> hitPoints;
				camera.GenerateRay(focus, srcX, srcY, srcZ, pupil, hitPoints, kns[channel]);
				float np = hitPoints.size();
				float dx = (float)screen.cols / 2.0;
				float dy = (float)screen.rows / 2.0;
				float scale = (float)screen.cols / 24;

				for (auto& p : hitPoints)
				{
					// учет наклона стенки сетчатки
					float ang = 90.0 * sqrt(p.x * p.x + p.y * p.y) / 12.0;
					p.x = round(p.x * scale + dx);
					p.y = round(p.y * scale + dy);
					if (p.x >= 0 && p.x < screen.cols && p.y >= 0 && p.y < screen.rows)
					{
						screen.at<cv::Vec3f>(p.y, p.x) += colorf;// * cos(ang*M_PI/180.0);
						screenHitCount.at<cv::Vec3f>(p.y, p.x)[channel] += 1.0;
					}
				}
				++col;
			}
			++row;
		}
	}
	// не знаю почему 2.0, видимо линзы собирают местами больше
	//cv::normalize(screen, screen,0,1.5, cv::NORM_MINMAX);
	std::vector<cv::Mat> ch;
	cv::split(screen, ch);
	std::vector<cv::Mat> ch2;
	cv::split(screenHitCount, ch2);
	int ii = 0;
	for (auto& c : ch)
	{
		cv::divide(c, ch2[ii], c);
		++ii;
	}
	cv::merge(ch, screen);


	cv::circle(screen, cv::Point(cx, cy), rad, cv::Scalar(255, 255, 255), 2);
	cv::imshow("screen", screen);
	cv::imwrite("res2.png", screen * 255);
	cv::waitKey();
	return 1;
}
#else
int main()
{
	float focus = 220;
	float srcZ = 120;
	float pupil = 6;
	cv::Mat screen = cv::Mat::zeros(1000, 1000, CV_32FC3);
	cv::Mat screenHitCount = cv::Mat::zeros(screen.size(), CV_32FC3);
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
	rad *= screen.cols / 24.0;
	
	float R = 200;
	std::vector<Point> hitPoints;
	
	float kns[3] = { 0.997,1,1.003 };

	for (int channel = 0; channel < 3; channel++)
	{
		for (float r = 0; r < R * R; r += R * R / 10)
		{
			for (float a = 0; a < 2 * M_PI; a += M_PI / 12)
			{
				float x = sqrt(r) * cos(a);
				float y = sqrt(r) * sin(a);

				camera.GenerateRay(focus, x, y, srcZ, pupil, hitPoints, kns[channel]);
				if (hitPoints.size() == 0) continue;

				float np = hitPoints.size();
				float dx = (float)screen.cols / 2.0;
				float dy = (float)screen.rows / 2.0;
				float scale = (float)screen.cols / 24;
				

				cv::Vec3f colorf = cv::Vec3f(0, 0, 0);
				colorf[channel] = 1;
				

				for (auto& p : hitPoints)
				{
					p.x = round(p.x * scale + dx);
					p.y = round(p.y * scale + dy);
					if (p.x >= 0 && p.x < screen.cols && p.y >= 0 && p.y < screen.rows)
					{
						//screen.at<cv::Vec3f>(p.y, p.x) = cv::Vec3f(1, 1, 1);
						screen.at<cv::Vec3f>(p.y, p.x) += colorf;
						screenHitCount.at<cv::Vec3f>(p.y, p.x) += colorf;
					}
				}
			}
		}
	}

	std::vector<cv::Mat> ch,ch2;
	cv::split(screen, ch);
	cv::split(screenHitCount, ch2);
	int ii = 0;
	for (auto& c : ch)
	{
		cv::divide(c, ch2[ii], c);
		++ii;
	}
	cv::merge(ch, screen);

	cv::circle(screen, cv::Point(cx, cy), rad, cv::Scalar(255, 255, 255), 2);
	cv::imshow("screen", screen);
	cv::imwrite("res1.png", screen * 255);
	cv::waitKey();
	return 1;
}
#endif
