#include "eye.h"

#include <fstream>
#include <algorithm>
#include <vector>


using namespace std;

HumanEye::HumanEye()
{
	// Функция чтения файла параметров линз и создание оптической системы.
	ParseSpecfile("data/Navarro.dat");
	float dz = -2.674404474;
	discZ = lenses[lenses.size() - 1].zPos + dz;
	testPoint = Point(0, 100, 200);
	phi = 0;
}

bool HumanEye::TraceLenses(const Ray& inRay, Ray& outRay, int start, int end) const
{
	Ray tempRay = inRay;
	
	if (inRay.o.z > 10)
	{
		// move start point close enough
		float t = (Point(0, 0, 10) - inRay.o).z / inRay.d.z;
		tempRay.o = inRay(t);
	}
	
	int increment = (end >= start) ? 1 : -1;
	int i = start;
	for (; i != end; i += increment)
	{
		if (!lenses[i].RefractRay(tempRay, outRay))
		{
			break;
		}
		tempRay = outRay;
	}
	if (i == end)
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool HumanEye::TraceLenses(const Ray& inRay, Ray& outRay, int start, int end, std::vector<Point>& points) const
{
	Ray tempRay = inRay;
	points.clear();
	points.push_back(inRay.o);
	if (inRay.o.z > 10)
	{
		float t = (Point(0, 0, 10) - inRay.o).z / inRay.d.z;
		tempRay.o = inRay(t);
	}
	int increment = (end >= start) ? 1 : -1;
	int i = start;
	for (; i != end; i += increment)
	{
		if (!lenses[i].RefractRay(tempRay, outRay))
		{
			points.push_back(outRay.o);
			break;
		}
		points.push_back(outRay.o);
		tempRay = outRay;
	}
	if (i == end)
	{
		//		points.push_back(outRay(100));
		return true;
	}
	else
	{
		return false;
	}
}

void ConcentricSampleDisk(float u1, float u2, float* dx, float* dy)
{
	float r, theta;
	// Map uniform random numbers to $[-1,1]^2$
	float sx = 2 * u1 - 1;
	float sy = 2 * u2 - 1;
	// Map square to $(r,\theta)$
	// Handle degeneracy at the origin
	if (sx == 0.0 && sy == 0.0)
	{
		*dx = 0.0;
		*dy = 0.0;
		return;
	}
	if (sx >= -sy)
	{
		if (sx > sy)
		{
			// Handle first region of disk
			r = sx;
			if (sy > 0.0)
			{
				theta = sy / r;
			}
			else
			{
				theta = 8.0f + sy / r;
			}
		}
		else
		{
			// Handle second region of disk
			r = sy;
			theta = 2.0f - sx / r;
		}
	}
	else
	{
		if (sx <= sy)
		{
			// Handle third region of disk
			r = -sx;
			theta = 4.0f - sy / r;
		}
		else
		{
			// Handle fourth region of disk
			r = -sy;
			theta = 6.0f + sx / r;
		}
	}
	theta *= M_PI / 4.f;
	*dx = r * cosf(theta);
	*dy = r * sinf(theta);
}

Point HumanEye::getRandomPointOnLens() const
{
	Point p;
	float ux, uy, dx, dy;
	ux = ((float)rand() / (RAND_MAX));
	uy = ((float)rand() / (RAND_MAX));
	ConcentricSampleDisk(ux, uy, &dx, &dy);
	p.x = dx * lenses[3].aperture / 2.0;
	p.y = dy * lenses[3].aperture / 2.0;
	p.z = discZ;
	return p;
}

void HumanEye::ParseSpecfile(const string& specfile)
{
	ifstream infile(specfile.c_str());
	if (!infile)
	{
		fprintf(stderr, "Cannot open spec file %s\n", specfile.c_str());
		exit(-1);
	}
	char line[512];
	float rad, zpos = 0, asph, axpos, lastRefr = 1.0f, refr, aper;
	while (!infile.eof())
	{
		infile.getline(line, 512);
		if (line[0] != '\0' && line[0] != '#' && line[0] != '\n')
		{
			sscanf(line, "%f\t%f\t%f\t%f\t%f\n", &rad, &asph, &axpos, &refr, &aper);
			Lens lens(-rad, asph, zpos, refr / lastRefr, aper);
			lenses.insert(lenses.begin(), lens);
			zpos -= axpos;
			if (refr != 0)
			{
				lastRefr = refr;
			}
		}
	}
}

HumanEye::~HumanEye()
{
	lenses.clear();
}

float HumanEye::GenerateRay(float focus, float x, float y, float z, float pupil, cv::Mat& screen, cv::Mat& screenHitCount, cv::Vec3f c) const
{
	Point objectSample = Point(x, y, z);
	std::vector<Point> points;
	
	// диаметр зрачка
	float aperture = lenses.back().aperture / 2.0;
	float nsteps = 10;// aperture* screen.cols / 2 / 30.0 / 10;
	float r2 = aperture * aperture / 4.0;
	float s = -aperture / 2.0;
	float e = aperture / 2.0;
	float st = 2.0*M_PI / nsteps;
			
	for (float i = s; i <= e; i += st)
	{
		for (float j = s; j <= e; j += st)
		{
			if (i * i + j * j > r2)
			{
				continue;
			}
			Point lensSample = Point(i, j, discZ);
			Ray filmLensRay(objectSample, Normalize(lensSample - objectSample), 0), primaryRay;
			if (TraceLenses(filmLensRay, primaryRay, (int)lenses.size() - 1, -1, points))
			{
				//primaryRay.d = Normalize(primaryRay.o - Point(0, 0, lenses[0].zPos + 2 * lenses[0].radius));
				//float t = (lenses[0].zPos - primaryRay.o.z) / primaryRay.d.z;
				//Point hit = primaryRay(t);
				Point hit = primaryRay.o;

				hit.x *= (float)screen.cols / 24.0; // assume diameter=30 mm
				hit.y *= (float)screen.rows / 24.0; // assume diameter=30 mm
				hit.x += (float)screen.cols / 2.0;
				hit.y += (float)screen.rows / 2.0;
				if (hit.x >= 0 && hit.x < screen.cols &&
				        hit.y >= 0 && hit.y < screen.rows
				   )
				{
					screen.at<cv::Vec3f>(round(hit.y), round(hit.x)) += c;
					screenHitCount.at<float>(round(hit.y), round(hit.x) ) += 1.0;
				}				
			}
		}
	}
	
	return 0;
}

void HumanEye::SetDP(float A)
{
	lenses[2].radius = -(10.2 - 1.75 * log(A + 1));
	lenses[1].radius = -(-6 + 0.2294 * log(A + 1));
	lenses[3].zPos = lenses[4].zPos - (3.05 - 0.05 * log(A + 1));
	lenses[2].zPos = lenses[3].zPos;
	lenses[1].zPos = lenses[2].zPos - (4 + 0.1 * log(A + 1));
	lenses[0].zPos = lenses[1].zPos - 16.3203;
	lenses[2].refractionRatio = (1.42 + 9e-5 * (10 * A + A * A)) / 1.3374;
	lenses[2].asphericity = -3.1316 - 0.34 * log(A + 1);
	lenses[1].asphericity = -1 - 0.125 * log(A + 1);
}

bool HumanEye::SetFocalLength(float f)
{
	Ray inRay, outRay;
	float start = 0, end = 13, current = start;
	int iter = 500;
	inRay.o = Point(0, 0, 0);
	while (iter--)
	{
		SetDP(current);
		inRay.o.z = f;
		float passed = 0;
		float focusSum = 0;
		int iter = 50;
		for (int i = 0; i < iter; ++i)
		{
			// trace 50 rays through eye model
			Point p = getRandomPointOnLens();
			p.x = 0;
			//p.y = 0.4;
			//p.z = 0;
			inRay.d = Normalize(p - inRay.o);
			//inRay.o = Point(0,1,1);
			//inRay.d = Vector(0,0,-1);
			if (TraceLenses(inRay, outRay, lenses.size() - 1, 0))
			{
				if (fabs(outRay.d.y) < 1e-6) continue;

				float t = -outRay.o.y / outRay.d.y;
				if (t > 0)
				{
					focusSum += outRay(t).z;
				}
				passed++;
			}
		}
		float focus = focusSum / passed;
		if (fabsf((float)(focus - lenses[0].zPos)) < 0.002)
		{
			//fprintf(fp, "%f", log(current));
			printf("%f %f\n", focus, current);
			return true;
		}
		if (focus > lenses[0].zPos)
		{
			end = current;
		}
		else
		{
			start = current;
		}
		current = (start + end) / 2;
	}
	printf("failed %f\n", f);
	return false;
}

// Add code for lens
Lens::Lens(float rad, float asph, float zpos, float refr, float aper)
	: radius(rad), asphericity(asph), zPos(zpos), refractionRatio(refr), aperture(aper)
{
}

Vector Lens::GetNormal(const Ray* ray, const Point& p) const
{
	Vector o = p - Point(0, 0, zPos);
	Vector result = Normalize(Vector(o.x,
	                                 o.y,
	                                 (1 + asphericity) * o.z - radius));
	if (result.z < 0)
	{
		result = -1 * result;
	}
	return result;
}

bool Lens::RefractRay(const Ray& inRay, Ray& outRay) const
{
	// точка попадантя луча в линзу
	Point phit;
	// -----------
	// Это апертура
	// ------------
	if (radius == 0)
	{
		// параметр точки на луче (расстояние от начала)
		float thit = (zPos - inRay.o.z) / inRay.d.z;
		// не попали в плоскость апертуры (луч в другую стоорону или недостаточеая длина луча)
		if (thit < 0 || thit > inRay.maxt)
		{
			return false;
		}
		// точка попадания луча в плоскость линзы
		phit = inRay(thit);
		// проверим попали ли в апертуру
		if (phit.x * phit.x + phit.y * phit.y > aperture * aperture * 0.25)
		{
			// не попали
			return false;
		}
		// новый луч из точки попадания,, параллельный сам себе.
		outRay.o = phit;
		outRay.d = inRay.d;
		return true;
	}
	// ---------
	// Это линза
	// // -------
	// Пересчитаем координаты начала луча относительно центра линзы
	Vector tOrigin = inRay.o - Point(0, 0, zPos);
	Vector tDir = inRay.d;//Normalize(inRay.d - Vector(0,0,zPos));
	// Compute quadratic sphere coefficients
	//float A = inRay.d.x * inRay.d.x + inRay.d.y * inRay.d.y + inRay.d.z * inRay.d.z;
	//float A = 1;
	float A = 1 + asphericity * tDir.z * tDir.z;
	float B = 2 * (tDir.x * tOrigin.x + tDir.y * tOrigin.y + (1 + asphericity) * tDir.z * tOrigin.z - radius * tDir.z);
	float C = tOrigin.x * tOrigin.x + tOrigin.y * tOrigin.y + (1 + asphericity) * tOrigin.z * tOrigin.z - 2 * radius * tOrigin.z;
	if (A == 0)
	{
		float t = -C / (B + 1e-6);
		phit = inRay(t);
	}
	else
	{
		const float discriminant = B * B - 4 * A * C;
		if (discriminant < 0)
		{
			return false;
		}
		const float t0 = (-B - sqrt(discriminant)) * 0.5f / A;
		const float t1 = (-B + sqrt(discriminant)) * 0.5f / A;
		//const float & tt = t0 > 0 ? t0 : t1;
		float tt;
		//if (asphericity <= -1) // first 2 lenses
		//	tt = t0 > t1? t0 : t1;
		//else
		//	tt = t0 > t1? t1 : t0;
		if (t0 > 0 && t1 > 0)
		{
			tt = t0 > t1 ? t1 : t0;
		}
		else if (t0 > 0 && t1 < 0)
		{
			tt = t0;
		}
		else if (t0 < 0 && t1 > 0)
		{
			tt = t1;
		}
		else if (t0 == 0)
		{
			tt = 0;
		}
		else
		{
			tt = -1;
		}
		if (tt < 0 || tt != tt)
		{
			return false;
		}
		phit = inRay(tt);
	}



	if (refractionRatio == 0)
	{
		// ------------------
		// Попали в сетчатку
		// ------------------
		outRay.o = phit;
		outRay.d = inRay.d;
		
		return true;
	}
	if (phit.HasNaNs() || !OnLens(phit))
	{
		return false;
	}
	// -----------------
	// нормальная линза
	// // -----------------
	//Vector n  = Normalize(phit - center);
	Vector n = GetNormal(&inRay, phit);
	Vector in = Normalize(tDir);
	float mu = refractionRatio;
	if (tDir.z < 0)
	{
		mu = 1.0 / mu;
	}
	float cosI = -Dot(n, in);
	float sinT2 = mu * mu * (1 - cosI * cosI);
	if (sinT2 > 1.0)
	{
		return false;
	}
	float cosT = sqrtf(1 - sinT2);
	outRay.d = mu * in + (mu * cosI - cosT) * n;
	outRay.o = phit;
	return true;
}

bool Lens::OnLens(Point p) const
{
	Vector pp = p - Point(0, 0, zPos);
	if (pp.x * pp.x + pp.y * pp.y > aperture * aperture * 0.25)
	{
		return false;
	}
	//if ((pp.z < -1e-6) != (radius < -1e-6))
	//	return false;
	return true;
}

float Lens::yToZ(float y)
{
	float posz1, posz2;
	if (refractionRatio == 0)
	{
		return zPos + radius - sqrt(radius * radius - y * y + 1e-6);
	}
	if (!Quadratic(1 + asphericity, -2 * radius, y * y, &posz1, &posz2))
	{
		return 1;
	}
	if (asphericity == -1)
	{
		return zPos + posz1;
	}
	if (asphericity < -1 && asphericity > -3)
	{
		return zPos + posz2;
	}
	if (asphericity < -3)
	{
		return zPos + posz1;
	}
	return zPos - ((radius > 0) ? posz1 : -posz2);
	//return zPos + posz1;
}

// Add code end