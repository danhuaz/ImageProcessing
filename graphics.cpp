#include "graphics.h"

Point operator+ (const Point& p, const Point& q) {
	Point plus;
	plus.x = p.x + q.x;
	plus.y = p.y + q.y;
	plus.z = p.z + q.z;
	return plus;
}

Point operator- (const Point& p, const Point& q) {
	Point plus;
	plus.x = p.x - q.x;
	plus.y = p.y - q.y;
	plus.z = p.z - q.z;
	return plus;
}

Point operator* (const float& a, const Point& p) {
	Point q;
	q.x = p.x*a;
	q.y = p.y*a;
	q.z = p.z*a;
	return q;
}

Point:: Point(const Point& p) {
	x = p.x;
	y = p.y;
	z = p.z;
}

float Point::scale() {
	return(sqrt(x*x + y*y + z*z));
}

void Point:: unit() {
	float s = this->scale();
	x = x / s;
	y = y / s;
	z = z / s;
}

float Point::Distance(Point p) {
	return sqrt((x - p.x)*(x - p.x) + (y - p.y)*(y - p.y) + (z - p.z)*(z - p.z));
}

Point Point::minus() {
	Point p;
	p.x = -x;
	p.y = -y;
	p.z = -z;
	return p;
}

float dot(Point a, Point b) {
	return (a.x*b.x + a.y*b.y + a.z*b.z);
}

Point crossproduct(Point a, Point b) {
	Point p;
	p.x = a.y*b.z - a.z*b.y;
	p.y = a.z*b.x - a.x*b.z;
	p.z = a.x*b.y - a.y*b.x;
	return p;
}

void FloatPixel::Set(float r_, float g_, float b_) {
	r = r_;
	g = g_;
	b = b_;
}


Pixel FloatPixel::ToPixel() {
	Pixel p;
	p.SetClamp(r, g, b);
	return p;
}

FloatPixel operator+ (const FloatPixel& p, const FloatPixel& q) {
	FloatPixel plus;
	plus.r = q.r + p.r;
	plus.g = q.g + p.g;
	plus.b = q.b + p.b;
	return plus;
}

FloatPixel operator* (const FloatPixel& p, const FloatPixel& q) {
	FloatPixel time;
	time.r = q.r*p.r;
	time.g = q.g*p.g;
	time.b = q.b*p.b;
	return time;
}

FloatPixel operator* (const float& a, const FloatPixel& p) {
	FloatPixel q;
	q.r = a*p.r;
	q.g = a*p.g;
	q.b = a*p.b;
	return q;
}

FloatPixel operator* (const FloatPixel& p, const float& a) {
	FloatPixel q;
	q.r = a*p.r;
	q.g = a*p.g;
	q.b = a*p.b;
	return q;
}

float Sphere:: Discriminant(Point ray_ori, Point ray_dire) {
	Point p = ray_ori - center;
	float discriminant = dot(ray_dire, p)*dot(ray_dire, p) - (dot(ray_dire, ray_dire)*(dot(p, p) - radius*radius));
	return discriminant;
}

bool Sphere::IfhitSphere(Point ray_ori, Point ray_dire) {
	if (Discriminant(ray_ori, ray_dire) >= 0) {
		return true;
	}
	else {
		return false;
	}
}

float Sphere::Hitpoint(Point ray_ori, Point ray_dire) {
	Point p = ray_ori - center;
	float discriminant = Discriminant(ray_ori, ray_dire);
	if (discriminant >= 0) {
		float b = dot(ray_dire, p);
		float a = dot(ray_dire, ray_dire);
		float t1 = (-b + sqrt(discriminant)) / a;
		float t2 = (-b - sqrt(discriminant)) / a;
		return fminf(t1, t2);
	}
}

