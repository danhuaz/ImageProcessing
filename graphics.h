#include <SDL.h>
#include <SDL_opengl.h>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <math.h>
#include "image.h"


#define MAX_CHARACTER 1024
#define Pi 3.141592657535897932

using namespace std; 

class Point {
public:
	float x;
	float y;
	float z;

	Point minus();
	void unit(); //unitize the vector
	float scale();
	Point(float x_ = 0, float y_ = 0, float z_ = 0) : x(x_), y(y_), z(z_) {};
	Point(const Point& p);
	float Distance(Point p);
};
Point operator+ (const Point& p, const Point& q);
Point operator- (const Point& p, const Point& q);
Point operator* (const float& a, const Point& p);
Point crossproduct(Point a, Point b);
float dot(Point a, Point b);

class FloatPixel {
public:
	float r, g, b;
	FloatPixel(float r_ = 0, float g_ = 0, float b_ = 0) : r(r_), g(g_), b(b_) {}
	void Set(float r, float g, float b);
	Pixel ToPixel();
};
FloatPixel operator+ (const FloatPixel& p, const FloatPixel& q);
FloatPixel operator* (const FloatPixel& p, const FloatPixel& q);
FloatPixel operator* (const float& a, const FloatPixel& p);
FloatPixel operator* (const FloatPixel& p, const float& a);

class Sphere {
public:
	float radius;
	Point center;
	int material;
	Sphere(float radius_ = 0, Point center_ = 0, int material_ = 0): radius(radius_), center(center_), material(material_) {};
	float Discriminant(Point ray_ori, Point ray_dire);
	bool IfhitSphere(Point ray_ori, Point ray_dire);
	float Hitpoint(Point ray_ori, Point ray_dire);
	void Setmaterial(int i) { material = i; };
};

class Directional_light {
public:
	FloatPixel color;
	Point direction;
};

class Point_light {
public:
	FloatPixel color;
	Point position;
};

class Spot_light {
public:
	FloatPixel color;
	Point location;
	Point direction;
	float angle1, angle2;
};

class Material {
public:
	FloatPixel ambient;
	FloatPixel diffuse;
	FloatPixel specular;
	FloatPixel transmissive;
	float ns;
	float ior;
	//Material(FloatPixel ambient_ = 0, FloatPixel diffuse_ = 1, FloatPixel specular_ = 0, FloatPixel transmissive_ = 0, \
	//	float ns_ = 5, float ior_ = 1) : ambient(ambient_), diffuse(diffuse_), specular(specular_), transmissive(transmissive_), \
	//	ns(ns_), ior(ior_) {};
	Material() {
		ambient = 0;
		diffuse.Set(1,1,1);
		specular = 0;
		transmissive = 0;
		ns = 5;
		ior = 1;
	}

};
