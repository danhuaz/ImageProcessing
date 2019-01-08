#include "image.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <time.h>
#define Pi 3.141592657535897932

/**
 * Image
 **/
Image::Image (int width_, int height_){

    assert(width_ > 0);
    assert(height_ > 0);

    width           = width_;
    height          = height_;
    num_pixels      = width * height;
    sampling_method = IMAGE_SAMPLING_POINT;
    
    data.raw = new uint8_t[num_pixels*4];
	int b = 0; //which byte to write to
	for (int j = 0; j < height; j++){
		for (int i = 0; i < width; i++){
			data.raw[b++] = 0;
			data.raw[b++] = 0;
			data.raw[b++] = 0;
			data.raw[b++] = 0;
		}
	}

    assert(data.raw != NULL);
}

Image::Image (const Image& src){
	
	width           = src.width;
    height          = src.height;
    num_pixels      = width * height;
    sampling_method = IMAGE_SAMPLING_POINT;
    
    data.raw = new uint8_t[num_pixels*4];
    
    //memcpy(data.raw, src.data.raw, num_pixels);
    *data.raw = *src.data.raw;
}

Image::Image (char* fname){

	int numComponents; //(e.g., Y, YA, RGB, or RGBA)
	data.raw = stbi_load(fname, &width, &height, &numComponents, 4);
	
	if (data.raw == NULL){
		printf("Error loading image: %s", fname);
		exit(-1);
	}
	

	num_pixels = width * height;
	sampling_method = IMAGE_SAMPLING_POINT;
	
}

Image::~Image (){
    delete data.raw;
    data.raw = NULL;
}

void Image::Write(char* fname){
	
	int lastc = strlen(fname);

	switch (fname[lastc-1]){
	   case 'g': //jpeg (or jpg) or png
	     if (fname[lastc-2] == 'p' || fname[lastc-2] == 'e') //jpeg or jpg
	        stbi_write_jpg(fname, width, height, 4, data.raw, 95);  //95% jpeg quality
	     else //png
	        stbi_write_png(fname, width, height, 4, data.raw, width*4);
	     break;
	   case 'a': //tga (targa)
	     stbi_write_tga(fname, width, height, 4, data.raw);
	     break;
	   case 'p': //bmp
	   default:
	     stbi_write_bmp(fname, width, height, 4, data.raw);
	}
}

void Image::AddNoise(double intensity)
{
	int x, y;
	for (x = 0; x < Width(); x++)
	{
		for (y = 0; y < Height(); y++)
		{
			Pixel p = GetPixel(x, y);
			double r, g, b;
			if ((rand() % 10) < 5) {
				r = (double)p.r + intensity * (double)PixelRandom().r;
				g = (double)p.g + intensity * (double)PixelRandom().g;
				b = (double)p.b + intensity * (double)PixelRandom().b;
			}
			else {
				r = (double)p.r - intensity * (double)PixelRandom().r;
				g = (double)p.g - intensity * (double)PixelRandom().g;
				b = (double)p.b - intensity * (double)PixelRandom().b;
			}
			Pixel noised_p;
			noised_p.SetClamp(r, g, b);
			GetPixel(x, y) = noised_p;
		}
	}
}

void Image::Brighten (double factor)
{
	int x,y;
	for (x = 0 ; x < Width() ; x++)
	{
		for (y = 0 ; y < Height() ; y++)
		{
			Pixel p = GetPixel(x, y);
			Pixel scaled_p = p*factor;
			GetPixel(x,y) = scaled_p;
		}
	}
}

void Image::ChangeContrast (double factor)
{
	int x, y;
	double r, g, b;
	for (x = 0; x < Width(); x++)
	{
		for (y = 0; y < Height(); y++)
		{
			Pixel p = GetPixel(x, y);
			Pixel contrasted_p;
			double luninance = 0.3*(double)p.r + 0.59*(double)p.g + 0.11*(double)p.b;
			r = ((double)p.r - luninance) * factor + luninance;
			g = ((double)p.g - luninance) * factor + luninance;
			b = ((double)p.b - luninance) * factor + luninance;
			contrasted_p.SetClamp(r, g, b);
			GetPixel(x, y) = contrasted_p;
		}
	}
}

void Image::ChangeSaturation(double factor)
{
	int x, y;
	double r, g, b;
	for (x = 0; x < Width(); x++)
	{
		for (y = 0; y < Height(); y++)
		{
			Pixel p = GetPixel(x, y);
			Pixel satured_p;
			double gray = ((double)p.r + (double)p.g + (double)p.b)/3.0;
			r = ((double)p.r - gray) * factor + gray;
			g = ((double)p.g - gray) * factor + gray;
			b = ((double)p.b - gray) * factor + gray;
			satured_p.SetClamp(r, g, b);
			GetPixel(x, y) = satured_p;
		}
	}
}

Image* Image::Crop(int x_, int y_, int w, int h)
{
	int x, y;
	double r, g, b;
	Image *subimage = new Image(w,h);
	subimage->width = w;
	subimage->height = h;
	if (x_ + w <= Width() || y_ + h <= Height()) {
		for (x = x_; x < w + x_; x++)
		{
			for (y = y_; y < h + y_; y++)
			{
				Pixel p = GetPixel(x, y);
				Pixel scaled_p = p;
				subimage->GetPixel(x-x_, y-y_) = scaled_p;
			}
		}
	}
	else {
		printf("Error! The range overflows!");
	}
	return subimage;
}

void Image::ExtractChannel(int channel)
{
	int x, y;
	//double r, g, b;
	for (x = 0; x < Width(); x++)
	{
		for (y = 0; y < Height(); y++)
		{
			Pixel p = GetPixel(x, y);
			Pixel output_p(0,0,0);

			if (channel == 1) {
				output_p.r = p.r;
			}
			else if (channel == 2) {
				output_p.g = p.g;
			}
			else if (channel == 3) {
				output_p.b = p.b;
			}
			else {
				printf("Error! The channel doesn't exist!");
			}
			GetPixel(x, y) = output_p;
		}
	}

}

void Image::Quantize (int nbits)
{
	if (nbits > 8) {
		printf("Error! The number of bits overflows!");
	}
	else {
		/*int x, y, i;
		double r, g, b;
		double dec = pow(2, 8-nbits);
		Pixel p, quantized_p;
		for (x = 0; x < Width(); x++)
		{
			for (y = 0; y < Height(); y++)
			{
				p = GetPixel(x, y);
				r = trunc(((double)p.r + 1) / dec) * dec;
				g = trunc(((double)p.g + 1) / dec) * dec;
				b = trunc(((double)p.b + 1) / dec) * dec;
				quantized_p.SetClamp(r, g, b);
				GetPixel(x, y) = quantized_p;
			}
		}*/
		int x, y, i;
		double r, g, b;
		Pixel p, quantized_p;
		for (x = 0; x < Width(); x++)
		{
			for (y = 0; y < Height(); y++)
			{
				p = GetPixel(x, y);
				quantized_p = PixelQuant(p,nbits);
				GetPixel(x, y) = quantized_p;
			}
		}
	}
}

void Image::RandomDither (int nbits)
{
	if (nbits > 8) {
		printf("Error! The number of bits overflows!");
	}
	else {
		int x, y, i;
		double r, g, b;
		Pixel p, quantized_p;
		for (x = 0; x < Width(); x++)
		{
			for (y = 0; y < Height(); y++)
			{
				p = GetPixel(x, y);
				p = PixelQuant(p, nbits) + PixelRandom();
				r = trunc((double)p.r + 0.5);
				g = trunc((double)p.g + 0.5);
				b = trunc((double)p.b + 0.5);
				quantized_p.SetClamp(r, g, b);
				GetPixel(x, y) = quantized_p;
			}
		}
	}
}

static int Bayer4[4][4] =
{
    {15,  7, 13,  5},
    { 3, 11,  1,  9},
    {12,  4, 14,  6},
    { 0,  8,  2, 10}
};

void Image::OrderedDither(int nbits)
{
	if (nbits > 8) {
		printf("Error! The number of bits overflows!");
	}
	else {
		int x, y, i, j;
		/*Pixel test(130, 0, 0);
		printf("%lf", (double)PixelQuant(test,2).r);
		*/
		double r, g, b, tr, tg, tb;
		double interval = pow(2, 8 - nbits);
		double interval_map = 255 / (pow(2, nbits) - 1);
		Pixel p, q, quantized_p;
		double t = 192;
		printf("%lf", (floor(t / interval)*interval_map));
		printf("%lf", ((floor(t / interval) + 1)*interval_map));
		for (x = 0; x < Width(); x++)
		{
			for (y = 0; y < Height(); y++)
			{
				p = GetPixel(x, y);
				q = PixelQuant(p, nbits);
				i = x % 4; j = y % 4;
				r = (double)p.r - (double)q.r;
				g = (double)p.g - (double)q.g;
				b = (double)p.b - (double)q.b;
				/*if (r > Bayer4[i][j]) {
					tr = (floor((double)q.r / interval) + 1)*interval_map;
				}
				else {
					tr = floor((double)q.r / interval)*interval_map; 
				}
				if (g > Bayer4[i][j]) {
					tg = (floor((double)q.g / interval) + 1)*interval_map;
				}
				else {
					tg = floor((double)q.g / interval)*interval_map; 
				}
				if (b > Bayer4[i][j]) {
					tb = (floor((double)q.b / interval) + 1)*interval_map;
				}
				else {
					tb = floor((double)q.b / interval)*interval_map; 
				}*/
				tr = floor((double)q.r + 256 / pow(2, nbits)*(Bayer4[i][j] - 0.5));
				tg = floor((double)q.g + 256 / pow(2, nbits)*(Bayer4[i][j] - 0.5));
				tb = floor((double)q.b + 256 / pow(2, nbits)*(Bayer4[i][j] - 0.5));
				quantized_p.SetClamp(tr, tg, tb);
				GetPixel(x, y) = quantized_p;
			}
		}
	}
}

/* Error-diffusion parameters */
const double
    ALPHA = 7.0 / 16.0,
    BETA  = 3.0 / 16.0,
    GAMMA = 5.0 / 16.0,
    DELTA = 1.0 / 16.0;

void Image::FloydSteinbergDither(int nbits)
{
	if (nbits > 8) {
		printf("Error! The number of bits overflows!");
	}
	else {
		int x, y;
		double r, g, b;
		Pixel p, q, floy_p;
		/*for (x = 0; x < Width(); x++) {
			SetPixel(x,1,GetPixel(x, 1));
			SetPixel(x, Width(),GetPixel(x, Width()));
		}
		for (y = 1; y < Height() - 1; y++) {
			SetPixel(1,y, GetPixel(1, y));
		}*/
		double e1r, e1g, e1b, e2r, e2g, e2b, e3r, e3g, e3b, e4r, e4g, e4b;
		for (x = 1; x < Width(); x++)
		{
			for (y = 1; y < Height()-1; y++)
			{
				p = GetPixel(x, y);
				q = PixelQuant(p, nbits);

				Pixel p1 = GetPixel(x, y - 1); 
				Pixel p2 = GetPixel(x - 1, y + 1);
				Pixel p3 = GetPixel(x - 1, y);
				Pixel p4 = GetPixel(x - 1, y - 1);

				e1r = (double)p1.r - (double)PixelQuant(p1, nbits).r;
				e1g = (double)p1.g - (double)PixelQuant(p1, nbits).g;
				e1b = (double)p1.b - (double)PixelQuant(p1, nbits).b;

				e2r = (double)p2.r - (double)PixelQuant(p2, nbits).r;
				e2g = (double)p2.g - (double)PixelQuant(p2, nbits).g;
				e2b = (double)p2.b - (double)PixelQuant(p2, nbits).b;

				e3r = (double)p3.r - (double)PixelQuant(p3, nbits).r;
				e3g = (double)p3.g - (double)PixelQuant(p3, nbits).g;
				e3b = (double)p3.b - (double)PixelQuant(p3, nbits).b;

				e4r = (double)p4.r - (double)PixelQuant(p4, nbits).r;
				e4g = (double)p4.g - (double)PixelQuant(p4, nbits).g;
				e4b = (double)p4.b - (double)PixelQuant(p4, nbits).b;

				r = (double)q.r + ALPHA*e1r + BETA*e2r + GAMMA*e3r + DELTA*e4r;
				g = (double)q.g + ALPHA*e1g + BETA*e2g + GAMMA*e3g + DELTA*e4g;
				b = (double)q.b + ALPHA*e1b + BETA*e2b + GAMMA*e3b + DELTA*e4b;
				floy_p.SetClamp(r, g, b);
				GetPixel(x, y) = floy_p;
			}
		}
	}
}

void Image::Blur(int n)
{
	int x, y, i, j;
	double r = 0, g = 0, b = 0;
	int radius = (int)(((double)n-1) / 2);
	Pixel p(0,0,0), blurred_p;

	for (x = 0; x < radius; x++) {
		for (y = 0; y < Height(); y++) {
			SetPixel(x,y,GetPixel(x,y));
		}
	}
	for (x = Width()-radius; x < Width(); x++) {
		for (y = 0; y < Height(); y++) {
			SetPixel(x, y, GetPixel(x, y));
		}
	}
	for (x = radius; x < Width() - radius+1; x++) {
		for (y = 0; y < radius; y++) {
			SetPixel(x, y, GetPixel(x, y));
		}
		for (y = Height()-radius; y < Height(); y++) {
			SetPixel(x, y, GetPixel(x, y));
		}
	}
	for (x = radius; x < Width()-radius; x++) {
		for (y = radius; y < Height()-radius; y++) {
			r = 0, g = 0, b = 0;
			for (i = -radius; i < radius+1; i++) {
				for (j = -radius; j < radius+1; j++) {
					p = GetPixel(x-i, y-j);
					r += (double)p.r * exp(-(i*i + j*j) / 2) / 2 / Pi;
					g += (double)p.g * exp(-(i*i + j*j) / 2) / 2 / Pi;
					b += (double)p.b * exp(-(i*i + j*j) / 2) / 2 / Pi;
				}
			}
			blurred_p.SetClamp(r, g, b);
			GetPixel(x, y) = blurred_p;
		}
	}
}

void Image::Sharpen(int n)
{
	int x, y, i, j;
	double r = 0, g = 0, b = 0;
	int radius = (int)(((double)n - 1) / 2);
	Pixel p(0, 0, 0), blurred_p;

	for (x = 0; x < radius; x++) {
		for (y = 0; y < Height(); y++) {
			SetPixel(x, y, GetPixel(x, y));
		}
	}
	for (x = Width() - radius; x < Width(); x++) {
		for (y = 0; y < Height(); y++) {
			SetPixel(x, y, GetPixel(x, y));
		}
	}
	for (x = radius; x < Width() - radius + 1; x++) {
		for (y = 0; y < radius; y++) {
			SetPixel(x, y, GetPixel(x, y));
		}
		for (y = Height() - radius; y < Height(); y++) {
			SetPixel(x, y, GetPixel(x, y));
		}
	}
	for (x = radius; x < Width() - radius; x++) {
		for (y = radius; y < Height() - radius; y++) {
			r = 0, g = 0, b = 0;
			for (i = -radius; i < radius + 1; i++) {
				for (j = -radius; j < radius + 1; j++) {
					p = GetPixel(x - i, y - j);
					r += (double)p.r * exp(-(i*i + j*j) / 2) / 2 / Pi;
					g += (double)p.g * exp(-(i*i + j*j) / 2) / 2 / Pi;
					b += (double)p.b * exp(-(i*i + j*j) / 2) / 2 / Pi;
				}
			}
			p = GetPixel(x, y);
			r = 1.5*(double)p.r - 0.5*r;
			g = 1.5*(double)p.g - 0.5*g;
			b = 1.5*(double)p.b - 0.5*b;
			blurred_p.SetClamp(r, g, b);
			GetPixel(x, y) = blurred_p;
		}
	}
}

void Image::EdgeDetect()
{
	int x, y, i, j;
	double r = 0, g = 0, b = 0;

	double f[3][3] = {
		{-1,-1,-1},
		{-1,8,-1},
		{-1,-1,-1}
	};

	Image* edge = new Image(Width(), Height());

	for (x = 1; x < Width() - 1; x++) {
		for (y = 1; y < Height() - 1; y++) {
			r = 0, g = 0, b = 0;
			for (i = -1; i < 1 + 1; i++) {
				for (j = -1; j < 1 + 1; j++) {
					Pixel p = GetPixel(x - i, y - j);
					r += (double)p.r * f[i+1][j+1];
					g += (double)p.g * f[i+1][j+1];
					b += (double)p.b * f[i+1][j+1];
				}
			}
			Pixel edge_p;
			edge_p.SetClamp(r, g, b);
			edge->GetPixel(x, y) = edge_p;
		}
	}
	for (x = 1; x < Width() - 1; x++) {
		for (y = 1; y < Height() - 1; y++) {
			Pixel edge_p = edge->GetPixel(x, y);
			GetPixel(x, y) = edge_p;
		}
	}
}

Image* Image::Scale(double sx, double sy)
{
	int x, y, u, v;
	double r, g, b;
	double factorx = (double)Width() / sx;
	double factory = (double)Height() / sy;
	Image *subimage = new Image(sx, sy);
	subimage->width = sx;
	subimage->height = sy;
	for (x = 1; x < sx-1; x++)
	{
		for (y = 1; y < sy-1; y++)
		{
			int i = x*factorx;
			int j = y*factory;
			Pixel p = Sample(i, j);
			Pixel scaled_p = p;
			subimage->GetPixel(x, y) = scaled_p;
		}
	}
	return subimage;
}

Image* Image::Rotate(double angle)
{
	////forward interpolation
	//int x, y, u, v;
	//Image *subimage = new Image(Width(), Height());
	//subimage->width = Width();
	//subimage->height = Height();
	//for (x = 0; x < Width(); x++)
	//{
	//	for (y = 0; y < Height(); y++)
	//	{
	//		u = trunc(x*cos(angle) - y*sin(angle)-0.5*Width()*cos(angle)+0.5*Height()*sin(angle)+0.5*Width());
	//		v = trunc(x*sin(angle) + y*cos(angle)-0.5*Width()*sin(angle)-0.5*Height()*cos(angle)+0.5*Height());
	//		if (u<0 || u>=Width() || v<0 || v>=Height()) {
	//			continue;
	//		}
	//		else {
	//			Pixel p = GetPixel(x, y);
	//			subimage->GetPixel(u, v) = p;
	//		}
	//	}
	//}
	//return subimage;

	//backward interpolation
	int x, y, u, v, i, j;
	double r, g, b;
	int size = trunc(sqrt(Width()*Width()+Height()*Height()))+1;
	Image *subimage = new Image(size, size);
	subimage->width = size;
	subimage->height = size;
	for (x = 0; x < size; x++)
	{
		for (y = 0; y < size; y++)
		{
			i = x*cos(angle) + y*sin(angle) - 0.5*size*cos(angle) - 0.5*size*sin(angle) + 0.5*Width();
			j = -x*sin(angle) + y*cos(angle) + 0.5*size*sin(angle) - 0.5*size*cos(angle) + 0.5*Height();
			if (i<0 || i >= Width() || j<0 || j >= Height()) {
				continue;
			}
			else {
				Pixel p = Sample(i, j);
				subimage->GetPixel(x, y) = p;
			}
		}
	}
	return subimage;
}

void Image::Fun()
{	
	int x, y, u, v;
	double r = Height() / 2, x0 = Width() / 2,y0= Height() / 2;/*fmin(Width(), Height())/2*/;
	double dx, dy, dis, r_, sign;
	double pr, pb, pg;
	for (x= - x0; x < Width() - x0 ; x++)
	{
		for (y =  - y0; y < Height() - y0; y++)
		{
			dx = 1.3*x; dy = 1.3*y;
			dis = dx*dx + dy*dy;
			r_ = sqrt(dx*dx + dy*dy);
			if (r_ <= r) {
				sign = asin(sqrt(r*r - dx*dx - dy*dy) / r);
				if (sign > 0) {
					u = trunc(2 * r*dx*acos(sqrt(r*r - dx*dx - dy*dy) / r) / Pi / sqrt(dx*dx + dy*dy) + x0);
					v = trunc(2 * r*dy*acos(sqrt(r*r - dx*dx - dy*dy) / r) / Pi / sqrt(dx*dx + dy*dy) + y0);
				}
				else {
					u = trunc(2 * r*dx*(-acos(sqrt(r*r - dx*dx - dy*dy) / r)) / Pi / sqrt(dx*dx + dy*dy) + x0);
					v = trunc(2 * r*dy*(-acos(sqrt(r*r - dx*dx - dy*dy) / r)) / Pi / sqrt(dx*dx + dy*dy) + y0);
				}
				
				/*u = trunc(2 * r*dx*acos(sqrt(r*r - dx*dx - dy*dy) / r) / Pi / sqrt(dx*dx + dy*dy) + x0);
				v = trunc(2 * r*dy*acos(sqrt(r*r - dx*dx - dy*dy) / r) / Pi / sqrt(dx*dx + dy*dy) + y0);*/
				if (u >= 1 && u < Width() && v >=1 && v < Height()) {
					//printf("hh");
					GetPixel(x+x0, y+y0) = Sample(u,v);
				}
				else {
					Pixel p = GetPixel(x + x0, y + x0);
					GetPixel(x + x0, y + y0) = p;
				}
			}			
		}
	}
}

/**
 * Image Sample
 **/
void Image::SetSamplingMethod(int method)
{
    assert((method >= 0) && (method < IMAGE_N_SAMPLING_METHODS));
    sampling_method = method;
}

Pixel Image::Sample (double u, double v){
	int i, j;
	double r, g, b;
	i = u - trunc(u);
	j = v - trunc(v);
	if (sampling_method == 0) {
		//point sampling
		//printf("0");
		Pixel p;
		p = GetPixel(trunc(u), trunc(v));
		return p;
	}
	else if (sampling_method == 1) {
		//bilinear sampling
		Pixel p1 = GetPixel(trunc(u), trunc(v));
		Pixel p2 = GetPixel(trunc(u), trunc(v) + 1);
		Pixel p3 = GetPixel(trunc(u) + 1, trunc(v));
		Pixel p4 = GetPixel(trunc(u) + 1, trunc(v) + 1);
		r = (1 - i) * (1 - j) * (double)p1.r + (1 - i) * j * (double)p2.r + i * (1 - j) * (double)p3.r + i * j * (double)p4.r;
		g = (1 - i) * (1 - j) * (double)p1.g + (1 - i) * j * (double)p2.g + i * (1 - j) * (double)p3.g + i * j * (double)p4.g;
		b = (1 - i) * (1 - j) * (double)p1.b + (1 - i) * j * (double)p2.b + i * (1 - j) * (double)p3.b + i * j * (double)p4.b;
		//printf("1");
		Pixel p;
		p.SetClamp(r, g, b);
		return p;
	}
	else if (sampling_method == 2) {
		//Gaussian sampling 
		if (u<1 || u >= Width()-1 || v<1 || v >= Height()-1) {
			return (Pixel(0, 0, 0));
		}
		r = 0; g = 0; b = 0;
		for (i = -1; i < 1 + 1; i++) {
			for (j = -1; j < 1 + 1; j++) {
				//printf("hh");
				Pixel t = GetPixel(trunc(u) - i, trunc(v) - j);
				//printf("hh1");
				r += (double)t.r * exp(-(i*i + j*j) / 2) / 2 / Pi;
				g += (double)t.g * exp(-(i*i + j*j) / 2) / 2 / Pi;
				b += (double)t.b * exp(-(i*i + j*j) / 2) / 2 / Pi;
			}
		}
		//printf("2");
		Pixel p;
		p.SetClamp(r, g, b);
		return p;
	}
	else {
		Pixel p;
		printf("Error! No such method!");
		return p;
	}
}