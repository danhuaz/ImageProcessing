#include "graphics.h"
#include <SDL.h>
#include <SDL_opengl.h>
#include <vector>
#include <stdlib.h>

#define INF_FLOAT 1.e38

#define STB_IMAGE_IMPLEMENTATION //only place once in one .cpp file
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION //only place once in one .cpp files
#include "stb_image_write.h"

using namespace std;

int main(int argc, char *argv[])
{
	SDL_Init(SDL_INIT_VIDEO);  //Initialize Graphics (for OpenGL)

							   //Ask SDL to get a recent version of OpenGL (3.2 or greater)
	SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
	SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 4);
	SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 4);
	SDL_GL_SetAttribute(SDL_GL_STENCIL_SIZE, 8);

	//Set default parameters for the output_image
	//camera
	float px = 0, py = 0, pz = 0, dx = 0, dy = 0, dz = 1, ux = 0, uy = 1, uz = 0, ha = 45;
	//film_resolution
	float width = 640, height = 480;
	//output_image
	char outFile[1024] = "spot.bmp";
	char outFile_bs[1024] = "basicsampling.bmp";
	//background
	FloatPixel background;
	//ambient_light
	FloatPixel ambient;
	//max_depth
	int max_depth = 5;

	//Store all data from the file
	vector <Material> material;
	Material m;
	material.push_back(m);
	vector <Sphere> sphere;
	Sphere s;
	sphere.push_back(s);
	vector <Directional_light> directionlight;
	vector <Point_light> pointlight;
	vector <Spot_light> spotlight;
	

	FILE *fp;
	long length;
	char line[MAX_CHARACTER]; //Assumes no line is longer than 1024 characters!

	//string fileName = "spot_.scn";
	string fileName = "ambient_sphere.scn";

	// open the file containing the scene description
	fp = fopen(fileName.c_str(), "r");

	// check for errors in opening the file
	if (fp == NULL) {
		printf("Can't open file '%s'\n", fileName.c_str());
		return 0;  //Exit
	}

	// determine the file size (this is optional -- feel free to delete the 4 lines below)
	fseek(fp, 0, SEEK_END); // move position indicator to the end of the file;
	length = ftell(fp);  // return the value of the current position
	printf("File '%s' is %ld bytes long.\n\n", fileName.c_str(), length);
	fseek(fp, 0, SEEK_SET);  // move position indicator to the start of the file

	int index = 0; //store the index of material in current loop
	               //material[0] is has default parameter values

							 //Loop through reading each line
	while (fgets(line, MAX_CHARACTER, fp)) { //Assumes no line is longer than 1024 characters!
		if (line[0] == '#') {
			printf("Skipping comment: %s", line);
			continue;
		}

		char command[100];
		int fieldsRead = sscanf(line, "%s ", command); //Read first word in the line (i.e., the command type)

		if (fieldsRead < 1) { //No command read
							  //Blank line
			continue;
		}

		if (strcmp(command, "camera") == 0) { //If the command is a camera command
			sscanf(line, "camera %f %f %f %f %f %f %f %f %f %f", &px, &py, &pz, &dx, &dy, &dz, &ux, &uy, &uz, &ha);
			printf("Camera position (%f,%f,%f) with viewing direction (%f,%f,%f)\n", px, py, pz, dx, dy, dz);
			printf("Up vector: (%f,%f,%f)\tHeight angle: %f\n", ux, uy, uz, ha);

			//Create a window (offsetx, offsety, width, height, flags)
			//SDL_Window* window = SDL_CreateWindow("homework2", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, width, height, SDL_WINDOW_OPENGL);

		}
		else if (strcmp(command, "film_resolution") == 0) { //If the command is a film_resolution command
			sscanf(line, "film_resolution %f %f", &width, &height);
			printf("Width: %f, Height: %f\n", width, height);			
		}
		else if (strcmp(command, "output_image") == 0) { //If the command is an output_image command
			memset(outFile, 0, MAX_CHARACTER);
			sscanf(line, "output_image %s", outFile);
			printf("Render to file named: %s\n", outFile);
		}
		else if (strcmp(command, "sphere") == 0) { //If the command is a sphere command
			Sphere s;
			sscanf(line, "sphere %f %f %f %f", &s.center.x, &s.center.y, &s.center.z, &s.radius);
			printf("Sphere at position (%f,%f,%f) with radius %f\n", s.center.x, s.center.y, s.center.z, s.radius);
			s.material = index;
			sphere.push_back(s);
		}
		else if (strcmp(command, "background") == 0) { //If the command is a background command
			sscanf(line, "background %f %f %f", &background.r, &background.g, &background.b);
			printf("Background color of (%f,%f,%f)\n", background.r, background.g, background.b);
		}
		else if (strcmp(command, "material") == 0) { //If the command is a material command
			Material m;
			sscanf(line, "material %f %f %f %f %f %f %f %f %f %f %f %f %f %f", &m.ambient.r, &m.ambient.g, &m.ambient.b,\
				&m.diffuse.r, &m.diffuse.g, &m.diffuse.b, &m.specular.r, &m.specular.g, &m.specular.b, \
				&m.ns, &m.transmissive.r, &m.transmissive.g, &m.transmissive.b, &m.ior);
			printf("Ambient color: (%f,%f,%f)\tDiffuse color: (%f,%f,%f)\n", m.ambient.r, m.ambient.g, m.ambient.b, m.diffuse.r, m.diffuse.g, m.diffuse.b);
			printf("Specular color: (%f,%f,%f)\tPhong cosine power: %f\n", m.specular.r, m.specular.g, m.specular.b, m.ns);
			printf("Transmisive color: (%f,%f,%f)\tIndex of refraction: %f\n", m.transmissive.r, m.transmissive.g, m.transmissive.b, m.ior);
			material.push_back(m);
			index++;
		}
		else if (strcmp(command, "directional_light") == 0) { //If the command is a directional_light command
			Directional_light dl;
			sscanf(line, "directional_light %f %f %f %f %f %f", &dl.color.r, &dl.color.g, &dl.color.b, &dl.direction.x, &dl.direction.y, &dl.direction.z);
			printf("Color: (%f,%f,%f)\tDirection: (%f,%f,%f)\n", dl.color.r, dl.color.b, dl.color.g, dl.direction.x, dl.direction.y, dl.direction.z);
			directionlight.push_back(dl);
		}
		else if (strcmp(command, "point_light") == 0) { //If the command is a point_light command
			Point_light pl;
			sscanf(line, "point_light %f %f %f %f %f %f", &pl.color.r, &pl.color.g, &pl.color.b, &pl.position.x, &pl.position.y, &pl.position.z);
			printf("Color: (%f,%f,%f)\tLocation: (%f,%f,%f)\n", pl.color.r, pl.color.b, pl.color.g, pl.position.x, pl.position.y, pl.position.z);
			pointlight.push_back(pl);
		}
		else if (strcmp(command, "spot_light") == 0) { //If the command is a spot_light command
			Spot_light sl;
			sscanf(line, "spot_light %f %f %f %f %f %f %f %f %f %f %f", &sl.color.r, &sl.color.g, &sl.color.b, \
				&sl.location.x, &sl.location.y, &sl.location.z, &sl.direction.x, &sl.direction.y, &sl.direction.z, &sl.angle1, &sl.angle2);
			printf("Color: (%f,%f,%f)\tPosition: (%f,%f,%f)\n", sl.color.r, sl.color.g, sl.color.b, sl.location.x, sl.location.y, sl.location.z);
			printf("Direction: (%f,%f,%f)\tangle1: %f\tangle2 %f\n", sl.direction.x, sl.direction.y, sl.direction.z, sl.angle1, sl.angle2);
			spotlight.push_back(sl);
		}
		else if (strcmp(command, "ambient_light") == 0) { //If the command is an ambient_light command
			sscanf(line, "ambient_light %f %f %f", &ambient.r, &ambient.g, &ambient.b);
			printf("Global ambient light: (%f,%f,%f)\n", ambient.r, ambient.g, ambient.b);
		}
		else if (strcmp(command, "max_depth") == 0) { //If the command is a max_depth command
			sscanf(line, "max_depth %d", &max_depth);
			printf("Maximum recursion depth: %d\n", max_depth);
		}
		else {
			printf("WARNING. Do not know command: %s\n", command);
		}
	}

	//the parameter for basic sampling
	//int large = 2;
	/*width = large * width;
	height = large * height;*/

	Image *img = NULL;
	img = new Image(width, height);
	background = 255 * background;
	for (int i = 0; i < img->Width(); i++) {
		for (int j = 0; j < img->Height(); j++) {	
			img->SetPixel(i, j, background.ToPixel());
		}
	}
	img->Write(outFile);

	FloatPixel pixel;

	//create camera coordinate
	float distance = height / 2 / tan(ha*Pi / 180);

	//viewing direction -w
	Point view_direction(dx, dy, dz);
	view_direction.unit();

	//upward vector
	Point UP_vector(ux, uy, uz);
	UP_vector.unit();

	//rightward vector u
	Point rightward_vector = crossproduct(view_direction, UP_vector);
	rightward_vector.unit();

	//upward from the camera v
	Point upward_vector = crossproduct(view_direction.minus(),rightward_vector);
	upward_vector.unit();

	//ray origin
	Point ray_ori(px,py,pz);

	// the pixel at position(i,j) in the image has position(u,v)
	float u, v;
	// process each pixel
	float temp, r;

	Point hitpoint, normal, lightdirection, halfvector, view;
	for (int i = 0; i < img->Width(); i++) {
		for (int j = 0; j < img->Height(); j++) {
			float t = INF_FLOAT;
			bool flag = false;
			int front = 0, mat_idx = 0;

			//compute viewing ray
			u = width / 2 - i + 0.5;
			v = height / 2 - j + 0.5;

			//ray direction
			Point ray_direction = distance*view_direction+u*rightward_vector+v*upward_vector;
			ray_direction.unit();

			//for each sphere calculate t
			for (int sph_idx = 1; sph_idx < sphere.size(); sph_idx++) {
				if (sphere[sph_idx].IfhitSphere(ray_ori, ray_direction)) {//shading on the pixel
					temp = sphere[sph_idx].Hitpoint(ray_ori, ray_direction);
					if (temp >= 0) {
						if (t >= temp) {
							t = temp;
							front = sph_idx;
							flag = true;
						}
					}
				}
			}

			//if hit an sphere
			view = ray_direction.minus();
			view.unit();

			if (flag) {
				//printf("t:%f\n", t);
				hitpoint = ray_ori + t*ray_direction;
				normal = hitpoint - sphere[front].center;
				normal.unit();
				mat_idx = sphere[front].material;

				//add ambient light
				pixel = material[mat_idx].ambient*ambient;
				
				//add point light
				for (int lidx = 0; lidx < pointlight.size(); lidx++) {
					//lightdirection = hitpoint - pointlight[lidx].position;
					lightdirection = pointlight[lidx].position - hitpoint;
					lightdirection.unit();
					halfvector = view + lightdirection;
					halfvector.unit();
					r = pointlight[lidx].position.Distance(hitpoint);
					pixel = pixel + (fmaxf(0, dot(normal, lightdirection))*material[mat_idx].diffuse + \
						powf(fmaxf(0, dot(normal, halfvector)), material[mat_idx].ns)*material[mat_idx].specular)* \
						pointlight[lidx].color *(1 / r)*(1 / r);
				}

				//add directional light
				for (int lidx = 0; lidx < directionlight.size(); lidx++) {
					lightdirection = directionlight[lidx].direction.minus();
					lightdirection.unit();
					float cos_theta = dot(normal, lightdirection);
					if (cos_theta >= 0) {
						halfvector = view + lightdirection;
						halfvector.unit();
						pixel = pixel + (fmaxf(0, cos_theta)*material[mat_idx].diffuse + \
							powf(fmaxf(0, dot(normal, halfvector)), material[mat_idx].ns)*material[mat_idx].specular)*directionlight[lidx].color;
					}
				}

				//add spot light
				for (int lidx = 0; lidx < spotlight.size(); lidx++) {
					spotlight[lidx].direction.unit();
					lightdirection = hitpoint - spotlight[lidx].location;
					lightdirection.unit();

					float cos_angle = dot(spotlight[lidx].direction, lightdirection);
					float angle = acos(cos_angle);
					float angle1 = spotlight[lidx].angle1 * Pi / 180;
					float angle2 = spotlight[lidx].angle2 * Pi / 180;
					lightdirection = lightdirection.minus();
					halfvector = view + lightdirection;
					halfvector.unit();
					r = spotlight[lidx].location.Distance(hitpoint);
					if (cos_angle >= cos(angle1)) {
						pixel = pixel + (fmaxf(0, dot(normal, lightdirection))*material[mat_idx].diffuse + \
							powf(fmaxf(0, dot(normal, halfvector)), material[mat_idx].ns)*material[mat_idx].specular)* \
							spotlight[lidx].color *(1 / r)*(1 / r);
					}
					else if (cos_angle >= cos(angle2)) {
						pixel = pixel + (fmaxf(0, dot(normal, lightdirection))*material[mat_idx].diffuse + \
							powf(fmaxf(0, dot(normal, halfvector)), material[mat_idx].ns)*material[mat_idx].specular)* \
							spotlight[lidx].color *((angle - angle2) / (angle1 - angle2))*(1 / r)*(1 / r);
					}
				}				
				pixel = 255 * pixel;
				img->SetPixel(i, j, pixel.ToPixel());
			}
		}
	}
	
	////basic sampling
	//Image *img_new = NULL;
	//width = width / large;
	//height = height / large;
	//img_new = new Image(width, height);
	//Pixel p;
	//int x, y, i, j;
	//int count = 0;
	//double k = 1.0 / large / large;
	//for (x = 0; x < img_new->Width(); x++) {
	//	for (y = 0; y < img_new->Height(); y++) {
	//		float r = 0, g = 0, b = 0;
	//		for (i = large*x; i < large*(x + 1); i++) {
	//			for (j = large*y; j < large*(y + 1); j++) {
	//				p = img->GetPixel(i, j);
	//				r += (double)p.r * k;
	//				g += (double)p.g * k;
	//				b += (double)p.b * k;
	//				count++;
	//			}
	//		}
	//		p.SetClamp(r, g, b);
	//		img_new->GetPixel(x, y) = p;
	//	}
	//}
	img->Write(outFile);
	//img_new->Write(outFile_bs);
	delete img;
	//delete img_new;

	//SDL_Delay(1000);
	SDL_Quit();
	return 0;
}