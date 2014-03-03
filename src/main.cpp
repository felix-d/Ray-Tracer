#include <main.h>

#include <fstream>
#include <iostream>
#include <iomanip>

#include <basic_structs.h>
#include <scene.h>
#include <material.h>

bool use_fresnel = false;

inline int discrete(decimal v, decimal max_val)
{
	return int(min(v * max_val, 255.0));
}

int main(int argc, const char* argv[])
{
	/////////////////////////////////
	// Step 0: Parse run arguments //
	/////////////////////////////////

	std::string outfilename = "image.ppm";
	std::string infilename = "../../scenes/phong.scn";

	uint width = 1024;
	uint height = 768;
	uint samples = 1;

	// Simple tokenization scheme
	{
		bool in_defined = false, out_defined = false;
		for (int x = 1; x < argc; x++)
		{
			if (0 == _stricmp(argv[x], "-w") || 0 == _stricmp(argv[x], "--width"))
				width = atoi(argv[++x]);
			else if (0 == _stricmp(argv[x], "-h") || 0 == _stricmp(argv[x], "--height"))
				height = atoi(argv[++x]);
			else if (0 == _stricmp(argv[x], "-s") || 0 == _stricmp(argv[x], "--samples"))
				samples = atoi(argv[++x]);
			else if (0 == _stricmp(argv[x], "-f") || 0 == _stricmp(argv[x], "--fresnel"))
				use_fresnel = true;
			else if (0 == _stricmp(argv[x], "-l") || 0 == _stricmp(argv[x], "--log"))
				Log::SetFile(argv[++x]);
			else
			{
				if (!in_defined)
				{
					infilename = argv[x];
					in_defined = true;
				}
				else if (!out_defined)
				{
					outfilename = argv[x];
					out_defined = true;
				}
				else
					_LOG_WARN() << "Unrecognized token '" << argv[x] << "'" << std::endl;
			}
		}
	}

	std::cout << "Using settings:" << std::endl
		 << "width = " << width << std::endl
		 << "height = " << height << std::endl
		 << "samples = " << samples << std::endl
		 << "fresnel = " << (use_fresnel ? "true" : "false") << std::endl << std::endl
		 << "Loading scene '" << infilename << "'." << std::endl << std::endl;

	/////////////////////////////////////
	// Step 1: Initialize image, scene //
	/////////////////////////////////////

	std::vector<vec3> image(width*height);
	Scene scene(infilename.c_str());
	std::cout << "Scene initialized successfully." << std::endl;
	std::cout.precision(2);
	std::cout << std::fixed;

	////////////////////////////////////
	// Step 2: Initialize render data //
	////////////////////////////////////

	uint8_t max_depth = scene.maxDepth();
	vec3 origin = glm::vec3(0.0f);
	vec3 direction = glm::vec3(1.0f);
	decimal maxdist = 1;
	decimal mindist = 1;
	// compteur indiquant la position a laquelle on est rendu dans le vecteur image
	uint image_pos = 0;

	////////////////////////////
	// Step 3: Perform render //
	////////////////////////////

	for (uint x_pixel = 0; x_pixel < width; x_pixel++) {
		for (uint y_pixel = 0; y_pixel < width; y_pixel++) {
			Ray ray = Ray{ origin, direction };
			std::unique_ptr<Intersection> isect = scene.trace(ray, max_depth, maxdist, mindist);
			if (isect == nullptr)
				image[image_pos] = scene.background();
			//else
				//image[image_pos] = isect->material->shade(const_cast<const Intersection*>(isect), max_depth);
		}
	}

	//	loop over pixels
	//		loop over subpixels
	//			generate ray
	//			scene.trace
	//			if intersection: material->shade
	//			else: scene.background
	//			add color to subpixel

	decimal sample_norm = 1.0 / (samples * samples);

	std::ofstream outfile(outfilename);
	outfile << "P3" << std::endl << width << " " << height << std::endl << 255 << std::endl;
	for (uint i = 0; i < width * height; i++)
		outfile << " " << discrete(image[i].r * sample_norm, scene.discretization()) << " " << discrete(image[i].g * sample_norm, scene.discretization()) << " " << discrete(image[i].b * sample_norm, scene.discretization()) << std::endl;

	std::cout << "Wrote file '" << outfilename << "'." << std::endl;

	return 0;
}