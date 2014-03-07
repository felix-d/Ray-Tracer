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

std::vector<vec3> superSampling(const uint& samples, const uint& x_pixel, const uint& y_pixel, 
	const Scene& scene, const decimal& width, const decimal& height);

int main(int argc, const char* argv[])
{
	/////////////////////////////////
	// Step 0: Parse run arguments //
	/////////////////////////////////

	std::string outfilename = "image.ppm";
	std::string infilename = "../../scenes/interreflect.scn";

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

	uint8_t max_depth = 5;//scene.maxDepth();
	vec3 origin = vec3(scene.cameraMatrix() * glm::vec4(0.0f, 0, 0, 1));
	uint image_pos = 0;
	
	////////////////////////////
	// Step 3: Perform render //
	////////////////////////////

	//VOIR http://www.scratchapixel.com/lessons/3d-basic-lessons/lesson-6-rays-cameras-and-images/building-primary-rays-and-rendering-an-image/
    
	for (uint y_pixel = 0; y_pixel < height; y_pixel++) {
		for (uint x_pixel = 0; x_pixel < width; x_pixel++) {
			std::vector<vec3> ps = superSampling(samples, x_pixel, y_pixel, scene, width, height);
			std::vector<vec3> rgbs;
			for (uint i = 0; i < ps.size(); i++){
				vec3 p = ps[i];
				vec3 direction = glm::normalize(p - origin);
				Ray ray = Ray{ origin, direction }; 
				//std::cout << (int)max_depth;
				std::unique_ptr<Intersection> isect = scene.trace(ray, 0);
				//std::cout << (int)max_depth;
				if (isect == nullptr) rgbs.push_back(scene.background());
				else {
					rgbs.push_back(isect->material->shade(isect.get(), 0));
					//if((int)currentDepth()<10)std::cout << (int)currentDepth() << std::endl;
				}
				
			}
			image[image_pos] = averageVec3s(rgbs);
			image_pos++;
		}
	}

	decimal sample_norm = 1.0 / (samples * samples);
	std::ofstream outfile(outfilename);
	outfile << "P3" << std::endl << width << " " << height << std::endl << 255 << std::endl;
	for (uint i = 0; i < width * height; i++)
		outfile << " " << discrete(image[i].r * sample_norm, scene.discretization()) << " " << discrete(image[i].g * sample_norm, scene.discretization()) << " " << discrete(image[i].b * sample_norm, scene.discretization()) << std::endl;

	std::cout << "Wrote file '" << outfilename << "'." << std::endl;
	std::string quit;
	std::cin >> quit;

	return 0;
}


std::vector<vec3> superSampling(const uint& samples, const uint& x_pixel, const uint& y_pixel, 
	const Scene& scene, const decimal& width, const decimal& height) {
	decimal angle = tan(scene.fov() / 2);
	decimal invWidth = 1 / (decimal)width, invHeight = 1 / (decimal)height;
	decimal aspectratio = (decimal)width / (decimal)height;
	std::vector<vec3> points;
	std::vector<decimal>xxs, yys;
	float samples2 = pow(samples, 2);
	for (uint k = 0; k < samples; k++){
		decimal indice;
		if (k != 0)
			indice = (k*samples + (k + 1)*samples) / (2 * samples2);
		else indice = samples / (2 * samples2);
		xxs.push_back((2 * ((x_pixel + indice) * invWidth) - 1) * angle * aspectratio);
		yys.push_back((1 - 2 * ((y_pixel + indice) * invHeight)) * angle);
	}
	for (uint k = 0; k < xxs.size(); k++){
		for (uint h = 0; h < yys.size(); h++){
			vec4 p_homog = scene.cameraMatrix() * vec4{ xxs[k], yys[h], -1, 1 };
			points.push_back(vec3(p_homog));
		}
	}
	return points;
}