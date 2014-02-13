#include <material.h>
#include <scene.h>

vec3 Material::shade(const Intersection* isect, uint8_t depth) const
{
	//	for all lights
	//		test for shadow if needed
	//		call shadeLight if not in shadow
	//		accumulate contribution
}