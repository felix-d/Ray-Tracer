#include <material.h>
#include <scene.h>

vec3 Material::shade(const Intersection* isect, uint8_t depth) const
{
	//	for all lights
	//		test for shadow if needed
	//		call shadeLight if not in shadow
	//		accumulate contribution
	const std::vector<std::unique_ptr<Light>>& lights = (*(*isect).scene).lights();
	//Pour l'accumulation de la contribution
	vec3 total_light (0);
	//initialisation de la position du shadow ray
	vec3 position = (*isect).position;
	
	for (uint i = 0; i < lights.size(); i++){
		
		bool inShadow;

		//Direction du shadow ray
		vec3 direction = glm::normalize((*lights.at(i)).positionOrDirection - position);
		Ray ray = Ray{ position, direction };

		if ((*(*isect).scene).trace(ray, 1) == nullptr) inShadow = false;
		else inShadow = true;
		
		if (!inShadow)total_light+=this->shadeLight(isect, lights[i].get(), depth);
	}
	return total_light;
}


vec3 Material::shadeLight(const Intersection* isect, const Light* l, uint8_t depth) const{
	return vec3(1);
}


vec3 MaterialLambert::shadeLight(const Intersection* isect, const Light* l, uint8_t depth) const{
	return vec3(1);
}


vec3 MaterialBlinnPhong::shadeLight(const Intersection* isect, const Light* l, uint8_t depth) const{
	return vec3(1);
}


vec3 MaterialCombined::shade(const Intersection* isect, uint8_t depth) const{
	return vec3(1);
}


vec3 MaterialCombined::shadeLight(const Intersection* isect, const Light* l, uint8_t depth) const{
	return vec3(1);
}


vec3 MaterialReflective::shade(const Intersection* isect, uint8_t depth) const{
	return vec3(1);
}


vec3 MaterialRefractive::shade(const Intersection* isect, uint8_t depth) const{
	return vec3(1);
}