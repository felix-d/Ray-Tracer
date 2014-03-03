#include <material.h>
#include <scene.h>
#include <iostream>

vec3 Material::shade(const Intersection* isect, uint8_t depth) const {
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


vec3 Material::shadeLight(const Intersection* isect, const Light* l, uint8_t depth) const {
	// facteur d'echelle pour le damier
	float s = 0.5f;
	float u = isect->uv.x;
	//std::cout << "u is" << u << std::endl;
	float v = isect->uv.y;
	//std::cout << "v is" << v << std::endl;
	vec3 color;

	if (static_cast<int>(floorf(s * u) + floorf(s * v)) % 2 == 1)
		color = glm::vec3(0.0f);
	else color = glm::vec3(1.0f);

	return color;
}


vec3 MaterialLambert::shadeLight(const Intersection* isect, const Light* l, uint8_t depth) const {
	return vec3(1);
}


vec3 MaterialBlinnPhong::shadeLight(const Intersection* isect, const Light* l, uint8_t depth) const {
	return vec3(1);
}


vec3 MaterialCombined::shade(const Intersection* isect, uint8_t depth) const {
	return vec3(1);
}


vec3 MaterialCombined::shadeLight(const Intersection* isect, const Light* l, uint8_t depth) const {
	return vec3(1);
}


vec3 MaterialReflective::shade(const Intersection* isect, uint8_t depth) const {
	return vec3(1);
}


vec3 MaterialRefractive::shade(const Intersection* isect, uint8_t depth) const {
	return vec3(1);
}