#include <material.h>
#include <scene.h>
#include <iostream>

vec3 Material::shade(const Intersection* isect, uint8_t depth) const {

	
	decimal bias = 1e-4;
	//	for all lights
	//		test for shadow if needed
	//		call shadeLight if not in shadow
	//		accumulate contribution
	const std::vector<std::unique_ptr<Light>>& lights = isect->scene->lights();
	////Pour l'accumulation de la contribution
	vec3 total_light (0);
	total_light += vec3(0.2);
	////initialisation de la position du shadow ray
	vec3 shadow_ray_origin = isect->position + bias * isect->normal;
	
	for (uint i = 0; i < lights.size(); i++){
		bool inShadow;
		//Direction du shadow ray
		vec3 shadow_ray_direction = glm::normalize(lights.at(i)->positionOrDirection-shadow_ray_origin);
		Ray shadow_ray = Ray{ shadow_ray_origin, shadow_ray_direction };

		if (isect->scene->trace(shadow_ray, 1) == nullptr)
			inShadow = false;
		else
			inShadow = true;
		
		if (!inShadow)
			total_light += this->shadeLight(isect, lights[i].get(), depth);
			
			
	}
	return total_light;
	
}

//NE PAS OUBLIER QUE LES LUMIERES DIRECTIONNELLES SONT NORMALISEES DANS SCENE.cpp

vec3 Material::shadeLight(const Intersection* isect, const Light* l, uint8_t depth) const {
	// facteur d'echelle pour le damier
	float scale = 10.0f;
	float u = isect->uv.x;
	//std::cout << "u is" << u << std::endl;
	float v = isect->uv.y;
	//std::cout << "v is" << v << std::endl;
	vec3 color;
	//std::cout << "floor u is" << floorf(scale * u) << std::endl;
	//std::cout << "floor v is" << floorf(scale * v) << std::endl;
	if ((int)(floorf(scale * u) + floorf(scale * v)) % 2 == 1){
		color = vec3(0.0f);
	}
	else
		color = vec3(1.0f);
	//vec3 color = glm::vec3(1.0f, 0.5f, 0.7f);
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