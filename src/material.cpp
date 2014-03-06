#include <material.h>
#include <scene.h>
#include <iostream>

vec3 Material::shade(const Intersection* isect, uint8_t depth) const {

	decimal bias = 1e-4;
	const std::vector<std::unique_ptr<Light>>& lights = isect->scene->lights();
	vec3 total_light (0.0f);
	vec3 shadow_ray_origin = isect->position + bias * isect->normal;
	uint nb_lights = lights.size();

	for (uint i = 0; i < nb_lights; i++){
		bool inShadow;
		vec3 shadow_ray_direction = glm::normalize(lights.at(i)->positionOrDirection-shadow_ray_origin);
		Ray shadow_ray = Ray{ shadow_ray_origin, shadow_ray_direction };
		if (isect->scene->trace(shadow_ray, 1) == nullptr)
			inShadow = false;
		else
			inShadow = true;
		//if (!inShadow)
			total_light += this->shadeLight(isect, lights[i].get(), 1);
	}
	return total_light;
}

//NE PAS OUBLIER QUE LES LUMIERES DIRECTIONNELLES SONT NORMALISEES DANS SCENE.cpp

vec3 Material::shadeLight(const Intersection* isect, const Light* l, uint8_t depth) const {
	float scale = 10.0f;
	float u = isect->uv.x;
	float v = isect->uv.y;
	vec3 color;
	if (_texture != nullptr){
		color = _texture->get(_texture->width() * (1-v), _texture->width()*u);
	}
	else{
		if ((int)(floorf(scale * u) + floorf(scale * v)) % 2 == 1){
			color = vec3(1.0f);
		}
		else
			color = vec3(0.0f);
	}
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

//
vec3 MaterialReflective::shade(const Intersection* isect, uint8_t depth) const {

	decrementCurrentDepth();
	if (depth == 0){
		return _color;
	}
	decimal bias = 1e-4;
	const std::vector<std::unique_ptr<Light>>& lights = isect->scene->lights();
	vec3 total_light(0.0f);
	vec3 shadow_ray_origin = isect->position + bias * isect->normal;

	for (uint i = 0; i < lights.size(); i++){
		
		
			bool inShadow;
			vec3 shadow_ray_direction = glm::normalize(lights.at(i)->positionOrDirection - shadow_ray_origin);
			Ray shadow_ray = Ray{ shadow_ray_origin, shadow_ray_direction };
			if (isect->scene->trace(shadow_ray, depth) == nullptr)
				inShadow = false;
			else
				inShadow = true;
			depth++;
			//if (!inShadow){
			vec3 v = isect->ray.direction;
			vec3 n = isect->normal;
			vec3 r = v - n * 2.0 * dot(v,n);
			Ray ray_reflex = Ray{ isect->position + bias*isect->normal, r };
			
			std::unique_ptr<Intersection> isect2 = isect->scene->trace(ray_reflex, currentDepth());
			if (isect2 != nullptr)
				total_light += vec3(isect2->material->shade(isect2.get(), currentDepth()));
			
	}
	return total_light;
}


vec3 MaterialRefractive::shade(const Intersection* isect, uint8_t depth) const {
	return vec3(1);
}