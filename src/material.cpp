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
	vec3 total_light (0.0f);
	////initialisation de la position du shadow ray
	vec3 shadow_ray_origin = isect->position + bias * isect->normal;

	uint nb_lights = lights.size();

	for (uint i = 0; i < nb_lights; i++){
		bool inShadow;
		//Direction du shadow ray
		vec3 shadow_ray_direction = glm::normalize(lights.at(i)->positionOrDirection-shadow_ray_origin);
		Ray shadow_ray = Ray{ shadow_ray_origin, shadow_ray_direction };

		if (isect->scene->trace(shadow_ray, 1) == nullptr || lights.at(i)->type >= lights.at(i)->NO_SHADOWS)
			inShadow = false;
		else
			inShadow = true;
		
		if (!inShadow)
			total_light += this->shadeLight(isect, lights[i].get(), depth);		
	}
	return total_light / (double)nb_lights;
	
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
	double lambert = max(glm::dot(l->positionOrDirection, isect->normal), 0.0);
	vec3 color;

	if (l->directional())
		color = lambert * (_color / glm::pi<double>()) * l->color;
	else {
		double dist = glm::length(l->positionOrDirection - isect->position);
		color = lambert * (_color / glm::pi<double>()) * (l->color / pow(dist, 2));
	}

	return color;
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
	decimal bias = 1e-4;
	//	for all lights
	//		test for shadow if needed
	//		call shadeLight if not in shadow
	//		accumulate contribution
	const std::vector<std::unique_ptr<Light>>& lights = isect->scene->lights();
	////Pour l'accumulation de la contribution
	vec3 total_light(0.0f);
	////initialisation de la position du shadow ray
	vec3 shadow_ray_origin = isect->position + bias * isect->normal;

	for (uint i = 0; i < lights.size(); i++){
		bool inShadow;
		//Direction du shadow ray
		vec3 shadow_ray_direction = glm::normalize(lights.at(i)->positionOrDirection - shadow_ray_origin);
		Ray shadow_ray = Ray{ shadow_ray_origin, shadow_ray_direction };

		if (isect->scene->trace(shadow_ray, 1) == nullptr)
			inShadow = false;
		else
			inShadow = true;

		if (!inShadow)
			total_light += this->shadeLight(isect, lights[i].get(), depth);
			/*vec3 vi_minus_1 = isect->ray.origin - isect->position;
			vec3 ni = isect->normal;
			vec3 ri = normalize(2 * dot(vi_minus_1, ni)*ni - vi_minus_1);
			Ray ray_reflex = Ray{ isect->position, ri };
			std::unique_ptr<Intersection> isect2 = isect->scene->trace(ray_reflex, depth);
			total_light+= vec3(isect2->material->shade(isect2.get(), depth));*/

		//else total_light += vec3(0.2);

	}
	return total_light;
	/*vec3 vi_minus_1 =  isect->ray.origin - isect->position;
	vec3 ni = isect->normal;
	vec3 ri = normalize(2 * dot(vi_minus_1, ni)*ni - vi_minus_1);
	Ray ray_reflex = Ray{ isect->position, ri };
	std::unique_ptr<Intersection> isect2 = isect->scene->trace(ray_reflex, depth);
	return vec3(isect2->material->shade(isect2.get(),depth));*/
}


vec3 MaterialRefractive::shade(const Intersection* isect, uint8_t depth) const {
	return vec3(1);
}