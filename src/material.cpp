#include <material.h>
#include <scene.h>
#include <iostream>

vec3 Material::shade(const Intersection* isect, uint8_t depth) const {
	depth++;
	decimal bias = 1e-4;
	const std::vector<std::unique_ptr<Light>>& lights = isect->scene->lights();
	vec3 total_light(0.0f);
	vec3 shadow_ray_origin = isect->position + bias * isect->normal;
	uint nb_lights = lights.size();

	for (uint i = 0; i < nb_lights; i++){
		bool inShadow = false;
		bool checkForDropShadows = true;
		vec3 shadow_ray_direction;

		//Light type check
		if (lights[i]->type >= lights[i]->NO_SHADOWS)
			inShadow = false;
		else{
			if (lights[i]->directional())
				shadow_ray_direction = glm::normalize(-lights.at(i)->positionOrDirection);
			else
				shadow_ray_direction = glm::normalize(lights.at(i)->positionOrDirection - shadow_ray_origin);

			//Shadow ray
			Ray shadow_ray = Ray{ shadow_ray_origin, shadow_ray_direction };
			if (isect->scene->trace(shadow_ray, 1) == nullptr)
				inShadow = false;
			else
				inShadow = true;
		}

		if (!inShadow)
			total_light += this->shadeLight(isect, lights[i].get(), depth);
		
	}
	return total_light/ (double)nb_lights;
}

//NE PAS OUBLIER QUE LES LUMIERES DIRECTIONNELLES SONT NORMALISEES DANS SCENE.cpp

vec3 Material::shadeLight(const Intersection* isect, const Light* l, uint8_t depth) const {
	decimal coefficient = 0.35;
	
	float scale = 10.0f;
	float u = isect->uv.x;
	float v = isect->uv.y;
	vec3 color;
	if (_texture != nullptr){
		color = _texture->get(_texture->width() * (1 - v), _texture->width()*u);
	}
	else{
		if ((int)(floorf(scale * u) + floorf(scale * v)) % 2 == 1){
			color = vec3(1.0f);
		}
		else
			color = vec3(0.0f);
	}
	if (l->directional()){
		
		//std::cout << lambert << std::endl;
		//if (lambert < 0) std::cout << lambert<<std::endl;
		color = coefficient * color * l->color;
	}
	else {
		double dist = glm::length(l->positionOrDirection - isect->position);
		color = coefficient * (color)* ((l->color*pi()) / pow(dist, 2.0));
	}
    //IL FAUT IMPLEMENTER LATTRIBUTION DE LA LUMIERE SUR LE MATERIAU PAR DEFAUT
	//color *= (l->color)/pi();
	return color;
}


vec3 MaterialLambert::shadeLight(const Intersection* isect, const Light* l, uint8_t depth) const {
	//double lambert = max(glm::dot(l->positionOrDirection, isect->normal), 0.0);
	vec3 color;

	if (l->directional()){
		double lambert = max(glm::dot(-l->positionOrDirection, isect->normal), 0.0);
		color = lambert * (_color) * l->color;
	}
	else {
		double lambert = max(glm::dot(-normalize(isect->position-l->positionOrDirection), isect->normal), 0.0);
		double dist = glm::length(l->positionOrDirection - isect->position);
		color = lambert * (_color) * ((l->color*pi()) / pow(dist, 2.0));
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

//
vec3 MaterialReflective::shade(const Intersection* isect, uint8_t depth) const {
	depth++;
	//On peut ajouter un coefficient si on veut que le materiel reflechisse aussi la lumiere.
	decimal coefficient = 0.11;
	decimal bias = 1e-4;
	const std::vector<std::unique_ptr<Light>>& lights = isect->scene->lights();
	vec3 total_light(0.0);
	vec3 shadow_ray_origin = isect->position + bias * isect->normal;
	for (uint i = 0; i < lights.size(); i++){
		bool inShadow = false;
		bool checkForDropShadows = true;
		vec3 shadow_ray_direction;

		//Light type check
		if (lights[i]->type >= lights[i]->NO_SHADOWS)
			inShadow = false;
		else{
			if (lights[i]->directional())
				shadow_ray_direction = glm::normalize(-lights.at(i)->positionOrDirection);
			else
				shadow_ray_direction = glm::normalize(lights.at(i)->positionOrDirection - shadow_ray_origin);

			//Shadow ray
			Ray shadow_ray = Ray{ shadow_ray_origin, shadow_ray_direction };
			if (isect->scene->trace(shadow_ray, 1) == nullptr)
				inShadow = false;
			else
				inShadow = true;
		}

		//Si ce n'est pas dans l'ombre, sinon cest automatiquement noir
		//il est possible de rajouter une couleur d'ombre au lieu en rajoutant
		//un else
		if (!inShadow){
			//Si on a pas encore atteint la profondeur maximale
			if (depth != isect->scene->maxDepth()){
				vec3 v = isect->ray.direction;
				vec3 n = isect->normal;
				vec3 r = v - n * 2.0 * dot(v, n);
				Ray ray_reflex = Ray{ isect->position + bias*isect->normal, r };

				std::unique_ptr<Intersection> isect2 = isect->scene->trace(ray_reflex, depth);
				if (isect2 != nullptr){
					total_light = _color* vec3(isect2->material->shade(isect2.get(), depth));
				}
				else total_light = _color * isect->scene->background();
			}
			//sinon on attribue a lintersection la couleur de lobjet, tout simplement
			else total_light *= _color;

			if (lights[i]->directional()){

				//std::cout << lambert << std::endl;
				//if (lambert < 0) std::cout << lambert<<std::endl;
				total_light = coefficient * total_light * lights[i]->color;
			}
			else {
				double dist = glm::length(lights[i]->positionOrDirection - isect->position);
				total_light = coefficient * (total_light)* ((lights[i]->color*pi()) / pow(dist, 2.0));
			}
			
		}
	}
	return total_light/(decimal)lights.size();
}


vec3 MaterialRefractive::shade(const Intersection* isect, uint8_t depth) const {
	depth++;
	//early exit
	if (depth >= isect->scene->maxDepth())return _color;
	
	//initialization
	std::unique_ptr<Intersection> refract;
	vec3 total_light = vec3(0.0);
	vec3 refrdir = vec3(0.0);
	decimal offset = 1e-4;
	decimal index_of_refraction = _index;
	decimal eta;
	vec3 n = isect->normal;
	vec3 d = isect->ray.direction;
	bool inside = false;

	//calculs
	if (dot(d, n) > 0){
		n = -n;
		inside = true;
	}
	if (inside)
		eta = index_of_refraction;
	else
		eta = 1 / index_of_refraction;
	refrdir = (eta*(d - n*dot(d, n))) - n*sqrt(1 - (pow(eta, 2)*(1 - pow((dot(d, n)), 2))));
	refrdir = refrdir / length(refrdir);
	Ray refractionRay = Ray{ isect->position - n*offset, refrdir };
	refract = isect->scene->trace(refractionRay, depth);
	
	//attribution de la couleur
	if (refract != nullptr)
	total_light = (refract->material->shade(refract.get(),depth) * _color);
	else 
		total_light =  _color;

	return total_light;
}