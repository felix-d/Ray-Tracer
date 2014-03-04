#include <geom.h>
#include <basic_structs.h>
#include <glm/gtx/euler_angles.hpp>
#include <iostream>
//Jai juste tout ecrit ce qui a dans le header, jai fait les constructeurs etc..
Geometry::Geometry(vec3 position, vec3 orientation, vec3 scaling, Material* mtl)
:_position(position),
_orientation(orientation),
_scaling(scaling),
_material(mtl)
{
	
	mat4 translation_mat = glm::translate(mat4(), position);
	mat4 scaling_mat = glm::scale(mat4(), scaling);
	mat4 orientation_mat = glm::orientate4(orientation);
	_modelTransform ={
		
		translation_mat
	};
	
	
	
	
}

Sphere::Sphere(vec3 position, vec3 orientation, vec3 scaling, Material* mtl)
:Geometry(position, orientation, scaling, mtl){
    //TODO implementer constructeur Sphere
    //L'equation de la sphere est donnee par (X-C).(X-C)=r^2
    //C est est le centre et r est le rayon
	_radius = 1.0f;
	_center = position;
	
}

std::unique_ptr<struct Intersection> Sphere::intersect(const struct Ray& ray, decimal &currentdepth) const{
   
	vec3 m = ray.origin - _center;
	decimal b = glm::dot(m, ray.direction);
	decimal c = glm::dot(m, m) - _radius * _radius;
	if (c > 0.0f && b > 0.0f)
		return nullptr;
	decimal discr = b * b - c;
	if (discr < 0.0f)
		return nullptr;
	decimal t = -b - sqrt(discr);
	if (t < 0.0f)
		t = 0.0f;
	vec3 ray_isect = ray.origin + (decimal)t * ray.direction;
	//std::cout << "intersection at " << position.x << " " << position.y << " " << position.z << std::endl;
	vec3 normal = glm::normalize(ray_isect - _center);

	// Calcul des coordonees uv
	decimal x = ray_isect.x;
	decimal y = ray_isect.y;
	decimal z = ray_isect.z;

	/*float v = acosf(z / _radius) / glm::pi<float>();
	float u = acosf(x / (_radius * sinf(glm::pi<float>() * v))) / (2.0f * glm::pi<float>());*/

	decimal u;
	if (y > 0)
		u = acosf(x / sqrtf(powf(x, 2.0f) + powf(y, 2.0f))) / (2.0f * glm::pi<decimal>());
	else
		u = 1 - acosf(x / sqrtf(powf(x, 2.0f) + powf(y, 2.0f))) / (2.0f * glm::pi<decimal>());

	decimal v = acosf(z / _radius) / glm::pi<decimal>();
	
	vec2 uv = glm::vec2(u, v);
	
	std::unique_ptr<struct Intersection> isect(new Intersection{ ray, ray_isect, normal, uv, _material });
	return std::move(isect);
    
}

Box::Box(vec3 position, vec3 orientation, vec3 scaling, Material* mtl)
:Geometry(position, orientation, scaling, mtl){
	//position = vec3(0);
	//La transformation est appliquee dans la fonction de collision
	
	vec4 center4 = vec4(0.0f, 0.0f, 0.0f, 1.0f);
	center4 = _modelTransform * center4;
	_center = vec3(center4);

	_faces.push_back(vec3(1.0f, 0.0f, 0.0f));
	_faces.push_back(vec3(0.0f, 1.0f, 0.0f));
	_faces.push_back(vec3(0.0f, 0.0f, 1.0f));
	_faces.push_back(vec3(-1.0f, 0.0f, 0.0f));
	_faces.push_back(vec3(0.0f, -1.0f, 0.0f));
	_faces.push_back(vec3(0.0f, 0.0f, -1.0f));

	
	for (uint i = 0; i < _faces.size(); i++){
		vec4 current = vec4(_faces[i], 1.0f);
		current =  _modelTransform*current;
		_faces[i] = vec3(current);
	}
	
}

std::unique_ptr<struct Intersection> Box::intersect(const struct Ray& ray, decimal &currentdepth) const{
    //https://code.google.com/p/opengl-tutorial-org/source/browse/misc05_picking/misc05_picking_custom.cpp

	// Intersection method from Real-Time Rendering and Essential Mathematics for Games

	vec3  _min = vec3(-1.0f, -1.0f, -1.0f);
	vec3 _max = vec3(1.0f, 1.0f, 1.0f);
	// Intersection method from Real-Time Rendering and Essential Mathematics for Games

	float tMin = 0.0f;
	float tMax = 100000.0f;

	glm::vec3 OBBposition_worldspace(_modelTransform[3].x, _modelTransform[3].y, _modelTransform[3].z);

	glm::vec3 delta = OBBposition_worldspace - ray.origin;

	// Test intersection with the 2 planes perpendicular to the OBB's X axis
	{
		glm::vec3 xaxis(_modelTransform[0].x, _modelTransform[0].y, _modelTransform[0].z);
		float e = glm::dot(xaxis, delta);
		float f = glm::dot(ray.direction, xaxis);

		if (fabs(f) > 0.001f){ // Standard case

			float t1 = (e + _min.x) / f; // Intersection with the "left" plane
			float t2 = (e + _max.x) / f; // Intersection with the "right" plane
			// t1 and t2 now contain distances betwen ray origin and ray-plane intersections

			// We want t1 to represent the nearest intersection, 
			// so if it's not the case, invert t1 and t2
			if (t1>t2){
				float w = t1; t1 = t2; t2 = w; // swap t1 and t2
			}

			// tMax is the nearest "far" intersection (amongst the X,Y and Z planes pairs)
			if (t2 < tMax)
				tMax = t2;
			// tMin is the farthest "near" intersection (amongst the X,Y and Z planes pairs)
			if (t1 > tMin)
				tMin = t1;

			// And here's the trick :
			// If "far" is closer than "near", then there is NO intersection.
			// See the images in the tutorials for the visual explanation.
			if (tMax < tMin)
				return false;

		}
		else{ // Rare case : the ray is almost parallel to the planes, so they don't have any "intersection"
			if (-e + _min.x > 0.0f || -e + _max.x < 0.0f)
				return false;
		}
	}


	// Test intersection with the 2 planes perpendicular to the OBB's Y axis
	// Exactly the same thing than above.
	{
		glm::vec3 yaxis(_modelTransform[1].x, _modelTransform[1].y, _modelTransform[1].z);
		float e = glm::dot(yaxis, delta);
		float f = glm::dot(ray.direction, yaxis);

		if (fabs(f) > 0.001f){

			float t1 = (e + _min.y) / f;
			float t2 = (e + _max.y) / f;

			if (t1>t2){ float w = t1; t1 = t2; t2 = w; }

			if (t2 < tMax)
				tMax = t2;
			if (t1 > tMin)
				tMin = t1;
			if (tMin > tMax)
				return false;

		}
		else{
			if (-e + _min.y > 0.0f || -e + _max.y < 0.0f)
				return false;
		}
	}


	// Test intersection with the 2 planes perpendicular to the OBB's Z axis
	// Exactly the same thing than above.
	{
		glm::vec3 zaxis(_modelTransform[2].x, _modelTransform[2].y, _modelTransform[2].z);
		float e = glm::dot(zaxis, delta);
		float f = glm::dot(ray.direction, zaxis);

		if (fabs(f) > 0.001f){

			float t1 = (e + _min.z) / f;
			float t2 = (e + _max.z) / f;

			if (t1>t2){ float w = t1; t1 = t2; t2 = w; }

			if (t2 < tMax)
				tMax = t2;
			if (t1 > tMin)
				tMin = t1;
			if (tMin > tMax)
				return false;

		}
		else{
			if (-e + _min.z > 0.0f || -e + _max.z < 0.0f)
				return false;
		}
	}

	vec3 ray_isect = ray.origin + ray.direction* (decimal)tMin;
	
	decimal min_face_dist = INFINITY;
	uint index;
	for (uint i = 0; i < _faces.size(); i++){
		decimal curr_face_dist = length(_faces[i] - ray_isect);
		if (curr_face_dist < min_face_dist){
			min_face_dist = curr_face_dist;
			index = i;
		}
	}
	vec3 normal = normalize(_faces[index] - _center);
	std::unique_ptr<struct Intersection> isect(new Intersection{ ray, ray_isect, normal, vec2(0), _material });
	return std::move(isect);
	}

Cylinder::Cylinder(vec3 position, vec3 orientation, vec3 scaling, Material* mtl)
:Geometry(position, orientation, scaling, mtl){
	//TODO implementer constructeur cylindre
}

std::unique_ptr<struct Intersection> Cylinder::intersect(const struct Ray& ray, decimal &currentdepth) const{
	//TODO implementer intersection entre cylindre et rayon.
	return NULL;
}

Cone::Cone(vec3 position, vec3 orientation, vec3 scaling, Material* mtl)
:Geometry(position, orientation, scaling, mtl){
	//TODO implementer constructeur cone
}

std::unique_ptr<struct Intersection> Cone::intersect(const struct Ray& ray, decimal &currentdepth) const{
	//TODO implementer intersection entre cone et rayon.
	return NULL;
}

