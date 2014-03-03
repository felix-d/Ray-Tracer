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
{}

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
	float b = glm::dot(m, ray.direction);
	float c = glm::dot(m, m) - _radius * _radius;
	if (c > 0.0f && b > 0.0f)
		return nullptr;
	float discr = b * b - c;
	if (discr < 0.0f)
		return nullptr;
	float t = -b - sqrt(discr);
	if (t < 0.0f)
		t = 0.0f;
	vec3 position = ray.origin + (double)t * ray.direction;
	//std::cout << "intersection at " << position.x << " " << position.y << " " << position.z << std::endl;
	vec3 normal = glm::normalize(position - _center);

	// Calcul des coordonees uv
	float x = position.x;
	float y = position.y;
	float z = position.z;

	//float v = acosf(z / _radius) / glm::pi<float>();
	//float u = acosf(x / (_radius * sinf(glm::pi<float>() * v))) / (2.0f * glm::pi<float>());

	float u;
	if (y > 0)
		u = acosf(x / sqrtf(powf(x, 2.0f) + powf(y, 2.0f))) / (2.0f * glm::pi<float>());
	else
		u = 1 - acosf(x / sqrtf(powf(x, 2.0f) + powf(y, 2.0f))) / (2.0f * glm::pi<float>());

	float v = acosf(z / _radius) / glm::pi<float>();
	vec2 uv = glm::vec2(u, v);
	
	std::unique_ptr<struct Intersection> isect(new Intersection{ ray, position, normal, uv, _material });
	return std::move(isect);

}

Box::Box(vec3 position, vec3 orientation, vec3 scaling, Material* mtl)
:Geometry(position, orientation, scaling, mtl){
    //TODO implementer constructeur boite
}

std::unique_ptr<struct Intersection> Box::intersect(const struct Ray& ray, decimal &currentdepth) const{
    //TODO implementer intersection entre boite et rayon.
	return NULL;
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

