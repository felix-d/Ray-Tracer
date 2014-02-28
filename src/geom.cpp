#include <geom.h>
#include <basic_structs.h>
#include <glm/gtx/euler_angles.hpp>

Geometry::Geometry(vec3 position, vec3 orientation, vec3 scaling, Material* mtl)
:_position(position),
_orientation(orientation),
_scaling(scaling),
_material(mtl)
{}

Sphere::Sphere(vec3 position, vec3 orientation, vec3 scaling, Material* mtl)
:Geometry(position,orientation,scaling,mtl){
    //TODO implementer constructeur Sphere
    //L'equation de la sphere est donnee par (X-C).(X-C)=r^2
    //C est est le centre et r est le rayon
	_radius = 1.0f * (float)scaling.x;
	_center = position;
}

std::unique_ptr<struct Intersection> Sphere::intersect(const struct Ray& ray, decimal &currentdepth) const{
   

	vec3 m = ray.origin - _center;
	double b = glm::dot(m, ray.direction);
	double c = glm::dot(m, m) - _radius*_radius;
	if (c > 0.0f && b > 0.0f) return 0;
	double discr = b*b - c;
	if (discr < 0.0f) return 0;
	double t = -b - sqrt(discr);
	if (t < 0.0f) t = 0.0f;
	vec3 q = ray.origin + t*ray.direction;
	vec3 normal = glm::normalize(q - _center);
	
	//std::unique_ptr<Intersection> intersection = Intersection { ray, q, normal, uv, _material, scene };
	return NULL;

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

