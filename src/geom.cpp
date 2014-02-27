#include <geom.h>
#include <basic_structs.h>
#include <glm/gtx/euler_angles.hpp>


Sphere::Sphere(vec3 position, vec3 orientation, vec3 scaling, Material* mtl = new Material())
:Geometry(position,orientation,scaling,mtl){
    //TODO implementer constructeur Sphere
}

std::unique_ptr<struct Intersection> Sphere::intersect(const struct Ray& ray, decimal &currentdepth) const{
    //TODO implementer intersection entre une sphere et un rayon
	return NULL;
}

Box::Box(vec3 position, vec3 orientation, vec3 scaling, Material* mtl = new Material())
:Geometry(position, orientation, scaling, mtl){
    //TODO implementer constructeur boite
}

std::unique_ptr<struct Intersection> Box::intersect(const struct Ray& ray, decimal &currentdepth) const{
    //TODO implementer intersection entre boite et rayon.
	return NULL;
}

Cylinder::Cylinder(vec3 position, vec3 orientation, vec3 scaling, Material* mtl = new Material())
:Geometry(position, orientation, scaling, mtl){
	//TODO implementer constructeur cylindre
}

std::unique_ptr<struct Intersection> Cylinder::intersect(const struct Ray& ray, decimal &currentdepth) const{
	//TODO implementer intersection entre cylindre et rayon.
	return NULL;
}

Cone::Cone(vec3 position, vec3 orientation, vec3 scaling, Material* mtl = new Material())
:Geometry(position, orientation, scaling, mtl){
	//TODO implementer constructeur cone
}

std::unique_ptr<struct Intersection> Cone::intersect(const struct Ray& ray, decimal &currentdepth) const{
	//TODO implementer intersection entre cone et rayon.
	return NULL;
}

