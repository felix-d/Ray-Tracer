#include <geom.h>
#include <basic_structs.h>
#include <glm/gtx/euler_angles.hpp>
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
	
	/*Salut Félix,

		Tu peux faire beaucoup plus simple.Construis ton objet Intersection
		dans geom.cpp, retourne - le(sans le pointeur vers la scène), puis,
		dans Scene::trace, tu peux juste faire isect->scene = this; avant de
		retourner isect.

		Ça simplifie pas mal les choses!

		Jean - Philippe

		2014 - 02 - 28 10:32 GMT - 05 : 00 felix d <felix.descoteaux@hotmail.com> :
	> J'ai de la difficulte a comprendre comment on passe le pointeur vers la
	> scene au struct Intersection que lon doit retourner dans la definition de la
	> fonction intersect dans geom.cpp.
	> Doit on rajouter un pointeur vers la scene dans le constructeur de geometry ?
	> (Un peu comme le material ? )
	>
	> Merci!*/
	//std::unique_ptr<Intersection> intersection = Intersection { ray, q, normal, vec2(0), _material, NULL };
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

