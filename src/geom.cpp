﻿#include <geom.h>
#include <basic_structs.h>
#include <glm/gtx/euler_angles.hpp>
#include <iostream>
#include <array>




//CEST NORMAL SI LES TRANSFORMATIONS NE MARCHENT PLUS
//JE LES AI ENLEVER POUR QUON IMPLEMENTE AVEC LA BONNE TECHNIQUE

// 1.Utiliser linverse de la matrice pour calculer le nouveau ray.origin et ray.direction et ainsi trouver lintersection
//   (a la place de modifier l'objet)
//
// 2.A lendroit du hit, etant donner que le rayon transforme frappe l'objet canonique, 
//   il faut retransformer le point touche avec la matrice de transformation pour connaitre le point de l'objet transforme.
//
// 3.Utiliser l'inverse de la transposee pour calculer la normale a partir de la normale de l'objet canonique
//
// 4.Pour les uv, calculer tout simplement les uv de l'objet canonique, ils seront mappes automatiquement.




//https://code.google.com/p/pwsraytracer/source/browse/trunk/raytracer/
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
		translation_mat*
		orientation_mat*
		scaling_mat
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
	vec3 ray_isect = ray.origin + t * ray.direction;
	vec3 normal = glm::normalize(ray_isect - _center);
	vec2 uv = calculateUVSphere(normal);
	std::unique_ptr<struct Intersection> isect(new Intersection{ ray, ray_isect, normal, uv, _material });
	return std::move(isect);
    
}

Box::Box(vec3 position, vec3 orientation, vec3 scaling, Material* mtl)
:Geometry(position, orientation, scaling, mtl){
	_center = vec3(0, 0, 0);
	_min = vec3(-1.0, -1, -1);
	_max = vec3(1.0, 1, 1);

	vec3 corner_ftr{ 1, 1, 1 };//0
	vec3 corner_rbl{-1, -1, -1 };//1
	vec3 corner_fbr{ corner_ftr.x, corner_rbl.y, corner_ftr.z };//2
	vec3 corner_ftl{ corner_ftr.x, corner_ftr.y, corner_rbl.z };//3
	vec3 corner_fbl{ corner_ftr.x, corner_rbl.y, corner_rbl.z };//4
	vec3 corner_rtr{ corner_rbl.x, corner_ftr.y, corner_ftr.z };//5
	vec3 corner_rtl{ corner_rbl.x, corner_ftr.y, corner_rbl.z };//6
	vec3 corner_rbr{ corner_rbl.x, corner_rbl.y, corner_ftr.z };//7
	points = { { corner_ftr, corner_rbl, corner_fbr, corner_ftl,
		corner_fbl, corner_rtr, corner_rtl, corner_rbr } };
	_faces_points.push_back(vec3(1.0f, 0.0f, 0.0f));
	_faces_points.push_back(vec3(0.0f, 1.0f, 0.0f));
	_faces_points.push_back(vec3(0.0f, 0.0f, 1.0f));
	_faces_points.push_back(vec3(-1.0f, 0.0f, 0.0f));
	_faces_points.push_back(vec3(0.0f, -1.0f, 0.0f));
	_faces_points.push_back(vec3(0.0f, 0.0f, -1.0f));
	
}




std::unique_ptr<struct Intersection> Box::intersect(const struct Ray& ray, decimal &currentdepth) const{

	
	float t1 = (_min.x - ray.origin.x)*(1.0 / ray.direction.x);
	float t2 = (_max.x - ray.origin.x)*(1.0 / ray.direction.x);
	float t3 = (_min.y - ray.origin.y)*(1.0 / ray.direction.y);
	float t4 = (_max.y - ray.origin.y)*(1.0 / ray.direction.y);
	float t5 = (_min.z - ray.origin.z)*(1.0 / ray.direction.z);
	float t6 = (_max.z - ray.origin.z)*(1.0 / ray.direction.z);
	decimal t = max(max(min(t1, t2), min(t3, t4)), min(t5, t6));
	float tmax = min(min(max(t1, t2), max(t3, t4)), max(t5, t6));
	if (tmax < 0) return nullptr;
	if (t > tmax) return nullptr;

	vec3 ray_isect = ray.origin + (ray.direction * t);
	
	//Calcul des normales
	decimal min_face_dist = INFINITY;
	uint index;
	for (uint i = 0; i < _faces_points.size(); i++){
		decimal curr_face_dist = length(_faces_points[i] - ray_isect);
		if (curr_face_dist < min_face_dist){
			min_face_dist = curr_face_dist;
			index = i;
		}
	}
	vec3 normal = normalize(_faces_points[index] - _center);

	// Calcul des coordonees uv		
	vec3 uv_coord_0_0;
	vec3 uv_coord_1_0;
	vec3 uv_coord_0_1;

	// test pour savoir sur quel plan se trouve l'intersection
	if (index == 0) {
		uv_coord_0_0 = points[4];
		uv_coord_1_0 = points[2];
		uv_coord_0_1 = points[3];
	}

	if (index == 1) {
		uv_coord_0_0 = points[3];
		uv_coord_1_0 = points[0];
		uv_coord_0_1 = points[6];
	}

	if (index == 2) {
		uv_coord_0_0 = points[2];
		uv_coord_1_0 = points[7];
		uv_coord_0_1 = points[0];
	}

	if (index == 3) {
		uv_coord_0_0 = points[1];
		uv_coord_1_0 = points[7];
		uv_coord_0_1 = points[6];
	}

	if (index == 4) {
		uv_coord_0_0 = points[1];
		uv_coord_1_0 = points[7];
		uv_coord_0_1 = points[4];
	}

	if (index == 5) {
		uv_coord_0_0 = points[1];
		uv_coord_1_0 = points[4];
		uv_coord_0_1 = points[6];
	}

	decimal u = glm::length(glm::cross((ray_isect - uv_coord_0_0), (ray_isect - uv_coord_1_0))) / glm::length(uv_coord_1_0 - uv_coord_0_0);
	decimal v = glm::length(glm::cross((ray_isect - uv_coord_0_0), (ray_isect - uv_coord_0_1))) / glm::length(uv_coord_0_1 - uv_coord_0_0);

	vec2 uv = glm::vec2(u, v);

	std::unique_ptr<struct Intersection> isect(new Intersection{ ray, ray_isect, normal, uv, _material });
	return std::move(isect);
}

Cylinder::Cylinder(vec3 position, vec3 orientation, vec3 scaling, Material* mtl)
:Geometry(position, orientation, scaling, mtl){
	_center =vec3(0.0,0,0);
	_p = vec3(0.0, 1, 0);
	_q = vec3(0.0, -1, 0);
	_radius = 1;
	_height = 2;
}

std::unique_ptr<struct Intersection> Cylinder::intersect(const struct Ray& ray, decimal &currentdepth) const{

https://code.google.com/p/pwsraytracer/source/browse/trunk/raytracer/cylinder.cpp?r=160

	double tmin = INFINITY;
	double t = DBL_MAX, t1 = DBL_MAX, t2 = DBL_MAX, t3 = DBL_MAX;
	bool sides = false;
	double a = pow(ray.direction.x, 2) + pow(ray.direction.z,2);
	double b = 2.0 * ((ray.origin.x - _center.x) * ray.direction.x + (ray.origin.z - _center.z) * ray.direction.z);
	double c = pow(ray.origin.x - _center.x, 2) + pow(ray.origin.z - _center.z,2) - pow(_radius,2);
	double D = b * b - 4.0 * a * c;
	if (D < 0.0){
		return nullptr;
	}
	else{
		double sqrtD = sqrt(D);
		double aa = a * 2;
		double t1 = (-b - sqrtD) / aa;
		if (t1 > epsilon<double>()){
			double y = (ray.origin.y - _center.y) + t1 *  ray.direction.y;
			if (y > _q.y && y < _p.y){
				if (t1 < t) {
					t = t1;
					sides = true;
				}
			}
		}
		t2 = (-b + sqrtD) / aa;
		if (t2 >  epsilon<double>()){
			double y = (ray.origin.y - _center.y) + t2 * ray.direction.y;
			if (y > _q.y && y < _p.y){
				if (t2 < t){
					t = t2;
					sides = true;
				}
			}
		}
	}
	if (plane_Intersection(ray, vec3(0, -1.0, 0), _p, t3)) { 
		vec3 intersection = ray.origin + (ray.direction*t3);
		if (pow(intersection.x - _p.x, 2) + pow(intersection.z - _p.z, 2) <= _radius*_radius) { 
			if (t3 < t){
				t = t3;
				sides = false;
			}
		}
	}
	if (t != DBL_MAX && t>0.00000001){
		vec3 ray_isect = ray.origin + ray.direction*(decimal)t;
		vec3 normal = calculateNormal(ray_isect, sides);
		vec2 uv;
		if (sides)uv = calculateUVCylinder(ray_isect-_center);
		else  uv = calculateUVCircle(ray_isect-_p); 
		std::unique_ptr<struct Intersection> isect(new Intersection{ ray, ray_isect, normal, uv, _material });
		return std::move(isect);
	}
	return nullptr;
}

vec3 Cylinder::calculateNormal(vec3& hitPoint, bool sides)const{
	vec3 normal;
	if (!sides) normal = normalize(_p-_center);
	else {
		vec3 direction = normalize(_p - _q);
		vec3 x = _q + (dot((hitPoint - _q), direction))*direction;
		normal = normalize(hitPoint - x);
	}
	return normal;
}


Cone::Cone(vec3 position, vec3 orientation, vec3 scaling, Material* mtl)
:Geometry(position, orientation, scaling, mtl){
	_radius = 1.0;
	_apex = vec3(0.0, 0.00, 0);
	_base_center = vec3(0.0, -1, 0); 
	_direction = vec3(0.0,1,0);
	_height = 1.00;
}

std::unique_ptr<struct Intersection> Cone::intersect(const struct Ray& ray, decimal &currentdepth) const{
    //https://github.com/Penetra/CG-Project/blob/master/Cone.cpp
	
	double rh = -(_radius*_radius) / (_height*_height);
	double a = pow(ray.direction.x, 2) + pow(ray.direction.z, 2) +
		rh*pow(ray.direction.y, 2);
	double b = 2 * (ray.direction.x*(ray.origin.x - _base_center.x) +
		ray.direction.z*(ray.origin.z - _base_center.z) +
		rh*ray.direction.y*(ray.origin.y - _base_center.y - _height));
	double c = pow(ray.origin.x - _base_center.x, 2) +
		pow(ray.origin.z - _base_center.z, 2) +
		rh*pow(ray.origin.y - _base_center.y - _height, 2);
	double root = b*b - 4.0*a*c;
	if (root<0) return nullptr;
	double t = 0, t1 = DBL_MAX, t2 = DBL_MAX, t3 = DBL_MAX;
	t1 = (-b + sqrtf(root)) / (2.0*a);
	t2 = (-b - sqrtf(root)) / (2.0*a);
	vec3 intersection = (ray.origin) + ((ray.direction) * t1);
	vec3 intersection2 = ray.origin + ((ray.direction) *t2);

	if (intersection.y <= _apex.y && intersection.y >= _base_center.y) {
		if (intersection2.y <= _apex.y && intersection2.y >= _base_center.y) {
			if (t1<t2 && t1>0 && t2 >0)t = t1;
			else if (t2 > 0.0f) t = t2;
		}
		else if (t1 > 0.0f) t = t1;
	}
	else if (intersection2.y <= _apex.y &&
		intersection2.y >= _base_center.y)
	if (t2 > 0.0f) t = t2;


	if (plane_Intersection(ray, vec3(0, -1, 0), _base_center, t3)) {
		intersection = ray.origin + (ray.direction*t3);
		if (pow(intersection.x - _base_center.x, 2) + pow(intersection.z - _base_center.z, 2) <= _radius*_radius)
		if (t3<t) t = t3;
	}

	if (t != DBL_MAX && t>0.00000001){
		vec3 ray_isect = ray.origin + ray.direction*(decimal)t;
		vec3 normal = normalize(ray_isect - _base_center);
		std::unique_ptr<struct Intersection> isect(new Intersection{ ray, ray_isect, normal, vec2(0), _material });
		return std::move(isect);
	}
	return nullptr;
}
//https://github.com/Penetra/CG-Project/blob/master/Cone.cpp
void Cone::calculateNormal(vec3 &hitPoint, vec3 &normal){
	//PAS COMPLETER ENCORE, JAI JUSTE COPY PASTE UNE FONCTION, QUE JAI LEGEREMENT MODIFIE

	if (hitPoint.y == _base_center.y) {
		normal.x = 0;
		normal.y = -1;
		normal.z = 0;
	}
	else {
		double e = -(_radius*_radius) / (_height*_height);
		normal.x = hitPoint.x - _base_center.x;
		normal.y = hitPoint.y - _base_center.y - _height;
		normal.z = hitPoint.z - _base_center.z;
		normal.y *= e;
		normal = normalize(normal);
	}
}

int plane_Intersection(const Ray& ray, vec3 normal, vec3 point, double &t){
	double d = 0 - (dot(normal, point));
	double temp = (dot(normal, ray.direction));
	if (temp == 0){
		return 0;
	}
	t = (-d - (dot(normal, ray.origin))) / temp;
	if (t < 0.0000001)
		return 0;
	return 1;
}

vec2 calculateUVSphere(const vec3& point){
	decimal u = 0.5 + atan2(point.z, point.x) / (2 * glm::pi<decimal>());
	decimal v = 0.5 - asin(point.y) / glm::pi<decimal>();

	return vec2(u, v);
}


vec2 calculateUVCylinder(const vec3& point){
	vec2 uv;
	double phi = atan2(point.x, point.z);
	if (phi < 0.0)
		phi += 2*pi<decimal>();

	uv.x = phi * (1/(2 * pi<decimal>()));
	uv.y = (point.y + 1) / 2;
	
	
	
	return uv;
}

vec2 calculateUVCircle(const vec3& point){
	vec2 uv;
	
	uv.x = sqrt((pow(point.x, 2) + pow(point.z, 2)));
	
	double phi = (atan2(point.x, point.z));
	if (phi < 0.0){
		phi += 2 * pi<decimal>();
	}
	uv.y = phi * (1 / (2 * pi<decimal>()));
	
	return uv;
}
