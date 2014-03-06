#include <geom.h>
#include <basic_structs.h>
#include <glm/gtx/euler_angles.hpp>
#include <iostream>
#include <array>


//https://code.google.com/p/pwsraytracer/source/browse/trunk/raytracer/
Geometry::Geometry(vec3 position, vec3 orientation, vec3 scaling, Material* mtl)
:_position(position),
_orientation(orientation),
_scaling(scaling),
_material(mtl)
{	
	mat4 translation_mat = glm::translate(mat4(), position);
	mat4 scaling_mat = glm::scale(mat4(), scaling);
	mat4 rotateX = glm::rotate(mat4(), orientation.x, vec3(1,0,0));
	mat4 rotateY = glm::rotate(mat4(), orientation.y, vec3(0, 1, 0));
	mat4 rotateZ = glm::rotate(mat4(), orientation.z, vec3(0, 0, 1));
	_modelTransform ={
		translation_mat*
		
		scaling_mat*
		rotateX*rotateY*rotateZ
	};

	_inv_modelTransform = glm::inverse(_modelTransform);
}

Sphere::Sphere(vec3 position, vec3 orientation, vec3 scaling, Material* mtl)
:Geometry(position, orientation, scaling, mtl){
	_radius = 1;
	_center = vec3(0.0f);
}

std::unique_ptr<struct Intersection> Sphere::intersect(const struct Ray& ray, decimal &currentdepth) const{  
	vec3 ray_origin = vec3(_inv_modelTransform * vec4(ray.origin, 1.0f));
	vec3 ray_direction = vec3(_inv_modelTransform * vec4(ray.direction, 0.0f));	
		
	float a = glm::dot(ray_direction, ray_direction);
	float b = 2 * glm::dot(ray_origin, ray_direction);
	float c = glm::dot(ray_origin, ray_origin) - _radius * _radius;
	float t0, t1;
	if (!solveQuadratic(a, b, c, t0, t1) || t1 < 0)
		return nullptr;
	if (t1 < t0)
		std::swap(t0, t1);
	double t = (t0 < 0) ? (double)t1 : (double)t0;

	vec3 ray_isect = ray_origin + t * ray_direction;
	vec3 normal = glm::normalize(ray_isect);
	vec2 uv = calculateUVSphere(normal);
	normal = glm::normalize((vec3(glm::transpose(_inv_modelTransform) * vec4(normal, 0.0f))));
	//vec2 uv = calculateUVSphere(normal);
	ray_isect = vec3(_modelTransform * vec4(ray_isect, 1.0f));

	std::unique_ptr<struct Intersection> isect(new Intersection{ ray, ray_isect, normal, uv, _material });
	return std::move(isect);
}


Box::Box(vec3 position, vec3 orientation, vec3 scaling, Material* mtl)
:Geometry(position, orientation, scaling, mtl){
	_center = vec3(0.0f);
	_min = vec3(-1.0f);
	_max = vec3(1.0f);

	vec3 corner_ftr{ 1.0f, 1.0f, 1.0f };//0
	vec3 corner_rbl{ -1.0f, -1.0f, -1.0f };//1
	vec3 corner_fbr{ corner_ftr.x, corner_rbl.y, corner_ftr.z };//2
	vec3 corner_ftl{ corner_ftr.x, corner_ftr.y, corner_rbl.z };//3
	vec3 corner_fbl{ corner_ftr.x, corner_rbl.y, corner_rbl.z };//4
	vec3 corner_rtr{ corner_rbl.x, corner_ftr.y, corner_ftr.z };//5
	vec3 corner_rtl{ corner_rbl.x, corner_ftr.y, corner_rbl.z };//6
	vec3 corner_rbr{ corner_rbl.x, corner_rbl.y, corner_ftr.z };//7

	points = { { corner_ftr, corner_rbl, corner_fbr, corner_ftl, corner_fbl, corner_rtr, corner_rtl, corner_rbr } };

	_faces_points.push_back(vec3(1.0f, 0.0f, 0.0f));
	_faces_points.push_back(vec3(0.0f, 1.0f, 0.0f));
	_faces_points.push_back(vec3(0.0f, 0.0f, 1.0f));
	_faces_points.push_back(vec3(-1.0f, 0.0f, 0.0f));
	_faces_points.push_back(vec3(0.0f, -1.0f, 0.0f));
	_faces_points.push_back(vec3(0.0f, 0.0f, -1.0f));
}


std::unique_ptr<struct Intersection> Box::intersect(const struct Ray& ray, decimal &currentdepth) const{
	vec3 ray_origin = vec3(_inv_modelTransform * vec4(ray.origin, 1.0f));
	vec3 ray_direction = vec3(_inv_modelTransform * vec4(ray.direction, 0.0f));

	float t1 = (_min.x - ray_origin.x)*(1.0 / ray_direction.x);
	float t2 = (_max.x - ray_origin.x)*(1.0 / ray_direction.x);
	float t3 = (_min.y - ray_origin.y)*(1.0 / ray_direction.y);
	float t4 = (_max.y - ray_origin.y)*(1.0 / ray_direction.y);
	float t5 = (_min.z - ray_origin.z)*(1.0 / ray_direction.z);
	float t6 = (_max.z - ray_origin.z)*(1.0 / ray_direction.z);
	decimal t = max(max(min(t1, t2), min(t3, t4)), min(t5, t6));
	float tmax = min(min(max(t1, t2), max(t3, t4)), max(t5, t6));
	if (tmax < 0)
		return nullptr;
	if (t > tmax)
		return nullptr;

	vec3 ray_isect = ray_origin + t * ray_direction;
	
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

	vec3 normal = glm::normalize(_faces_points[index]);
	normal = glm::normalize((vec3(glm::transpose(_inv_modelTransform) * vec4(normal, 0.0f))));

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

	vec3 u_vec = glm::cross((ray_isect - uv_coord_0_0), (ray_isect - uv_coord_1_0));
	float u = glm::length(u_vec) / (2 * glm::length(uv_coord_1_0 - uv_coord_0_0));
	vec3 v_vec = glm::cross((ray_isect - uv_coord_0_0), (ray_isect - uv_coord_0_1));
	float v = glm::length(v_vec) / (2 * glm::length(uv_coord_0_1 - uv_coord_0_0));

	ray_isect = vec3(_modelTransform * vec4(ray_isect, 1.0f));
	vec2 uv = glm::vec2(u, v);
	std::unique_ptr<struct Intersection> isect(new Intersection{ ray, ray_isect, normal, uv, _material });
	return std::move(isect);
}

Cylinder::Cylinder(vec3 position, vec3 orientation, vec3 scaling, Material* mtl)
:Geometry(position, orientation, scaling, mtl){
	_center = vec3(0.0f);
	_p = vec3(0.0f, 1.0f, 0.0f);
	_q = vec3(0.0f, -1.0f, 0.0f);
	_radius = 1;
	_height = 2;
}

std::unique_ptr<struct Intersection> Cylinder::intersect(const struct Ray& ray, decimal &currentdepth) const{
	vec3 ray_origin = vec3(_inv_modelTransform * vec4(ray.origin, 1.0f));
	vec3 ray_direction = vec3(_inv_modelTransform * vec4(ray.direction, 0.0f));

https://code.google.com/p/pwsraytracer/source/browse/trunk/raytracer/cylinder.cpp?r=160

	double tmin = INFINITY;
	double t = DBL_MAX, t1 = DBL_MAX, t2 = DBL_MAX, t3 = DBL_MAX, t4 = DBL_MAX;
	bool sides = false;
	double a = pow(ray_direction.x, 2) + pow(ray_direction.z,2);
	double b = 2.0 * ((ray_origin.x - _center.x) * ray_direction.x + (ray_origin.z - _center.z) * ray_direction.z);
	double c = pow(ray_origin.x - _center.x, 2) + pow(ray_origin.z - _center.z,2) - pow(_radius,2);
	double D = b * b - 4.0 * a * c;
	if (D < 0.0) {
		return nullptr;
	}
	else{
		double sqrtD = sqrt(D);
		double aa = a * 2;
		double t1 = (-b - sqrtD) / aa;
		if (t1 > epsilon<double>()) {
			double y = (ray_origin.y - _center.y) + t1 * ray_direction.y;
			if (y > _q.y && y < _p.y) {
				if (t1 < t) {
					t = t1;
					sides = true;
				}
			}
		}
		t2 = (-b + sqrtD) / aa;
		if (t2 >  epsilon<double>()) {
			double y = (ray_origin.y - _center.y) + t2 * ray_direction.y;
			if (y > _q.y && y < _p.y) {
				if (t2 < t){
					t = t2;
					sides = true;
				}
			}
		}
	}
	if (plane_Intersection(ray_origin, ray_direction, vec3(0, -1.0, 0), _p, t3)) { 
		vec3 intersection = ray_origin + t3 * ray_direction;
		if (pow(intersection.x - _p.x, 2) + pow(intersection.z - _p.z, 2) <= _radius * _radius) { 
			if (t3 < t) {
				t = t3;
				sides = false;
			}
		}
	}
	if (plane_Intersection(ray_origin, ray_direction, vec3(0, 1.0, 0), _q, t4)) {
		vec3 intersection = ray_origin + t4 * ray_direction;
		if (pow(intersection.x - _q.x, 2) + pow(intersection.z - _q.z, 2) <= _radius * _radius) {
			if (t4 < t){
				t = t4;
				sides = false;
			}
		}
	}
	if (t != DBL_MAX && t > 0.00000001) {
		vec3 ray_isect = ray_origin + t * ray_direction;
		vec3 normal = calculateNormal(ray_isect, sides);
		normal = glm::normalize((vec3(glm::transpose(_inv_modelTransform) * vec4(normal, 0.0f))));
		vec2 uv;
		if (sides)
			uv = calculateUVCylinder(ray_isect);
		else 
			uv = calculateUVCircle(ray_isect - _p);
		ray_isect = vec3(_modelTransform * vec4(ray_isect, 1.0f));
		std::unique_ptr<struct Intersection> isect(new Intersection{ ray, ray_isect, normal, uv, _material });
		return std::move(isect);
	}
	return nullptr;
}

vec3 Cylinder::calculateNormal(vec3& hitPoint, bool sides)const{
	vec3 normal;
	if (!sides && hitPoint.y> _center.y)
		normal = normalize(_p);
	else if (!sides && hitPoint.y<_center.y)
		normal = normalize(_q);
	else {
		vec3 direction = normalize(_p - _q);
		vec3 x = _q + (dot((hitPoint - _q), direction)) * direction;
		normal = normalize(hitPoint - x);
	}
	return normal;
}


Cone::Cone(vec3 position, vec3 orientation, vec3 scaling, Material* mtl)
:Geometry(position, orientation, scaling, mtl){
	_radius = 1;
	_apex = vec3(0.0f);
	_base_center = vec3(0.0f, -1.0f, 0.0f); 
	_direction = vec3(0.0f, 1.0f, 0.0f);
	_height = 1;
	_center = vec3(0.0f);
}

std::unique_ptr<struct Intersection> Cone::intersect(const struct Ray& ray, decimal &currentdepth) const{
	vec3 ray_origin = vec3(_inv_modelTransform * vec4(ray.origin, 1.0f));
	vec3 ray_direction = vec3(_inv_modelTransform * vec4(ray.direction, 0.0f));

	// https://github.com/Penetra/CG-Project/blob/master/Cone.cpp
	double t = 0, t1 = DBL_MAX, t2 = DBL_MAX, t3 = DBL_MAX;
	bool sides = true;
	double rh = -(_radius * _radius) / (_height * _height);
	double a = pow(ray_direction.x, 2) + pow(ray_direction.z, 2) + rh * pow(ray_direction.y, 2);
	double b = 2 * (ray_direction.x * (ray_origin.x - _base_center.x) + ray_direction.z * (ray_origin.z - _base_center.z) + rh * ray_direction.y * (ray_origin.y - _base_center.y - _height));
	double c = pow(ray_origin.x - _base_center.x, 2) + pow(ray_origin.z - _base_center.z, 2) + rh * pow(ray_origin.y - _base_center.y - _height, 2);
	double root = b * b - 4.0 * a * c;
	if (root < 0)
		sides = false;
	if (sides){	
		t1 = (-b + sqrtf(root)) / (2.0 * a);
		t2 = (-b - sqrtf(root)) / (2.0 * a);
		vec3 isect = ray_origin + t1 * ray_direction;
		vec3 isect2 = ray_origin + t2 * ray_direction;
	
		if (isect.y <= _apex.y && isect.y >= _base_center.y) {
			if (isect2.y <= _apex.y && isect2.y >= _base_center.y) {
				if (t1 < t2 && t1 > 0.0 && t2 > 0.0)
					t = t1;
				else if (t2 > 0.0)
					t = t2;
			}
			else if (t1 > 0.0)
				t = t1;
		}
		else if (isect2.y <= _apex.y && isect2.y >= _base_center.y)
		if (t2 > 0.0)
			t = t2;
	}

	if (plane_Intersection(ray_origin, ray_direction, vec3(0, 1.0, 0), _base_center, t3)) {
		vec3 isect = ray_origin + t3 * ray_direction;
		if (pow(isect.x - _base_center.x, 2) + pow(isect.z - _base_center.z, 2) <= _radius * _radius) {
			if (t3 < t) {
				t = t3;
				sides = false;
			}
		}
	}

	if (t != DBL_MAX && t > 0.00000001) {
		vec3 ray_isect = ray_origin + t * ray_direction;
		vec3 normal = calculateNormal(ray_isect, sides);
		normal = glm::normalize((vec3(glm::transpose(_inv_modelTransform) * vec4(normal, 0.0f))));
		vec2 uv;
		if (sides)
			uv = calculateUVCylinder(ray_isect);
		else{	
			uv = calculateUVCircle(ray_isect - _base_center);
		}

		ray_isect = vec3(_modelTransform * vec4(ray_isect, 1.0f));
		std::unique_ptr<struct Intersection> isect(new Intersection{ ray, ray_isect, normal, uv, _material });
		return std::move(isect);
	}
	return nullptr;
}


//https://github.com/Penetra/CG-Project/blob/master/Cone.cpp
vec3 Cone::calculateNormal(vec3 &hitPoint, bool sides) const{
	//PAS COMPLETER ENCORE, JAI JUSTE COPY PASTE UNE FONCTION, QUE JAI LEGEREMENT MODIFIE
	vec3 normal;
	if (!sides)
		normal = normalize(_base_center - _apex);
	else {
		double e = -(_radius*_radius) / (_height * _height);
		normal.x = hitPoint.x - _base_center.x;
		normal.y = hitPoint.y - _base_center.y - _height;
		normal.z = hitPoint.z - _base_center.z;
		normal.y = e * normal.y;
		normal = normalize(normal);
	}
	return normal;
}

int plane_Intersection(vec3 ray_origin, vec3 ray_direction, vec3 normal, vec3 point, double &t){
	double d = 0 - (dot(normal, point));
	double temp = (dot(normal, ray_direction));
	if (temp == 0)
		return 0;
	t = (-d - (dot(normal, ray_origin))) / temp;
	if (t < 0.0000001)
		return 0;
	return 1;
}

vec2 calculateUVSphere(const vec3& point){
	//decimal u = 0.5 + atan2(point.z, point.x) / (2 * glm::pi<decimal>());
	
	float u = acos(point.y / 1.0f);
	float v = atan2(point.x, point.z);
	
	v += pi<decimal>();
	u = u / pi<decimal>();
	v = v / (2 * pi<decimal>());
	return vec2(u,v);
}

vec2 calculateUVCylinder(const vec3& point){
	vec2 uv;
	double phi = atan2(point.x, point.z);
	phi += pi<decimal>();

	uv.y = phi * (1 / (2 * pi<decimal>()));
	uv.x = (point.y + 1) / 2;

	return uv;
}


vec2 calculateUVCircle(const vec3& point){
	vec2 uv;
	uv.x = sqrt((pow(point.x, 2) + pow(point.z, 2)));
	double phi = (atan2(point.x, point.z));

	phi += pi<decimal>();

	uv.y = phi * (1 / (2 * pi<decimal>()));
	return uv;
}