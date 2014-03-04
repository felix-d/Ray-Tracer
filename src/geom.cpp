#include <geom.h>
#include <basic_structs.h>
#include <glm/gtx/euler_angles.hpp>
#include <iostream>
#include <array>
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
		scaling_mat*
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

	vec3 corner_ftr{ 1, 1, 1 };//0
	vec3 corner_rbl{-1, -1, -1 };//1
	vec3 corner_fbr{ corner_ftr.x, corner_rbl.y, corner_ftr.z };//2
	vec3 corner_ftl{ corner_ftr.x, corner_ftr.y, corner_rbl.z };//3
	vec3 corner_fbl{ corner_ftr.x, corner_rbl.y, corner_rbl.z };//4
	vec3 corner_rtr{ corner_rbl.x, corner_ftr.y, corner_ftr.z };//5
	vec3 corner_rtl{ corner_rbl.x, corner_ftr.y, corner_rbl.z };//6
	vec3 corner_rbr{ corner_rbl.x, corner_rbl.y, corner_ftr.z };//7
	init_points = { { corner_ftr, corner_rbl, corner_fbr, corner_ftl, corner_fbl, corner_rtr, corner_rtl, corner_rbr } };
	for (int i = 0; i < 8; i++){
		vec4 v4(init_points[i], 1);
		v4 = (_modelTransform*v4);
		vec3 v3(v4);
		points[i] = v3;
	}
	SetCenter();
	SetNormals();
	SetExtents();

	_faces_points.push_back(vec3(1.0f, 0.0f, 0.0f));
	_faces_points.push_back(vec3(0.0f, 1.0f, 0.0f));
	_faces_points.push_back(vec3(0.0f, 0.0f, 1.0f));
	_faces_points.push_back(vec3(-1.0f, 0.0f, 0.0f));
	_faces_points.push_back(vec3(0.0f, -1.0f, 0.0f));
	_faces_points.push_back(vec3(0.0f, 0.0f, -1.0f));

	for (uint i = 0; i < _faces_points.size(); i++){
		vec4 current = vec4(_faces_points[i], 1.0f);
		current =  _modelTransform*current;
		_faces_points[i] = vec3(current);
	}
	
}

void Box::SetCenter(){
	auto resultX = std::minmax_element(points.begin(), points.end(), [](const vec3& lhs, const vec3& rhs) {
		return lhs.x < rhs.x;
	});
	auto resultY = std::minmax_element(points.begin(), points.end(), [](const vec3& lhs, const vec3& rhs) {
		return lhs.y < rhs.y;
	});
	auto resultZ = std::minmax_element(points.begin(), points.end(), [](const vec3& lhs, const vec3& rhs) {
		return lhs.z < rhs.z;
	});
	_center.x = resultX.first->x + (resultX.second->x - resultX.first->x) / 2;
	_center.y = resultY.first->y + (resultY.second->y - resultY.first->y) / 2;
	_center.z = resultZ.first->z + (resultZ.second->z - resultZ.first->z) / 2;
}
void Box::SetNormals(){
	u.at(0) = normalize(points.at(4) - points.at(1)); //axe des x
	u.at(1) = normalize(points.at(6) - points.at(1)); //axe des y
	u.at(2) = normalize(points.at(7) - points.at(1)); //axe des z
}
void Box::SetExtents(){
	e.at(0) = glm::length(points.at(4) - points.at(1)) / 2; //fbl-rbl
	e.at(1) = glm::length(points.at(6) - points.at(1)) / 2; //axe des y
	e.at(2) = glm::length(points.at(7) - points.at(1)) / 2; //axe des z
}

std::unique_ptr<struct Intersection> Box::intersect(const struct Ray& ray, decimal &currentdepth) const{
	float maxS = -FLT_MAX;
	float minT = FLT_MAX;


	vec3 diff = _center - ray.origin;

	for (int i = 0; i < 3; ++i){
		vec3 axis = u[i];
		float et = dot(axis,diff);
		float f = dot(ray.direction,axis);
		if (f==0){
			if (-et - e[i] > 0.0f || -et + e[i] > 0.0f)
				return false;
			continue;
		}
		float s = (et - e[i]) / f;
		float t = (et + e[i]) / f;
		if (s > t){
			float temp = s;
			s = t;
			t = temp;
		}
		if (s > maxS)
			maxS = s;
		if (t < minT)
			minT = t;
		if (minT < 0.0f || maxS > minT)
			return false;
	}
	vec3 ray_isect = ray.origin + ray.direction* (decimal)maxS;
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
	std::unique_ptr<struct Intersection> isect(new Intersection{ ray, ray_isect, normal, vec2(0), _material });
	return std::move(isect);
}

Cylinder::Cylinder(vec3 position, vec3 orientation, vec3 scaling, Material* mtl)
:Geometry(position, orientation, scaling, mtl){
	_center = vec3(0,0,0);
	_p = vec3(0.0,1,0);
	_q = vec3(0.0, -1, 0);
	_height = 2;
	_radius = 1;
	mY0 = -1;
	mY1 = 1;
	
	//TODO implementer constructeur cylindre
}

std::unique_ptr<struct Intersection> Cylinder::intersect(const struct Ray& ray, decimal &currentdepth) const{

https://code.google.com/p/pwsraytracer/source/browse/trunk/raytracer/cylinder.cpp?r=160
	double tmin = INFINITY;
	double t;
	double ox = ray.origin.x - _center.x;
	double oy = ray.origin.y - _center.y;
	double oz = ray.origin.z - _center.z;

	double dx = ray.direction.x;
	double dy = ray.direction.y;
	double dz = ray.direction.z;

	double a = dx * dx + dz * dz;
	double b = 2.0 * (ox * dx + oz * dz);
	double c = ox * ox + oz * oz - _radius * _radius;
	double D = b * b - 4.0 * a * c;

	if (D < 0.0)
	{
		//No hitpoints
		return nullptr;
	}
	else
	{
		double sqrtD = sqrt(D);
		double aa = a * 2;
		double t = (-b - sqrtD) / aa;

		if (t > epsilon<double>())
		{
			double y = oy + t * dy;
			if (y > mY0 && y < mY1)
			{
				
				 vec3 ray_isect = ray.origin + ray.direction*(decimal)t;
				 std::unique_ptr<struct Intersection> isect(new Intersection{ ray, ray_isect, vec3(0), vec2(0), _material });
				 return std::move(isect);
				
			}
		}

		t = (-b + sqrtD) / aa;
		if (t >  epsilon<double>())
		{
			double y = oy + t * dy;
			if (y > mY0 && y < mY1)
			{
				
				vec3 ray_isect = ray.origin + ray.direction*(decimal)t;
				std::unique_ptr<struct Intersection> isect(new Intersection{ ray, ray_isect, vec3(0), vec2(0), _material });
				return std::move(isect);
			}
		}
	}
	return nullptr;
}

Cone::Cone(vec3 position, vec3 orientation, vec3 scaling, Material* mtl)
:Geometry(position, orientation, scaling, mtl){
	//TODO implementer constructeur cone
}

std::unique_ptr<struct Intersection> Cone::intersect(const struct Ray& ray, decimal &currentdepth) const{
	//TODO implementer intersection entre cone et rayon.
	return NULL;
}

