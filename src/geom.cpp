#include <geom.h>
#include <basic_structs.h>
#include <glm/gtx/euler_angles.hpp>
#include <iostream>
#include <array>

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
	//std::cout << "intersection at " << position.x << " " << position.y << " " << position.z << std::endl;
	vec3 normal = glm::normalize(ray_isect - _center);

	// Calcul des coordonees uv
	decimal x = ray_isect.x;
	decimal y = ray_isect.y;
	decimal z = ray_isect.z;
	
	decimal u = 0.5 + atan2(normal.z, normal.x) / (2 * glm::pi<decimal>());
	decimal v = 0.5 - asin(normal.y) / glm::pi<decimal>();

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
	//calcul de la normale
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
	_center = vec3(_modelTransform*vec4(0.0,0,0,1));
	_p = vec3(_modelTransform*vec4(0.0, 1, 0,1));
	_q = vec3(_modelTransform*vec4(0.0, -1, 0, 1));
	_height = scaling.y*2;
	_radius = 1*scaling.x;

	
	//TODO implementer constructeur cylindre
}

std::unique_ptr<struct Intersection> Cylinder::intersect(const struct Ray& ray, decimal &currentdepth) const{

https://code.google.com/p/pwsraytracer/source/browse/trunk/raytracer/cylinder.cpp?r=160

	//ps cest normal que la lumiere du cylindre soit fuckee, pcq il faut le fermer.
	//ce qui se passe cest que le rayon touche a linterieur le fond, ensuite, avec le bias sur la normale, 
	//il reussi quand mm a pogner la lumiere
	double tmin = INFINITY;
	double t = DBL_MAX, t1 = DBL_MAX, t2 = DBL_MAX, t3 = DBL_MAX;
	vec3 normal;
	double ttemp = 0;
	
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
		double t1 = (-b - sqrtD) / aa;

		if (t1 > epsilon<double>())
		{
			double y = oy + t1 * dy;
			if (y > _q.y && y < _p.y)
			{
				
				if (t1 < tmin) {
					t = t1;
					//WARNING
					//attention au calcul de la normal, il va falloir changer le y
					normal = normalize(vec3(((ox + dx*t) * (1 / _radius)), 0, ((oz + dz*t) * (1 / _radius))));
				}
			}
		}
		t2 = (-b + sqrtD) / aa;
		if (t2 >  epsilon<double>())
		{
			double y = oy + t2 * dy;
			if (y > _q.y && y < _p.y)
			{
				
				if (t2 < t){
					t = t2;
					//WARNING
					//attention au calcul de la normal, il va falloir changer le y
					normal = normalize(vec3(((ox + dx*t) * (1 / _radius)), 0, ((oz + dz*t) * (1 / _radius))));
				}
			}
		}
	}
	if (plane_Intersection(ray, vec3(0, -1, 0), _p, t3)) { /* intersects plane */
		vec3 intersection = ray.origin + (ray.direction*t3);
		if (pow(intersection.x - _p.x, 2) + pow(intersection.z - _p.z, 2) <= _radius*_radius) { /* intersects Cone */
			if (t3 < t){
				t = t3;
				normal = vec3(0, 1, 0);
			}
		}
	}

	if (t != DBL_MAX && t>0.00000001){
		vec3 ray_isect = ray.origin + ray.direction*(decimal)t;
		std::unique_ptr<struct Intersection> isect(new Intersection{ ray, ray_isect, normal, vec2(0), _material });
		return std::move(isect);
	}


	return nullptr;
}

Cone::Cone(vec3 position, vec3 orientation, vec3 scaling, Material* mtl)
:Geometry(position, orientation, scaling, mtl){
	_radius = 1.0;
	_apex = vec3(0.0, 2, 0);
	_base_center = vec3(0.0, 0, 0); 
	_direction = normalize(_base_center - _apex);
	_theta = atan(_radius / length(_base_center - _apex));
	_height = length(_base_center - _apex);
	
	std::cout << _theta;
}
inline bool Quadratic(float A, float B, float C, float *t0,
	float *t1) {
	// Find quadratic discriminant
	float discrim = B * B - 4.f * A * C;
	if (discrim < 0.) return false;
	float rootDiscrim = sqrtf(discrim);
	// Compute quadratic _t_ values
	float q;
	if (B < 0) q = -.5f * (B - rootDiscrim);
	else       q = -.5f * (B + rootDiscrim);
	*t0 = q / A;
	*t1 = C / q;
	if (*t0 > *t1) std::swap(*t0, *t1);
	return true;

}
std::unique_ptr<struct Intersection> Cone::intersect(const struct Ray& ray, decimal &currentdepth) const{
    //https://github.com/Penetra/CG-Project/blob/master/Cone.cpp

	double rh = -(_radius*_radius) / (_height*_height);
	double a = pow(ray.direction.x, 2) + pow(ray.direction.z, 2) + rh*pow(ray.direction.y, 2);
	double b = 2 * (ray.direction.x*(ray.origin.x - _base_center.x) + ray.direction.z*(ray.origin.z - _base_center.z) + rh*ray.direction.y*(ray.origin.y - _base_center.y - _height));
	double c = pow(ray.origin.x - _base_center.x, 2) + pow(ray.origin.z - _base_center.z, 2) + rh*pow(ray.origin.y - _base_center.y -_height, 2);
	double root = b*b - 4.0*a*c;
	if (root<0)
		return nullptr;
	double t = DBL_MAX, t1 = DBL_MAX, t2 = DBL_MAX, t3 = DBL_MAX;
	t1 = (-b + sqrtf(root)) / (2.0*a);
	t2 = (-b - sqrtf(root)) / (2.0*a);

	vec3 intersection = (ray.origin) + ((ray.direction) * t1);
	vec3 intersection2 = ray.origin + ((ray.direction) *t2);

	if (intersection.y <= _apex.y && intersection.y >= _base_center.y) {
		if (intersection2.y <= _apex.y && intersection2.y >= _base_center.y) {
			if (t1<t2 && t1>0 && t2 >0)
				t = t1;
			else {
				if (t2>0.0f)
					t = t2;
			}
		}
		else {
			if (t1>0.0f)
				t = t1;
		}
	}
	else {
		if (intersection2.y <= _apex.y && intersection2.y >=_base_center.y) {
			if (t2 > 0.0f)
				t = t2;
		}
	}
	

	if (plane_Intersection(ray, vec3(0,-1,0), _base_center, t3)) {
		intersection = ray.origin + (ray.direction*t3);
		if (pow(intersection.x - _base_center.x, 2) + pow(intersection.z - _base_center.z, 2) <= _radius*_radius) {
			if (t3<t)
				t = t3;
		}
	}

	if (t != DBL_MAX && t>0.00000001){
		vec3 ray_isect = ray.origin + ray.direction*(decimal)t;
		vec3 normal = normalize(ray_isect - _base_center);
		std::unique_ptr<struct Intersection> isect(new Intersection{ ray, ray_isect, normal, vec2(0), _material });
		return std::move(isect);
	}
	return nullptr;
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