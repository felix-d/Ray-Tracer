#pragma once
#include <main.h>
#include <material.h>
#include <array>
#include "../vs12/rtbase/utils.h"

int plane_Intersection(const Ray& ray, vec3 normal, vec3 point, double &t);
vec2 calculateUVSphere(const vec3& point);
vec2 calculateUVCylinder(const vec3& point);
vec2 calculateUVCircle(const vec3& point);

class Geometry
{
public:
	Geometry(vec3 position, vec3 orientation, vec3 scaling, Material* mtl = new Material());

	// An interesting use of this function definition is to only create an intersection
	// object when we are certain that this new intersection would be the closest.
	// Passing currentdepth (the linear distance between the ray's origin and the previous closest object)
	// lets us do just that. Moreover, the & lets us update that value straight in the intersect method.
	virtual std::unique_ptr<struct Intersection> intersect(const struct Ray& ray, decimal &currentdepth) const abstract;
	
protected:
	vec3 _position;
	vec3 _orientation;
	vec3 _scaling;
	mat4 _modelTransform;
	Material* _material;

	// Transform order: scaling, then rotation, then translation (use glm methods)
	// Hint: store the transforms in this class as a single matrix
	// Hint: preprocess any modified matrices you might need (like the inverse)
};

class Sphere : public Geometry
{
public:
	Sphere(vec3 position, vec3 orientation, vec3 scaling, Material* mtl = new Material());

	//const after a function declaration means that the function is not allowed to 
	//change any class members(except ones that are marked mutable)
	virtual std::unique_ptr<struct Intersection> intersect(const struct Ray& ray, decimal &currentdepth) const override;
protected:
	decimal _radius;
	vec3 _center;
	
};

class Box : public Geometry
{
public:
	Box(vec3 position, vec3 orientation, vec3 scaling, Material* mtl = new Material());

	virtual std::unique_ptr<struct Intersection> intersect(const struct Ray& ray, decimal &currentdepth) const override;
protected:
	vec3 _min;
	vec3 _max;
	vec3 _center;
	std::array<vec3, 8> points;
	std::vector<vec3>_faces_points;
	

};

class Cylinder : public Geometry
{
public:
	Cylinder(vec3 position, vec3 orientation, vec3 scaling, Material* mtl = new Material());

	virtual std::unique_ptr<struct Intersection> intersect(const struct Ray& ray, decimal &currentdepth) const override;
	vec3 calculateNormal(vec3& hitPoint, bool sides) const;
protected:
	vec3 _center;
	vec3 _p;
	vec3 _q;
	decimal _radius;
	decimal _height;
;
};

class Cone : public Geometry
{
public:
	Cone(vec3 position, vec3 orientation, vec3 scaling, Material* mtl = new Material());
	void Cone::calculateNormal(vec3& hitPoint, vec3& normal);
	virtual std::unique_ptr<struct Intersection> intersect(const struct Ray& ray, decimal &currentdepth) const override;
protected:
	vec3 _base_center;
	vec3 _apex;
	vec3 _direction;
	decimal _radius;
	decimal _height;
	
};

