#pragma once
#include <main.h>
#include <material.h>

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
	float _radius;
	vec3 _center;
	
};

class Box : public Geometry
{
public:
	Box(vec3 position, vec3 orientation, vec3 scaling, Material* mtl = new Material());

	virtual std::unique_ptr<struct Intersection> intersect(const struct Ray& ray, decimal &currentdepth) const override;
};

class Cylinder : public Geometry
{
public:
	Cylinder(vec3 position, vec3 orientation, vec3 scaling, Material* mtl = new Material());

	virtual std::unique_ptr<struct Intersection> intersect(const struct Ray& ray, decimal &currentdepth) const override;
};

class Cone : public Geometry
{
public:
	Cone(vec3 position, vec3 orientation, vec3 scaling, Material* mtl = new Material());

	virtual std::unique_ptr<struct Intersection> intersect(const struct Ray& ray, decimal &currentdepth) const override;
};