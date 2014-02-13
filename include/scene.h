#pragma once
#include <main.h>
#include <unordered_map>

#include <basic_structs.h>
#include <geom.h>

class Scene
{
public:
	Scene(const char* file);
	std::unique_ptr<Intersection> trace(const Ray& ray, uint8_t depth, decimal maxdist = 1e20, decimal mindist = 1e-5) const;

	const mat4& cameraMatrix() const { return _cameraMatrix; }
	decimal fov() const { return _fov; }
	uint8_t maxDepth() const { return _maxDepth; }
	vec3 background() const { return _background; }
	decimal discretization() const { return _discretization; }

	const std::vector<std::unique_ptr<Light>>& lights() const { return _lights; }

protected:
	std::vector<std::unique_ptr<Geometry>> _geometry;
	std::vector<std::unique_ptr<Light>> _lights;
	std::unordered_map<std::string, std::unique_ptr<Material>> _materials;
	vec3 _background;

	mat4 _cameraMatrix;
	decimal _fov;

	decimal _discretization;

	uint _maxDepth;
};