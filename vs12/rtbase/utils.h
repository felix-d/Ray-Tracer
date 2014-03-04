int planeIntersection(const Ray& ray, vec3 normal, vec3 point, double &t){
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