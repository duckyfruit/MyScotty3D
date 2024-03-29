
#include "shape.h"
#include "../geometry/util.h"

namespace Shapes {

Vec2 Sphere::uv(Vec3 dir) {
	float u = std::atan2(dir.z, dir.x) / (2.0f * PI_F);
	if (u < 0.0f) u += 1.0f;
	float v = std::acos(-1.0f * std::clamp(dir.y, -1.0f, 1.0f)) / PI_F;
	return Vec2{u, v};
}

BBox Sphere::bbox() const {
	BBox box;
	box.enclose(Vec3(-radius));
	box.enclose(Vec3(radius));
	return box;
}

PT::Trace Sphere::hit(Ray ray) const {
	//A3T2 - sphere hit

    // TODO (PathTracer): Task 2
    // Intersect this ray with a sphere of radius Sphere::radius centered at the origin.

    // If the ray intersects the sphere twice, ret should
    // represent the first intersection, but remember to respect
    // ray.dist_bounds! For example, if there are two intersections,
    // but only the _later_ one is within ray.dist_bounds, you should
    // return that one!

    PT::Trace ret;
    ret.origin = ray.point;

	//Vec3 s = ray.point - v_0.position;
	Vec3 d = ray.dir;
	Vec3 o = ray.point;
	Vec3 p1 = Vec3{};
	Vec3 p2 = Vec3{};
	Vec3 finp = Vec3{};
	float r = Sphere::radius;

	float tpos = 0.0f; //time to find these values!
	float tneg = 0.0f;
	float magd = sqrt(pow(d.x, 2.0f) + pow(d.y, 2.0f) + pow(d.z, 2.0f));
	float mago = sqrt(pow(o.x, 2.0f) + pow(o.y, 2.0f) + pow(o.z, 2.0f));

	if(magd == 0.0f)
	{
		ret.hit = false; 
		return ret;
	}

	float discr = 4.0f * pow(dot(o, d), 2.0f) - 4.0f * pow(magd, 2.0f) *  (pow(mago, 2.0f) - pow(r, 2.0f));
	if(discr < 0.0f)
	{
		ret.hit = false; 
		return ret;
	}

	float posnum = -2.0f * dot(o, d) + sqrt(4.0f * pow(dot(o, d), 2.0f) - 4.0f * pow(magd, 2.0f) *  (pow(mago, 2.0f) - pow(r, 2.0f)));
	float negnum = -2.0f * dot(o, d) - sqrt(4.0f * pow(dot(o, d), 2.0f) - 4.0f * pow(magd, 2.0f) *  (pow(mago, 2.0f) - pow(r, 2.0f)));
	float denom = 2.0f * pow(magd, 2.0f);

	tpos = posnum / denom;
	tneg = negnum / denom;

	p1 = o + tpos * d;
	p2 = o + tneg * d;
	
	finp = p1;

	float p1dist = sqrt(pow(p1.x - ray.point.x,2.0f) + pow(p1.y - ray.point.y,2.0f) + pow(p1.z - ray.point.z,2.0f));
	float p2dist = sqrt(pow(p2.x - ray.point.x,2.0f) + pow(p2.y - ray.point.y,2.0f) + pow(p2.z - ray.point.z,2.0f));

	float fpdist = p1dist;

	if(p1dist > p2dist)
	{
		if((p2dist < ray.dist_bounds.x) || (p2dist > ray.dist_bounds.y))
		{
			ret.hit = false;
			return ret;	
		}
		finp = p2;
		fpdist = p2dist;
	}
	else
	{
		if((p1dist < ray.dist_bounds.x) || (p1dist > ray.dist_bounds.y))
		{
			ret.hit = false;
			return ret;	
		}
	}
	
	/*if((p1dist < ray.dist_bounds.x) || (p1dist > ray.dist_bounds.y))
	{
		if((p2dist < ray.dist_bounds.x) || (p2dist > ray.dist_bounds.y))
		{
			ret.hit = false;
			return ret;	
		}
		finp = p2;
		fpdist = p2dist;
	} */
	
	ret.hit = true;       // was there an intersection?
    ret.distance = fpdist;   // at what distance did the intersection occur?
    ret.position = ray.at(fpdist); // where was the intersection?
    ret.normal = finp.unit();   // what was the surface normal at the intersection?
	ret.uv = Sphere::uv(ret.normal); 	   // what was the uv coordinates at the intersection? (you may find Sphere::uv to be useful)
    return ret;

	/*
    ret.hit = false;       // was there an intersection?
    ret.distance = 0.0f;   // at what distance did the intersection occur?
    ret.position = Vec3{}; // where was the intersection?
    ret.normal = Vec3{};   // what was the surface normal at the intersection?
	ret.uv = Vec2{}; 	   // what was the uv coordinates at the intersection? (you may find Sphere::uv to be useful)
    return ret; */
}

Vec3 Sphere::sample(RNG &rng, Vec3 from) const {
	die("Sampling sphere area lights is not implemented yet.");
}

float Sphere::pdf(Ray ray, Mat4 pdf_T, Mat4 pdf_iT) const {
	die("Sampling sphere area lights is not implemented yet.");
}

Indexed_Mesh Sphere::to_mesh() const {
	return Util::closed_sphere_mesh(radius, 2);
}

} // namespace Shapes

bool operator!=(const Shapes::Sphere& a, const Shapes::Sphere& b) {
	return a.radius != b.radius;
}

bool operator!=(const Shape& a, const Shape& b) {
	if (a.shape.index() != b.shape.index()) return false;
	return std::visit(
		[&](const auto& shape) {
			return shape != std::get<std::decay_t<decltype(shape)>>(b.shape);
		},
		a.shape);
}
