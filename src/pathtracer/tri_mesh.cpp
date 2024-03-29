
#include "../test.h"

#include "samplers.h"
#include "tri_mesh.h"

namespace PT {

BBox Triangle::bbox() const {
	//A3T2 / A3T3
	BBox box;
	Tri_Mesh_Vert v_0 = vertex_list[v0];
    Tri_Mesh_Vert v_1 = vertex_list[v1];
    Tri_Mesh_Vert v_2 = vertex_list[v2];

	Vec3 p1 = v_0.position;
	Vec3 p2 = v_1.position;
	Vec3 p3 = v_2.position;

	box.enclose(p1);
	box.enclose(p2);
	box.enclose(p3);
	/*box.min.x = p1.x;
	box.min.y = p1.y;
	box.min.z = p1.z;

	float xdistp2 = sqrt(pow(p2.x - p1.x,2.0f));
	float xdistp3 = sqrt(pow(p3.x - p1.x,2.0f));


	float ydistp2 = sqrt(pow(p2.y - p1.y,2.0f));
	float ydistp3 = sqrt(pow(p3.y - p1.y,2.0f));


	float zdistp2 = sqrt(pow(p2.z - p1.z,2.0f));
	float zdistp3 = sqrt(pow(p3.z - p1.z,2.0f));

	if(xdistp2 > xdistp3)
	box.max.x = p2.x;
	else 
	box.max.x = p3.x;

	if(ydistp2 > ydistp3)
	box.max.y = p2.y;
	else 
	box.max.y = p3.y;

	if(zdistp2 > zdistp3)
	box.max.z = p2.z;
	else 
	box.max.z = p3.z;

	float temp = 0.0f;

	if(box.max.x < box.min.x)
	{
		temp = box.min.x;
		box.min.x = box.max.x;
		box.max.x = temp;
	}
	if(box.max.y < box.min.y)
	{
		temp = box.min.y;
		box.min.y = box.max.y;
		box.max.y = temp;
	}
	if(box.max.z < box.min.z)
	{
		temp = box.min.z;
		box.min.z = box.max.z;
		box.max.z = temp;
	}
	 */
	
	
	// TODO (PathTracer): Task 2 or 3
    // Compute the bounding box of the triangle.

    // Beware of flat/zero-volume boxes! You may need to
    // account for that here, or later on in BBox::hit.

    
    return box;
}

Trace Triangle::hit(const Ray& ray) const {
	//A3T2
	
	// Each vertex contains a postion and surface normal
    Tri_Mesh_Vert v_0 = vertex_list[v0];
    Tri_Mesh_Vert v_1 = vertex_list[v1];
    Tri_Mesh_Vert v_2 = vertex_list[v2];
    (void)v_0;
    (void)v_1;
    (void)v_2;
	Trace ret;
    ret.origin = ray.point;

    // TODO (PathTracer): Task 2
    // Intersect the ray with the triangle defined by the three vertices.

	//first, we can find out if the ray intersects with the formula from rasterization!
	//we just need to get the point 
	//find out if the ray intersected! With the formula given!
	//barycentric coordinates must sum to 1
	//barycentric coordinates must be positive

	Vec3 s = ray.point - v_0.position;
	Vec3 d = ray.dir;
	Vec3 e1 = v_1.position - v_0.position;
	Vec3 e2 =  v_2.position - v_0.position;
	
	
	float u = 0.0f;
	float v = 0.0f;
	float w = 0.0f;
	float t = 0.0f; //time to find these values!

	
	float denom = dot(cross(e1,d),e2);
	if(denom == 0.0f) //determinant is 0, triangle has no area
	{
		ret.hit = false;
		return ret;
	}

	denom = 1/denom;
	u = denom * (dot(-1.0f * (cross(s,e2)), d));
	v = denom * (dot((cross(e1,d)), s));
	t = denom * (dot(-1.0f * (cross(s,e2)), e1));

	w = 1 - u - v;

	Vec3 P = ray.point + t * ray.dir;

	float dist = sqrt(pow(P.x - ray.point.x,2.0f) + pow(P.y - ray.point.y,2.0f) + pow(P.z - ray.point.z,2.0f));
	if((dist < ray.dist_bounds.x) || (dist > ray.dist_bounds.y))
	{
		ret.hit = false;
		return ret;		
	}

   /* ret.hit = false;       // was there an intersection?
    ret.distance = 0.0f;   // at what distance did the intersection occur?
    ret.position = Vec3{}; // where was the intersection?
    ret.normal = Vec3{};   // what was the surface normal at the intersection?
                           // (this should be interpolated between the three vertex normals)
	ret.uv = Vec2{};	   // What was the uv associated with the point of intersection?
						   // (this should be interpolated between the three vertex uvs) */
	
	if((u < 0.0f) || (v < 0.0f )||( w < 0.0f))
	{
		ret.hit = false;
    	return ret;
	}
	
	ret.hit = true;
	ret.distance = dist;
	ret.position = ray.at(dist);
	ret.normal = v_0.normal * w + v_1.normal * u + v_2.normal * v;
	ret.uv = v_0.uv * w + v_1.uv * u + v_2.uv * v;

	return ret;
}

Triangle::Triangle(Tri_Mesh_Vert* verts, uint32_t v0, uint32_t v1, uint32_t v2)
	: v0(v0), v1(v1), v2(v2), vertex_list(verts) {
}

Vec3 Triangle::sample(RNG &rng, Vec3 from) const {
	Tri_Mesh_Vert v_0 = vertex_list[v0];
	Tri_Mesh_Vert v_1 = vertex_list[v1];
	Tri_Mesh_Vert v_2 = vertex_list[v2];
	Samplers::Triangle sampler(v_0.position, v_1.position, v_2.position);
	Vec3 pos = sampler.sample(rng);
	return (pos - from).unit();
}

float Triangle::pdf(Ray wray, const Mat4& T, const Mat4& iT) const {

	Ray tray = wray;
	tray.transform(iT);

	Trace trace = hit(tray);
	if (trace.hit) {
		trace.transform(T, iT.T());
		Vec3 v_0 = T * vertex_list[v0].position;
		Vec3 v_1 = T * vertex_list[v1].position;
		Vec3 v_2 = T * vertex_list[v2].position;
		Samplers::Triangle sampler(v_0, v_1, v_2);
		float a = sampler.pdf(trace.position);
		float g = (trace.position - wray.point).norm_squared() / std::abs(dot(trace.normal, wray.dir));
		return a * g;
	}
	return 0.0f;
}

bool Triangle::operator==(const Triangle& rhs) const {
	if (Test::differs(vertex_list[v0].position, rhs.vertex_list[rhs.v0].position) ||
	    Test::differs(vertex_list[v0].normal, rhs.vertex_list[rhs.v0].normal) ||
	    Test::differs(vertex_list[v0].uv, rhs.vertex_list[rhs.v0].uv) ||
	    Test::differs(vertex_list[v1].position, rhs.vertex_list[rhs.v1].position) ||
	    Test::differs(vertex_list[v1].normal, rhs.vertex_list[rhs.v1].normal) ||
	    Test::differs(vertex_list[v1].uv, rhs.vertex_list[rhs.v1].uv) ||
	    Test::differs(vertex_list[v2].position, rhs.vertex_list[rhs.v2].position) ||
	    Test::differs(vertex_list[v2].normal, rhs.vertex_list[rhs.v2].normal) ||
	    Test::differs(vertex_list[v2].uv, rhs.vertex_list[rhs.v2].uv)) {
		return false;
	}
	return true;
}

Tri_Mesh::Tri_Mesh(const Indexed_Mesh& mesh, bool use_bvh_) : use_bvh(use_bvh_) {
	for (const auto& v : mesh.vertices()) {
		verts.push_back({v.pos, v.norm, v.uv});
	}

	const auto& idxs = mesh.indices();

	std::vector<Triangle> tris;
	for (size_t i = 0; i < idxs.size(); i += 3) {
		tris.push_back(Triangle(verts.data(), idxs[i], idxs[i + 1], idxs[i + 2]));
	}

	if (use_bvh) {
		triangle_bvh.build(std::move(tris), 4);
	} else {
		triangle_list = List<Triangle>(std::move(tris));
	}
}

Tri_Mesh Tri_Mesh::copy() const {
	Tri_Mesh ret;
	ret.verts = verts;
	ret.triangle_bvh = triangle_bvh.copy();
	ret.triangle_list = triangle_list.copy();
	ret.use_bvh = use_bvh;
	return ret;
}

BBox Tri_Mesh::bbox() const {
	if (use_bvh) return triangle_bvh.bbox();
	return triangle_list.bbox();
}

Trace Tri_Mesh::hit(const Ray& ray) const {
	if (use_bvh) return triangle_bvh.hit(ray);
	return triangle_list.hit(ray);
}

size_t Tri_Mesh::n_triangles() const {
	return use_bvh ? triangle_bvh.n_primitives() : triangle_list.n_primitives();
}

uint32_t Tri_Mesh::visualize(GL::Lines& lines, GL::Lines& active, uint32_t level,
                             const Mat4& trans) const {
	if (use_bvh) return triangle_bvh.visualize(lines, active, level, trans);
	return 0u;
}

Vec3 Tri_Mesh::sample(RNG &rng, Vec3 from) const {
	if (use_bvh) {
		return triangle_bvh.sample(rng, from);
	}
	return triangle_list.sample(rng, from);
}

float Tri_Mesh::pdf(Ray ray, const Mat4& T, const Mat4& iT) const {
	if (use_bvh) {
		return triangle_bvh.pdf(ray, T, iT);
	}
	return triangle_list.pdf(ray, T, iT);
}

} // namespace PT
