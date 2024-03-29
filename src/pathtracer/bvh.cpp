
#include "bvh.h"
#include "aggregate.h"
#include "instance.h"
#include "tri_mesh.h"

#include <stack>

namespace PT {

struct BVHBuildData {
	BVHBuildData(size_t start, size_t range, size_t dst) : start(start), range(range), node(dst) {
	}
	size_t start; ///< start index into the primitive array
	size_t range; ///< range of index into the primitive array
	size_t node;  ///< address to update
};

struct SAHBucketData {
	BBox bb;          ///< bbox of all primitives
	size_t num_prims; ///< number of primitives in the bucket
};

template<typename Primitive>
void BVH<Primitive>::build(std::vector<Primitive>&& prims, size_t max_leaf_size) {
	//A3T3 - build a bvh

	// Keep these
    
	nodes.clear(); //clear the nodes
    primitives = std::move(prims); //move the primitives
	BBox start;

	
	printf(" %f ",primitives[0].bbox().center().z );
	printf(" %f ",primitives[1].bbox().center().z );
	printf(" %f ",primitives[2].bbox().center().z );
	printf(" %f ",primitives[3].bbox().center().z );
	for(int x = 0; x< primitives.size(); x++)
	{
		start.enclose(primitives[x].bbox());
	}

	//new_node(start, 0, primitives.size(),1,2);
	// HELPER FUNCTION -----------------------------
	

	// HELPER FUNCTION ----------------------

	for(int a = 0; a < 3; a++)
	{

		makeTree(start, 0, primitives.size(), 0, max_leaf_size, a);
	}

}

template<typename Primitive>
size_t BVH<Primitive>::makeTree(BBox box, size_t start, size_t size, size_t leaf_num, size_t max_leaf_size, int axis)
	{
		int numBuckets = 10; //number of buckets
		float bucketSize = 0.0f; //size of buckets
		float len = 0.0f; //center of the current node bbox
		
		std::vector<SAHBucketData> buckets = std::vector<SAHBucketData>{10}; //vector of buckets
		std::vector<size_t> startingPrim = std::vector<size_t>{10}; //vector of starting primitives for buckets

		float SAHCost = 0.0f; //all the SAH variables 
		float sA = 0.0f; //surface area
		float sB = 0.0f;
		size_t nA = 0; //num primitives
		size_t nB = 0;
		float bestSAH = FLT_MAX; //best SAH 
		int SAHIndexL = 0; //index of the best SAH
		int SAHIndexR = 0; //index of the best SAH

		if(leaf_num == max_leaf_size)
		{ //return new_node(box, start, size,0,0); 
		 return leaf_num;
		}

		if(axis == 0 ) //x axis
		len = box.max.x - box.min.x;
		else if(axis == 1) //yaxis
		len = box.max.y - box.min.y;
		else //zaxis
		len = box.max.z - box.min.z;

		printf(" %f ", len);
		if(int(len) == 0)
		return new_node(box, start, size, leaf_num , leaf_num); //make new leaf node based on best SAH

		

		bucketSize = float((len)/numBuckets); //get the size of each buckets
		
		for(size_t i = start; i < size; i ++) //range of primitives
		{
			
			float buc = 0;
			if(axis == 0 )
			buc = ((primitives[i].bbox().center().x) /bucketSize);
			else if(axis == 1)
			buc = ((primitives[i].bbox().center().y) /bucketSize);
			else
			buc = ((primitives[i].bbox().center().z) /bucketSize);

			//printf(" %f ", buc);

			buckets[int(buc)].bb.enclose(primitives[i].bbox()); //enclose the primitive bbox in the bucket bbox

			if(buckets[int(buc)].num_prims == 0)
			{startingPrim[int(buc)] = i;}
			
			//printf(" %f ", buckets[int(buc)].bb.center().z);

			buckets[int(buc)].num_prims += 1; //increase num primitives

			if(buckets[int(buc)].num_prims == primitives.size())
			{
				return new_node(buckets[int(buc)].bb, start, size, leaf_num , leaf_num); //make new leaf node based on best SAH
			}
			

		}
		
		for(int s = 0; s < numBuckets; s++) //conduct SAH on each partition of bucket 
		{
			for(int b = 0; b < numBuckets; b++) //conduct SAH on each partition of bucket 
			{
				if((s < b) && !(buckets[s].bb.empty() || buckets[b].bb.empty())  ) //s will always be left side, b will always be right side
				{ 
					sA = buckets[s].bb.surface_area();
					sB = buckets[b].bb.surface_area();
					nA = buckets[s].num_prims;
					nB = buckets[b].num_prims;
					//printf(" %f ", sB);
					SAHCost = sA * float(nA) + sB * float(nB);
					if((SAHCost < bestSAH) && (SAHCost > 0.0f))
					{
						//printf("%f", SAHCost);
						bestSAH = SAHCost;
						SAHIndexL = s;
						SAHIndexR = b;
					}

		
					//SANA SBNB
					//conduct SAH equation 
					//update SAHIndex with the best SAH 
					//update indicies 
				}
			}
		}
		printf(" %d %d ", SAHIndexL, SAHIndexR);
		
		size_t indl = makeTree(buckets[SAHIndexL].bb, startingPrim[SAHIndexL],buckets[SAHIndexL].num_prims, leaf_num + 1, max_leaf_size,axis);
		size_t indr = makeTree(buckets[SAHIndexR].bb, startingPrim[SAHIndexR],buckets[SAHIndexR].num_prims, leaf_num + 1, max_leaf_size, axis);
		
		return new_node(box,start, size, indl, indr); //make new leaf node based on best SAH
}

template<typename Primitive> Trace BVH<Primitive>::hit(const Ray& ray) const {
	//A3T3 - traverse your BVH

    // Implement ray - BVH intersection test. A ray intersects
    // with a BVH aggregate if and only if it intersects a primitive in
    // the BVH that is not an aggregate.

    // The starter code simply iterates through all the primitives.
    // Again, remember you can use hit() on any Primitive value.

	//TODO: replace this code with a more efficient traversal:
    Trace ret;
    for(const Primitive& prim : primitives) {
        Trace hit = prim.hit(ray);
        ret = Trace::min(ret, hit);
    }
    return ret;
}

template<typename Primitive>
BVH<Primitive>::BVH(std::vector<Primitive>&& prims, size_t max_leaf_size) {
	build(std::move(prims), max_leaf_size);
}

template<typename Primitive> std::vector<Primitive> BVH<Primitive>::destructure() {
	nodes.clear();
	return std::move(primitives);
}

template<typename Primitive>
template<typename P>
typename std::enable_if<std::is_copy_assignable_v<P>, BVH<P>>::type BVH<Primitive>::copy() const {
	BVH<Primitive> ret;
	ret.nodes = nodes;
	ret.primitives = primitives;
	ret.root_idx = root_idx;
	return ret;
}

template<typename Primitive> Vec3 BVH<Primitive>::sample(RNG &rng, Vec3 from) const {
	if (primitives.empty()) return {};
	int32_t n = rng.integer(0, static_cast<int32_t>(primitives.size()));
	return primitives[n].sample(rng, from);
}

template<typename Primitive>
float BVH<Primitive>::pdf(Ray ray, const Mat4& T, const Mat4& iT) const {
	if (primitives.empty()) return 0.0f;
	float ret = 0.0f;
	for (auto& prim : primitives) ret += prim.pdf(ray, T, iT);
	return ret / primitives.size();
}

template<typename Primitive> void BVH<Primitive>::clear() {
	nodes.clear();
	primitives.clear();
}

template<typename Primitive> bool BVH<Primitive>::Node::is_leaf() const {
	// A node is a leaf if l == r, since all interior nodes must have distinct children
	return l == r;
}

template<typename Primitive>
size_t BVH<Primitive>::new_node(BBox box, size_t start, size_t size, size_t l, size_t r) {
	Node n;
	n.bbox = box;
	n.start = start;
	n.size = size;
	n.l = l;
	n.r = r;
	nodes.push_back(n);
	return nodes.size() - 1;
}
 
template<typename Primitive> BBox BVH<Primitive>::bbox() const {
	if(nodes.empty()) return BBox{Vec3{0.0f}, Vec3{0.0f}};
	return nodes[root_idx].bbox;
}

template<typename Primitive> size_t BVH<Primitive>::n_primitives() const {
	return primitives.size();
}

template<typename Primitive>
uint32_t BVH<Primitive>::visualize(GL::Lines& lines, GL::Lines& active, uint32_t level,
                                   const Mat4& trans) const {

	std::stack<std::pair<size_t, uint32_t>> tstack;
	tstack.push({root_idx, 0u});
	uint32_t max_level = 0u;

	if (nodes.empty()) return max_level;

	while (!tstack.empty()) {

		auto [idx, lvl] = tstack.top();
		max_level = std::max(max_level, lvl);
		const Node& node = nodes[idx];
		tstack.pop();

		Spectrum color = lvl == level ? Spectrum(1.0f, 0.0f, 0.0f) : Spectrum(1.0f);
		GL::Lines& add = lvl == level ? active : lines;

		BBox box = node.bbox;
		box.transform(trans);
		Vec3 min = box.min, max = box.max;

		auto edge = [&](Vec3 a, Vec3 b) { add.add(a, b, color); };

		edge(min, Vec3{max.x, min.y, min.z});
		edge(min, Vec3{min.x, max.y, min.z});
		edge(min, Vec3{min.x, min.y, max.z});
		edge(max, Vec3{min.x, max.y, max.z});
		edge(max, Vec3{max.x, min.y, max.z});
		edge(max, Vec3{max.x, max.y, min.z});
		edge(Vec3{min.x, max.y, min.z}, Vec3{max.x, max.y, min.z});
		edge(Vec3{min.x, max.y, min.z}, Vec3{min.x, max.y, max.z});
		edge(Vec3{min.x, min.y, max.z}, Vec3{max.x, min.y, max.z});
		edge(Vec3{min.x, min.y, max.z}, Vec3{min.x, max.y, max.z});
		edge(Vec3{max.x, min.y, min.z}, Vec3{max.x, max.y, min.z});
		edge(Vec3{max.x, min.y, min.z}, Vec3{max.x, min.y, max.z});

		if (!node.is_leaf()) {
			tstack.push({node.l, lvl + 1});
			tstack.push({node.r, lvl + 1});
		} else {
			for (size_t i = node.start; i < node.start + node.size; i++) {
				uint32_t c = primitives[i].visualize(lines, active, level - lvl, trans);
				max_level = std::max(c + lvl, max_level);
			}
		}
	}
	return max_level;
}

template class BVH<Triangle>;
template class BVH<Instance>;
template class BVH<Aggregate>;
template BVH<Triangle> BVH<Triangle>::copy<Triangle>() const;

} // namespace PT
