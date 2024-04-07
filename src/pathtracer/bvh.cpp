
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
	/*BBox start;

	std::vector<BBox> boxes;

	
	printf(" %f ",primitives[0].bbox().center().z );
	printf(" %f ",primitives[1].bbox().center().z );
	printf(" %f ",primitives[2].bbox().center().z );
	printf(" %f ",primitives[3].bbox().center().z ); */
	/*for(int x = 0; x< primitives.size(); x++)
	{
		start.enclose(primitives[x].bbox());
		boxes.push_back(primitives[x].bbox());
	}

	//new_node(start, 0, primitives.size(),1,2);
	// HELPER FUNCTION -----------------------------
	

	// HELPER FUNCTION ----------------------


	makeTree(start,primitives.size(), boxes, 0, max_leaf_size);

	for(int i = 0; i < nodes.size(); i++)
	{
		printf("%d\n", i);
		printf(" %f %f %f \n", nodes[i].bbox.min.x, nodes[i].bbox.min.y, nodes[i].bbox.min.z );
	} */
	
}

template<typename Primitive>
size_t BVH<Primitive>::makeTree(BBox box, size_t size, std::vector<BBox> boxes, size_t leaf_num, size_t max_leaf_size)
	{
		size_t root = new_node(); //node we are working on!
		
		if((leaf_num == max_leaf_size) || (size <= 1))
		{ 
			printf("%zu %zu\n", root, size);
			nodes[root].bbox = box;
			nodes[root].start = leaf_num;
			nodes[root].size = size;
			nodes[root].l = nodes.size() - 1;
			nodes[root].r = nodes.size() - 1;
			
			return root;
		}

		size_t numBuckets = 10; //number of buckets
		float bucketSize = 0.0f; //size of buckets
		float len = 0.0f; //center of the current node bbox
		int bestAxis = 0; //leave it as X for now

		float sA = 0.0f; //surface area
		float sB = 0.0f;
		size_t nA = 0; //num primitives
		size_t nB = 0;
		float bestSAH = FLT_MAX; //best SAH 
		float SAHCost = 0.0f; //all the SAH variables 

		BBox partitionL;
		BBox partitionR;
		size_t sizeL = 0;
		size_t sizeR = 0;
		size_t indl = 0;
		size_t indr = 0;


		std::vector<SAHBucketData> bucketsX = std::vector<SAHBucketData>{numBuckets}; //vector of buckets
		std::vector<SAHBucketData> bucketsY = std::vector<SAHBucketData>{numBuckets}; //vector of buckets
		std::vector<SAHBucketData> bucketsZ = std::vector<SAHBucketData>{numBuckets}; //vector of buckets

		std::vector<std::vector<BBox>> mapX(numBuckets); //maps of primitives
		std::vector<std::vector<BBox>> mapY(numBuckets); //maps of primitives
		std::vector<std::vector<BBox>> mapZ(numBuckets); //maps of primitives

		std::vector<BBox> primL; //left primitives
		std::vector<BBox> primR; //right primitives

		for(size_t i = 0; i < numBuckets; i++) //create the 2d array of bboxes
		{
			mapX[i] = std::vector<BBox>(size); 
			mapY[i] = std::vector<BBox>(size); 
			mapZ[i] = std::vector<BBox>(size); 
		}

		for(int axis = 0; axis < 3; axis ++)
		{
			if(axis == 0 ) //x axis
			len = box.max.x - box.min.x;
			else if(axis == 1) //yaxis
			len = box.max.y - box.min.y;
			else //zaxis
			len = box.max.z - box.min.z;
			
			if(!(int(len) == 0))
			{

					bucketSize = float((len)/(numBuckets - 1)); //get the size of each buckets
					//printf("%f\n",bucketSize);
					for(size_t i = 0; i < size; i ++) //range of primitives
					{
						int buc = 0;

						if(axis == 0 )
						{

							buc = int((boxes[i].center().x) /bucketSize);
							bucketsX[buc].bb.enclose(boxes[i]); //enclose the primitive bbox in the bucket bbox
							mapX[buc][bucketsX[buc].num_prims] = boxes[i];
							bucketsX[buc].num_prims += 1; //increase num primitives
							//printf("%d\n",buc);
						}
						else if(axis == 1)
						{
							buc = int((boxes[i].center().y) /bucketSize);
							bucketsY[buc].bb.enclose(boxes[i]); //enclose the primitive bbox in the bucket bbox
							mapY[buc][bucketsY[buc].num_prims] = boxes[i];
							bucketsY[buc].num_prims += 1; //increase num primitives
						}
						else
						{
							buc = int((boxes[i].center().z) /bucketSize);
							bucketsZ[buc].bb.enclose(boxes[i]); //enclose the primitive bbox in the bucket bbox
							mapZ[buc][bucketsZ[buc].num_prims] = boxes[i];
							bucketsZ[buc].num_prims += 1; //increase num primitives
						}

					}
				
				
				for(int s = 1; s < numBuckets - 1; s++) //conduct SAH on each partition of bucket 
				{
					sA = 0.0f;
					sB = 0.0f;
					nA = 0;
					nB = 0;
					if(axis == 0)
					{
						for(int x = 0; x < s; x ++)
						{
							sA += bucketsX[x].bb.surface_area();
							nA += bucketsX[x].num_prims;
						}
					
						for(int y = s; y < numBuckets; y ++)
						{
							sB += bucketsX[y].bb.surface_area();
							nB += bucketsX[y].num_prims;
						}
					
						SAHCost = sA * float(nA) + sB * float(nB);
						if((SAHCost < bestSAH) && (SAHCost != 0) )
						{
							//printf("%f %f %zu %zu\n", sA, sB, nA, nB);
							bestAxis = axis;
							bestSAH = SAHCost;
							partitionL.reset();
							sizeL = 0;
							primL.clear();
							partitionR.reset();
							sizeR = 0;
							primR.clear();
							for(int x = 0; x < s; x ++)
							{
								if(!bucketsX[x].bb.empty())
								{
									partitionL.enclose(bucketsX[x].bb);
									sizeL += bucketsX[x].num_prims;

									//printf("%f %f\n",bucketsX[x].bb.min.x, bucketsX[x].bb.max.x);

									for(int i = 0; i < mapX[x].size(); i++)
									{
										if(!mapX[x][i].empty())
										primL.push_back(mapX[x][i]);
									}
									
								}
							}
						
							for(int y = s; y < numBuckets; y ++)
							{
								if(!bucketsX[y].bb.empty())
								{
									partitionR.enclose(bucketsX[y].bb);
									sizeR += bucketsX[y].num_prims;
									//printf("%f %f\n",bucketsX[y].bb.min.x, bucketsX[y].bb.max.x);
									for(int i = 0; i < mapX[y].size(); i++)
									{
										if(!mapX[y][i].empty())
										primR.push_back(mapX[y][i]);
									}
								}
							}
					
						}
					}
					else if(axis == 1)
					{
						for(int x = 0; x < s; x ++)
						{
							sA += bucketsY[x].bb.surface_area();
							nA += bucketsY[x].num_prims;
						}
						
						for(int y = s; y < numBuckets; y ++)
						{
							sB += bucketsY[y].bb.surface_area();
							nB += bucketsY[y].num_prims;
						}
					
						SAHCost = sA * float(nA) + sB * float(nB);
						if((SAHCost < bestSAH) && (SAHCost != 0) )
						{
							
							bestAxis = axis;
							bestSAH = SAHCost;

							partitionL.reset();
							sizeL = 0;
							primL.clear();
							partitionR.reset();
							sizeR = 0;
							primR.clear();
							for(int x = 0; x < s; x ++)
							{
								if(!bucketsY[x].bb.empty())
								{
									partitionL.enclose(bucketsY[x].bb);
									sizeL += bucketsY[x].num_prims;

									for(int i = 0; i < mapY[x].size(); i++)
									{
										if(!mapY[x][i].empty())
										primL.push_back(mapY[x][i]);
									}
								}
							}
							for(int y = s; y < numBuckets; y ++)
							{
								if(!bucketsY[y].bb.empty())
								{
									partitionR.enclose(bucketsY[y].bb);
									sizeR += bucketsY[y].num_prims;

									for(int i = 0; i < mapY[y].size(); i++)
									{
										if(!mapY[y][i].empty())
										primR.push_back(mapY[y][i]);
									}
								}
							}
						}
					}
					else
					{
						for(int x = 0; x < s; x ++)
						{
							sA += bucketsZ[x].bb.surface_area();
							nA += bucketsZ[x].num_prims;
						}
						for(int y = s; y < numBuckets; y ++)
						{
							sB += bucketsZ[y].bb.surface_area();
							nB += bucketsZ[y].num_prims;
						}
						SAHCost = sA * float(nA) + sB * float(nB);
						
						if((SAHCost < bestSAH) && (SAHCost != 0) )
						{
							
							bestAxis = axis;
							bestSAH = SAHCost;

							partitionL.reset();
							sizeL = 0;
							primL.clear();
							partitionR.reset();
							sizeR = 0;
							primR.clear();
							for(int x = 0; x < s; x ++)
							{
								if(!bucketsZ[x].bb.empty())
								{
									partitionL.enclose(bucketsZ[x].bb);
									sizeL += bucketsZ[x].num_prims;

									for(int i = 0; i < mapZ[x].size(); i++)
									{
										if(!mapZ[x][i].empty())
										primL.push_back(mapZ[x][i]);
									}
								}
							}
							for(int y = s; y < numBuckets; y ++)
							{
								if(!bucketsZ[y].bb.empty())
								{
									partitionR.enclose(bucketsZ[y].bb);
									sizeR += bucketsZ[y].num_prims;

									for(int i = 0; i < mapZ[y].size(); i++)
									{
										if(!mapZ[y][i].empty())
										primR.push_back(mapZ[y][i]);
									}
								}
							}
						}
					}
				}
			
			}
		}

		//sort primitives at the end?
		for(size_t i = 0; i < sizeL; i ++)
		{
			primitives[i].bbox() = primL[i];
		}
		for(size_t i = sizeL; i < sizeL + sizeR; i ++)
		{
			primitives[i].bbox() = primR[i - sizeL];
		}
		

		if(bestAxis == 0)
		{
			indl = makeTree(partitionL, sizeL, primL, leaf_num + 1, max_leaf_size);
			indr = makeTree(partitionR, sizeR, primR, leaf_num + 2, max_leaf_size);
		}
		else if (bestAxis == 1)
		{
			indl = makeTree(partitionL, sizeL, primL,leaf_num + 1, max_leaf_size);
			indr = makeTree(partitionR, sizeR, primR,leaf_num + 2, max_leaf_size);
		}
		else
		{
			indl = makeTree(partitionL, sizeL, primL,leaf_num + 1, max_leaf_size);
			indr = makeTree(partitionR, sizeR, primR,leaf_num + 2, max_leaf_size);
		}

		printf("%zu %zu\n", root, size);
		nodes[root].bbox = box;
		nodes[root].start = 0;
		nodes[root].size = size;
		nodes[root].l = indl;
		nodes[root].r = indr;

		//printf("%zu %zu", indl, indr);
		//size_t node =  new_node(box, 0, size, indl, indr); //make new leaf node based on best SAH
		//printf(" %zu\n ", root);

		

		return root;
		
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
