#include "Primitive.h"
#include "Threading.hpp"
#include <cstring>

struct OctTree : IntersectionAccelerator {
	struct Node {
		BBox box;
		Node *children[8] = {nullptr, };
		std::vector<Intersectable*> primitives;
		bool isLeaf() const {
			return children[0] == nullptr;
		}
	};

	std::vector<Intersectable*> allPrimitives;
	Node *root = nullptr;
	int depth = 0;
	int leafSize = 0;
	int nodes = 0;
	int MAX_DEPTH = 35;
	int MIN_PRIMITIVES = 10;

	void clear(Node *n) {
		if (!n) {
			return;
		}

		for (int c = 0; c < 8; c++) {
			clear(n->children[c]);
			delete n->children[c];
		}
	}

	void clear() {
		clear(root);
		allPrimitives.clear();
	}

	void addPrimitive(Intersectable* prim) override {
		allPrimitives.push_back(prim);
	}

	void build(Node *n, int currentDepth = 0) {
		if (currentDepth >= MAX_DEPTH || n->primitives.size() <= MIN_PRIMITIVES) {
			leafSize = std::max(leafSize, int(n->primitives.size()));
			return;
		}

		depth = std::max(depth, currentDepth);
		BBox childBoxes[8];
		n->box.octSplit(childBoxes);

		for (int c = 0; c < 8; c++) {
			Node *& child = n->children[c];
			child = new Node;
			nodes++;
			memset(child->children, 0, sizeof(child->children));
			child->box = childBoxes[c];
			for (int r = 0; r < n->primitives.size(); r++) {
				if (n->primitives[r]->boxIntersect(child->box)) {
					child->primitives.push_back(n->primitives[r]);
				}
			}
			if (child->primitives.size() == n->primitives.size()) {
				build(child, MAX_DEPTH + 1);
			} else {
				build(child, currentDepth + 1);
			}
		}
		n->primitives.clear();
	}

	void build(Purpose purpose) override {
		const char *treePurpose = "";
		if (purpose == Purpose::Instances) {
			MAX_DEPTH = 5;
			MIN_PRIMITIVES = 4;
			treePurpose = " instances";
		} else if (purpose == Purpose::Mesh) {
			MAX_DEPTH = 35;
			MIN_PRIMITIVES = 20;
			treePurpose = " mesh";
		}

		if (root) {
			clear(root);
			delete root;
		}

		printf("Building%s oct tree with %d primitives... ", treePurpose, int(allPrimitives.size()));
		Timer timer;
		nodes = leafSize = depth = 0;
		root = new Node();
		root->primitives.swap(allPrimitives);
		for (int c = 0; c < root->primitives.size(); c++) {
			root->primitives[c]->expandBox(root->box);
		}
		build(root);
		printf(" done in %lldms, nodes %d, depth %d, %d leaf size\n", timer.toMs(timer.elapsedNs()), nodes, depth, leafSize);
	}

	bool intersect(Node *n, const Ray& ray, float tMin, float &tMax, Intersection& intersection) {
		bool hasHit = false;

		if (n->isLeaf()) {
			for (int c = 0; c < n->primitives.size(); c++) {
				if (n->primitives[c]->intersect(ray, tMin, tMax, intersection)) {
					tMax = intersection.t;
					hasHit = true;
				}
			}
		} else {
			for (int c = 0; c < 8; c++) {
				if (n->children[c]->box.testIntersect(ray)) {
					if (intersect(n->children[c], ray, tMin, tMax, intersection)) {
						tMax = intersection.t;
						hasHit = true;
					}
				}
			}
		}

		return hasHit;
	}

	bool intersect(const Ray& ray, float tMin, float tMax, Intersection& intersection) override {
		return intersect(root, ray, tMin, tMax, intersection);
	}

	bool isBuilt() const override {
		return root != nullptr;
	}

	~OctTree() override {
		clear();
	}
};

/// TODO: Implement one/both or any other acceleration structure and change makeDefaultAccelerator to create it
struct KDTree : IntersectionAccelerator {
	void addPrimitive(Intersectable *prim) override {}
	void clear() override {}
	void build(Purpose purpose) override {}
	bool isBuilt() const override { return false; }
	bool intersect(const Ray &ray, float tMin, float tMax, Intersection &intersection) override { return false; }
};


struct BVHTree : IntersectionAccelerator {

	struct BVHNode {
		BBox bounds;
		BVHNode* leftChild;
		BVHNode* rightChild;
		unsigned firstIndex, primitiveCount; //instead of primitive vector we will store the index of the first primitive in the array 
											 //and the count of primitives after it that are inside the box

		bool isLeaf() const {
			return !leftChild && !rightChild;
		}

		void initLeaf(unsigned firstIndex, unsigned primitiveCount, BBox bounds) {
			this->firstIndex = firstIndex;
			this->primitiveCount = primitiveCount;
			this->bounds = bounds;
			leftChild = rightChild = nullptr;
		}

		void initInterior(BVHNode* leftChild, BVHNode* rightChild) {
			this->leftChild = leftChild;
			this->rightChild = rightChild;
			bounds.add(leftChild->bounds);
			bounds.add(rightChild->bounds);

			primitiveCount = 0;
		}
	};

	struct LBVHTreelet {
		unsigned startIndex, primitivesCount;
		BVHNode* buildNodes;

		LBVHTreelet(unsigned startIndex, unsigned primCount, BVHNode* nodes) : startIndex(startIndex), primitivesCount(primCount), buildNodes(nodes)
		{}
	};

	struct MortonPrimitive {
		int primIndex;
		uint32_t mortonCode;

	public: //no need of the public but I do like to separate the functions from the data
		//Take the coordinates of a 3D point and left shift them all to match the Morton code pattern
		static uint32_t calculateMorton3D(float x, float y, float z) {
			return (LeftShift3(z) << 2) | (LeftShift3(y) << 1) | LeftShift3(x);
		}

		static void RadixSort(std::vector<MortonPrimitive>* v) {
			std::vector<MortonPrimitive> tempVector(v->size());

			constexpr int bitsPerPass = 6;
			constexpr int nBits = 30;
			constexpr int nPasses = nBits / bitsPerPass;

			for (int pass = 0; pass < nPasses; ++pass) {
				int lowBit = pass * bitsPerPass;

				std::vector<MortonPrimitive>& in = (pass & 1) ? tempVector : *v;
				std::vector<MortonPrimitive>& out = (pass & 1) ? *v : tempVector;

				constexpr int nBuckets = 1 << bitsPerPass;
				int bucketCount[nBuckets] = { 0 };
				constexpr int bitMask = (1 << bitsPerPass) - 1;

				for (const MortonPrimitive& mp : in) {
					int bucket = (mp.mortonCode >> lowBit) & bitMask;
					++bucketCount[bucket];
				}

				int outIndex[nBuckets];
				outIndex[0] = 0;

				for (int i = 1; i < nBuckets; ++i) {
					outIndex[i] = outIndex[i - 1] + bucketCount[i - 1];
				}

				for (const MortonPrimitive& mp : in) {
					int bucket = (mp.mortonCode >> lowBit) & bitMask;
					out[outIndex[bucket]++] = mp;
				}
			}

			if (nPasses & 1) {
				std::swap(*v, tempVector);
			}
		}

	private:
		//Distribute the bits of a coordinate like the masks to make them in Morton code
		static inline uint32_t LeftShift3(uint32_t x) {
			if (x == (1 << 10)) {
				--x;
			}
			x = (x | (x << 16)) & 0b00000011000000000000000011111111;
			x = (x | (x << 8)) & 0b00000011000000001111000000001111;
			x = (x | (x << 4)) & 0b00000011000011000011000011000011;
			x = (x | (x << 2)) & 0b00001001001001001001001001001001;
			return x;
		}
	};

	struct GenerateMortonCodes : Task {

		const std::vector<Intersectable*>& allPrimitives;
		const BBox& primBounds;
		std::vector<MortonPrimitive>& mortonPrims;

		GenerateMortonCodes(const std::vector<Intersectable*>& allPrimitives, const BBox& centerBounds, 
			std::vector<MortonPrimitive>& mortonPrims) : allPrimitives(allPrimitives), primBounds(centerBounds), mortonPrims(mortonPrims) {

			count = allPrimitives.size();
		}

		void run(int threadIndex, int threadCount) override {

			const unsigned perThread = count / threadCount;
			const unsigned first = threadIndex * perThread;
			const unsigned end = std::min((threadIndex + 1) * perThread, count);

			for (size_t i = first; i < end; i++) {

				initMortonFor(i, primBounds, mortonPrims);
			}
		}

	private:

		unsigned count;

		//Takes a primitive index and the bounding box of all primitives and writes the primitive's morton code inside the vector
		void initMortonFor(unsigned i, const BBox& primBounds, std::vector<MortonPrimitive>& mortonPrims) {

			//convert the centers to coordinates in range [0;1] (use Offset function)
			//scale the floating point number by 2^10 = 1024 which would make it a 10 bit number
			//store the 3 coordinates in a 32bit variable (10 bits for each coordinate and 2 bits leftover)

			constexpr int mortonBits = 10;
			constexpr int mortonScale = 1 << mortonBits; //scaling by 2^10 = 1024

			mortonPrims[i].primIndex = i; //check on that later
			vec3 centerOffset = primBounds.Offset(allPrimitives[i]->getCenter());
			vec3 offset = centerOffset * mortonScale;
			mortonPrims[i].mortonCode = MortonPrimitive::calculateMorton3D(offset.x, offset.y, offset.z);
		}
		
	};

	struct GenerateLBVHTreelet : Task {

		const std::vector<Intersectable*>& primitives;
		std::vector<LBVHTreelet>& treeletsToBuild;
		const std::vector<MortonPrimitive>& mortonPrims;
		std::vector<Intersectable*>& orderedPrims;
		std::atomic<int>* orderedPrimsOffset;
		std::atomic<int>* totalNodes;

		const unsigned MAX_PRIMS_IN_LEAF;

		GenerateLBVHTreelet(const std::vector<Intersectable*>& primitives, std::vector<LBVHTreelet>& treeletsToBuild,
							const std::vector<MortonPrimitive>& mortonPrims, std::vector<Intersectable*>& orderedPrims,
							std::atomic<int>* orderedPrimsOffset, std::atomic<int>* totalNodes, const unsigned MAX_PRIMS_IN_LEAF)
							: primitives(primitives), treeletsToBuild(treeletsToBuild), mortonPrims(mortonPrims), orderedPrims(orderedPrims), 
							  orderedPrimsOffset(orderedPrimsOffset), totalNodes(totalNodes), MAX_PRIMS_IN_LEAF(MAX_PRIMS_IN_LEAF) {
			count = treeletsToBuild.size();
		}

		void run(int threadIndex, int threadCount) override {

			const unsigned perThread = count / threadCount;
			const unsigned first = threadIndex * perThread;
			const unsigned end = std::min((threadIndex + 1) * perThread, count);

			for (size_t i = first; i < end; i++) {

				int nodesCreated = 0;
				const int firstBitIndex = 29 - 12;
				LBVHTreelet& tr = treeletsToBuild[i];

				tr.buildNodes = emitLBVH(tr.buildNodes, &mortonPrims[tr.startIndex], tr.primitivesCount, &nodesCreated, orderedPrimsOffset, firstBitIndex);
				*totalNodes += nodesCreated;
			}
		}

	private:
		unsigned count;

		BVHNode* emitLBVH(BVHNode*& buildNodes, const MortonPrimitive* mortonPrims, int primitivesCount, int* totalNodes,
						  std::atomic<int>* orderedPrimsOffset, int bitIndex) {

			//if we reach the end of the indices or the primitives are less than the max count for a leaf, we create a leaf
			if (bitIndex == -1 || primitivesCount < MAX_PRIMS_IN_LEAF) {

				++*totalNodes;
				BVHNode* node = buildNodes++;
				BBox bounds;
				int firstIndex = orderedPrimsOffset->fetch_add(primitivesCount);

				//When creating a leaf store its primitives continuous in the ordered array so they can be loaded at once in the CPU cache
				for (int i = 0; i < primitivesCount; ++i) {
					int primitiveIndex = mortonPrims[i].primIndex;
					orderedPrims[firstIndex + i] = primitives[primitiveIndex];
					primitives[primitiveIndex]->expandBox(bounds);
				}

				node->initLeaf(firstIndex, primitivesCount, bounds);
				return node;
			}
			else {
				int mask = 1 << bitIndex;

				//if all the primitives have the same bit at bitIndex, that means that all primitives are on the same side of the splitting plane
				//and we proceed without making an empty node
				if ((mortonPrims[0].mortonCode & mask) == (mortonPrims[primitivesCount - 1].mortonCode & mask)) {
					return emitLBVH(buildNodes, mortonPrims, primitivesCount, totalNodes, orderedPrimsOffset, bitIndex - 1);
				}

				//if we have primitives on both sides, search for the index where the i-th bit is different from the i+1-th bit
				unsigned splitOffset = findSplitOffset(primitivesCount, mortonPrims, mask);
				++splitOffset;

				//recursively create node for the first to splitOffset primitives and node for the rest
				(*totalNodes)++;
				BVHNode* node = buildNodes++;
				BVHNode* leftChild = emitLBVH(buildNodes, mortonPrims, splitOffset, totalNodes, orderedPrimsOffset, bitIndex - 1);
				BVHNode* rightChild = emitLBVH(buildNodes, &mortonPrims[splitOffset], primitivesCount - splitOffset, totalNodes, orderedPrimsOffset, bitIndex - 1);
				int axis = bitIndex % 3;
				node->initInterior(leftChild, rightChild); //need to include axis
				return node;
			}

		}

		//binary search to find the i-th primitive where the masked bit differes from the i+1-th primitive
		unsigned findSplitOffset(unsigned primCount, const MortonPrimitive* mortonPrims, int mask)
		{
			unsigned first = 0;
			unsigned last = primCount - 1;

			while((last - first) > 1)
			{
				unsigned half = (last - first) / 2;

				if ((mortonPrims[first].mortonCode & mask) == (mortonPrims[first + half].mortonCode & mask))
				{
					first = first + half;
				}
				else
				{
					last = first + half;
				}
			}

			return first;
		}
	};

	struct LinearBVHNode {
		BBox bounds;
		unsigned primOffset;
		unsigned secondChildOffset;

		unsigned primitiveCount;
	};

	unsigned MAX_PRIMS_IN_NODE = 20;
	std::vector<Intersectable*> allPrimitives;

	//keeping track of allocated memory, remove it if possible
	std::vector<BVHNode*> allocatedNodes;
	unsigned allocatedArrays = 0;

	BVHNode* root = nullptr;
	LinearBVHNode* nodes = nullptr;

	std::atomic<int> totalNodes{0};
	std::atomic<int> orderedPrimsOffset{0};

	void addPrimitive(Intersectable *prim) override {
		allPrimitives.push_back(prim);
	}

	void clear() override {
		for (size_t i = 0; i < allocatedArrays; i++)
		{
			delete[] allocatedNodes[i];
		}

		for (size_t i = allocatedArrays; i < allocatedNodes.size(); i++)
		{
			delete allocatedNodes[i];
		}

		delete[] nodes;
	}

	void build(Purpose purpose) override {
	
		const char* treePurpose = "";

		//get the threads count and start the thread manager that will later handle the morton codes and LBVHTreelets generation
		const int threadCount = std::max<unsigned>(std::thread::hardware_concurrency() - 1, 1);
		ThreadManager tm(threadCount);
		tm.start();


		if (purpose == Purpose::Instances) {

			//properties when building a tree for instances
			treePurpose = "Instances";
			MAX_PRIMS_IN_NODE = 5;

		} else if (purpose == Purpose::Mesh){

			//properties when building a tree for meshes
			treePurpose = "Mesh";
			MAX_PRIMS_IN_NODE = 20;
		}

		printf("Building%s BVH tree with %d primitives... ", treePurpose, int(allPrimitives.size()));
		Timer timer;

		//initialize the array that will hold the morton codes of the primitives
		unsigned primCount = allPrimitives.size();
		std::vector<Intersectable*> orderedPrims(primCount);
		std::vector<MortonPrimitive> mortonPrims(primCount);

		//get the bounding box of all centers
		BBox centersBBox;
		for (size_t i = 0; i < primCount; i++)
		{
			centersBBox.add(allPrimitives[i]->getCenter());
		}

		//Generate morton codes from the primitives using multithreading
		GenerateMortonCodes mortonTask(allPrimitives, centersBBox, mortonPrims);
		mortonTask.runOn(tm);

		//radix sort the array with morton codes
		MortonPrimitive::RadixSort(&mortonPrims);

		//storage for LBVHTreelets
		std::vector<LBVHTreelet> treeletsToBuild;
		
		//take the sorted mortion codes and traverse them linearly by taking the start and end index of primitives inside the same cluster
		for (size_t start = 0, end = 1; end <= mortonPrims.size(); ++end) {
			uint32_t mask = 0b00111111111111000000000000000000;
			if (end == (int)mortonPrims.size() || ((mortonPrims[start].mortonCode & mask) != (mortonPrims[end].mortonCode & mask))) {
				size_t primitivesCount = end - start;

				//there is upper limit for the count of treenodes based on the primitives count - 2*N - 1
				int maxBVHNodes = 2 * primitivesCount - 1;
				BVHNode* nodes = new BVHNode[maxBVHNodes];

				//for each cluster create a structure that holds the index of the first primitive in it, the number of primitives and the allocated treenodes
				treeletsToBuild.push_back(LBVHTreelet(start, primitivesCount, nodes));
				start = end;
			}
		}

		//Create LBVHs from the Treelets with multithreading
		GenerateLBVHTreelet treeletTask(allPrimitives, treeletsToBuild, mortonPrims, orderedPrims, &orderedPrimsOffset, &totalNodes, MAX_PRIMS_IN_NODE);
		treeletTask.runOn(tm);

		std::vector<BVHNode*> finishedTreelets(treeletsToBuild.size());
		for (size_t i = 0; i < treeletsToBuild.size(); i++)
		{
			finishedTreelets[i] = treeletsToBuild[i].buildNodes;
		}

		allocatedNodes = finishedTreelets;
		allocatedArrays = finishedTreelets.size();

		//after getting 16x16x16 grid with treelets finish the BVH with top to bottom SAH build
		root = buildTopToBottomSAH(finishedTreelets, 0, finishedTreelets.size() - 1);

		//make the ordered primitives the main array for primitives and remove the other one
		allPrimitives.swap(orderedPrims);
		orderedPrims.resize(0);
		orderedPrims.shrink_to_fit();

		//finally flatten the tree into array for faster traversal
		nodes = new LinearBVHNode[totalNodes];
		unsigned offset = 0;
		flattenBVH(root, offset);

		printf(" done in %lldms, nodes %d\n", timer.toMs(timer.elapsedNs()), (unsigned)totalNodes);
		tm.stop();
	}

	bool isBuilt() const override { 
		return root != nullptr;
	}

	bool intersect(const Ray &ray, float tMin, float tMax, Intersection &intersection) override { 
		//later check 7.3.5 in https://www.pbr-book.org/4ed/Primitives_and_Intersection_Acceleration/Bounding_Volume_Hierarchies#fragment-CreateleafmonoBVHBuildNode-0

		if(!nodes->bounds.testIntersect(ray)) {
			return false;
		}

		return false;
	}

private:

	struct BVHSplitBucket {
		unsigned count = 0;
		BBox bounds;
	};

	void partitionPrimitivesWithSAH(std::vector<BVHNode*>& bvhNodes, unsigned startIndx, unsigned endIndx, const BBox& centersBounds, char splitAxis, unsigned mid)
	{
		unsigned nodeCount = (endIndx - startIndx) + 1;

		if (nodeCount == 2) {
			//partitioning the elements into equally sized subsets, the nth_element works for O(n) time so it would be the same as doing a for loop
			//but the code can be reused for partitioning bigger arrays

			mid = startIndx;
			if(bvhNodes[startIndx]->bounds.getCenter()[splitAxis] > bvhNodes[endIndx]->bounds.getCenter()[splitAxis])
			{
				std::swap(bvhNodes[startIndx], bvhNodes[endIndx]);
			}

			/*mid = startIndx + nodeCount / 2;
			std::nth_element(bvhNodes[startIndx], bvhNodes[startIndx + mid], bvhNodes[endIndx], 
				[splitAxis](const BVHNode& a, const BVHNode& b) { return a.bounds.getCenter()[splitAxis] < a.bounds.getCenter()[splitAxis]; });*/
		}
		else {
			//initialize buckets with which we will check the best splitting position
			constexpr unsigned BUCKET_COUNT = 12;
			BVHSplitBucket buckets[BUCKET_COUNT];

			for (size_t i = startIndx; i < endIndx + 1; i++) {
				unsigned boxIndex = BUCKET_COUNT * centersBounds.Offset(bvhNodes[i]->bounds.getCenter())[splitAxis];
				if (boxIndex == BUCKET_COUNT) {
					boxIndex = BUCKET_COUNT - 1;
				}
				buckets[boxIndex].count++;
				buckets[boxIndex].bounds.add(bvhNodes[i]->bounds);
			}

			//get the cost for splitting after each bucket without the last one
			constexpr unsigned SPLITS = BUCKET_COUNT - 1;
			float costs[SPLITS] = {};

			unsigned countBelow = 0;
			BBox boundsBelow;

			for (size_t i = 0; i < SPLITS; i++) {
				boundsBelow.add(buckets[i].bounds);
				countBelow += buckets[i].count;
				costs[i] += countBelow * boundsBelow.getVolume();
			}

			unsigned countAbove = 0;
			BBox boundsAbove;

			for (size_t i = SPLITS; i >= 1; i--) {
				boundsAbove.add(buckets[i].bounds);
				countAbove += buckets[i].count;
				costs[i - 1] += countAbove * boundsAbove.getVolume();
			}

			//find bucket with minimum cost to split at
			int minCostBucket = -1;
			float minCost = INFINITY;

			for (size_t i = 0; i < SPLITS; i++) {
				if(costs[i] < minCost) {
					minCost = costs[i];
					minCostBucket = i;
				}
			}

			float leafCost = nodeCount;
			minCost = 1.f / 2.f + minCost / centersBounds.getVolume();

			if(nodeCount > 1) {
				auto midNode = std::partition(bvhNodes.begin() + startIndx, bvhNodes.begin() + endIndx, [=](BVHNode* node) {
					int bucket = BUCKET_COUNT * centersBounds.Offset(node->bounds.getCenter())[splitAxis];
					if (bucket == BUCKET_COUNT) bucket = BUCKET_COUNT - 1;
					return bucket <= minCostBucket;
					});
				mid = midNode - bvhNodes.begin();
			}
			//this might not get through as the nodeCount is checked before calling the function
			else {
				mid = startIndx;
			}
		}
	}

	BVHNode* buildTopToBottomSAH(std::vector<BVHNode*>& bvhNodes, unsigned startIndx, unsigned endIndx) {

		unsigned nodeCount = (endIndx - startIndx) + 1;

		//create a new node that will hold the set of nodes from startIndx to endIndx
		BVHNode* node = new BVHNode;
		allocatedNodes.push_back(node);
		++totalNodes;

		//get bounds of the nodes
		BBox bounds;
		for (size_t i = startIndx; i < endIndx + 1; i++) {
			bounds.add(bvhNodes[i]->bounds);
		}

		//if nodes size is 1 make a leaf
		if (nodeCount <= 1) {
			return bvhNodes[startIndx]; 
			//if nodeCount is 1 that would mean that start and end indices are the same, so doesn't matter which we return
		}
		else {
		//get the bounds of the node centers and choose a split axis based on the largest dimention
			BBox centerBounds;
			for (size_t i = startIndx; i < endIndx + 1; i++)
			{
				centerBounds.add(bvhNodes[i]->bounds.getCenter());
			}
			char splitAxis = centerBounds.getLongestSide();

			unsigned mid = startIndx + nodeCount / 2;

			//rearanges the bvhNodes based on the chosen mid and modifies the mid value
			partitionPrimitivesWithSAH(bvhNodes, startIndx, endIndx, centerBounds, splitAxis, mid);

			//initialize both child nodes and call parallel recursion
			BVHNode* leftChild = nullptr;
			BVHNode* rightChild = nullptr;

			leftChild = buildTopToBottomSAH(bvhNodes, startIndx, mid - 1);
			rightChild = buildTopToBottomSAH(bvhNodes, mid, endIndx);

			node->initInterior(leftChild, rightChild);
		}
		
	}

	unsigned flattenBVH(BVHNode* node, unsigned& offset) { 
		LinearBVHNode* linearNode = &nodes[offset];
		linearNode->bounds = node->bounds;
		unsigned nodeOffset = offset++;

		if(node->primitiveCount > 0) {
			linearNode->primOffset = node->firstIndex;
			linearNode->primitiveCount = node->primitiveCount;
		}
		else {
			linearNode->primitiveCount = 0;
			flattenBVH(node->leftChild, offset);
			linearNode->secondChildOffset = flattenBVH(node->rightChild, offset);
		}
		return nodeOffset;
	}

};

AcceleratorPtr makeDefaultAccelerator() {
	// TODO: uncomment or add the acceleration structure you have implemented
	//return AcceleratorPtr(new KDTree());
	return AcceleratorPtr(new BVHTree());
	//return AcceleratorPtr(new OctTree());
}