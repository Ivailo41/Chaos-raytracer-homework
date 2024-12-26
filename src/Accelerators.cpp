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

	unsigned MAX_PRIMS_IN_NODE = 20;

	struct BVHNode {
		BBox box;
		BVHNode* leftChild;
		BVHNode* rightChild;
		unsigned firstIndex, primitiveCount; //instead of primitive vector we will store the index of the first primitive in the array 
											 //and the count of primitives after it that are inside the box

		bool isLeaf() const {
			return !leftChild && !rightChild;
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

	struct MakeMortonCodes : Task {

		const std::vector<Intersectable*>& allPrimitives;
		const BBox& primBounds;
		std::vector<MortonPrimitive>& mortonPrims;

		MakeMortonCodes(const std::vector<Intersectable*>& allPrimitives, const BBox& centerBounds, 
			std::vector<MortonPrimitive>& mortonPrims) : allPrimitives(allPrimitives), primBounds(centerBounds), mortonPrims(mortonPrims)
		{
			count = allPrimitives.size();
		}

		void run(int threadIndex, int threadCount) override 
		{
			const unsigned perThread = count / threadCount;
			const unsigned first = threadIndex * perThread;
			const unsigned end = std::min((threadIndex + 1) * perThread, count);

			for (size_t i = first; i < end; i++)
			{
				initMortonFor(i, primBounds, mortonPrims);
			}
		}

	private:

		unsigned count;

		//Takes a primitive index and the bounding box of all primitives and writes the primitive's morton code inside the vector
		void initMortonFor(unsigned i, const BBox& primBounds, std::vector<MortonPrimitive>& mortonPrims)
		{

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

	std::vector<Intersectable*> allPrimitives;

	void addPrimitive(Intersectable *prim) override {
		allPrimitives.push_back(prim);
	}

	void clear() override {}

	void build(Purpose purpose) override {
	
		const char* treePurpose = "";

		//get the threads count and start the thread manager that will later handle the morton codes and LBVHTreelets generation
		const int threadCount = std::max<unsigned>(std::thread::hardware_concurrency() - 1, 1);
		ThreadManager tm(threadCount);
		tm.start();


		if (purpose == Purpose::Instances) {

			//properties when building a tree for instances
			treePurpose = "Instances";

		} else if (purpose == Purpose::Mesh){

			//properties when building a tree for meshes
			treePurpose = "Mesh";
			MAX_PRIMS_IN_NODE = 20;
		}

		printf("Building%s oct tree with %d primitives... ", treePurpose, int(allPrimitives.size()));
		Timer timer;

		//initialize the array that will hold the morton codes of the primitives
		unsigned primCount = allPrimitives.size();
		std::vector<unsigned> primIndices(primCount);
		std::vector<MortonPrimitive> mortonPrims(primCount);

		//initialize array of the primitive indices
		for (size_t i = 0; i < primCount; i++)
		{
			primIndices[i] = i;
		}

		//get the bounding box of all centers
		BBox centersBBox;
		for (size_t i = 0; i < primCount; i++)
		{
			centersBBox.add(allPrimitives[i]->getCenter());
		}

		//wanted to use parallel for each here but there are some issues
		//std::for_each(primIndices.begin(), primIndices.end(), [&](unsigned i) { initMortonFor(i, centersBBox, mortonPrims); });

		MakeMortonCodes mortonTask(allPrimitives, centersBBox, mortonPrims);
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

		//Create LBVHs from the Treelets (Could be parallel)
		for (size_t i = 0; i < treeletsToBuild.size(); i++)
		{
			int nodesCreated = 0;
			const int firstBitIndex = 29 - 12;
			LBVHTreelet& tr = treeletsToBuild[i];

			tr.buildNodes = emitLBVH(tr.buildNodes, allPrimitives, &mortonPrims[tr.startIndex], tr.primitivesCount, &nodesCreated, firstBitIndex);
		}

		//do multithreaded build of the LBVHTreelets

		//after getting 16x16x16 grid with treelets finish the BVH with top to bottom SAH build

		//finally flatten the tree into array for faster traversal

		//printf(" done in %lldms, nodes %d, depth %d, %d leaf size\n", timer.toMs(timer.elapsedNs()), nodes, depth, leafSize);
		tm.stop();
	}
	bool isBuilt() const override { 
		return false; //later return if the root is not nullptr
	}

	bool intersect(const Ray &ray, float tMin, float tMax, Intersection &intersection) override { 
		return false; //later check 7.3.5 in https://www.pbr-book.org/4ed/Primitives_and_Intersection_Acceleration/Bounding_Volume_Hierarchies#fragment-CreateleafmonoBVHBuildNode-0
	}

private:

	BVHNode* emitLBVH(BVHNode*& buildNodes, const std::vector<Intersectable*>& bvhPrimitives, MortonPrimitive* mortonPrims,
		int nPrimitives, int* totalNodes, int bitIndex) {

		if (bitIndex == -1 || nPrimitives < MAX_PRIMS_IN_NODE) {

			//init leaf
		}
		else {
			int mask = 1 << bitIndex;

			//check if the first and the last morton codes have the same bit which would mean that all primitives are on the same side of the splitting plane
			//if so proceed to the next bit creating a node
			if ((mortonPrims[0].mortonCode & mask) == (mortonPrims[nPrimitives - 1].mortonCode & mask))
				return emitLBVH(buildNodes, bvhPrimitives, mortonPrims, nPrimitives, totalNodes, bitIndex - 1);

			//if there are primitives on theboth sides then we do a binary search to get the point where the bit goes from 0 to 1
			/*int splitOffset = FindInterval(nPrimitives, [&](int index) {
				return ((mortonPrims[0].mortonCode & mask) == (mortonPrims[index].mortonCode & mask));
				});
			++splitOffset;*/
		}
	}
};

AcceleratorPtr makeDefaultAccelerator() {
	// TODO: uncomment or add the acceleration structure you have implemented
	//return AcceleratorPtr(new KDTree());
	return AcceleratorPtr(new BVHTree());
	//return AcceleratorPtr(new OctTree());
}