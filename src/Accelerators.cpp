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
		BBox box;
		BVHNode* leftChild;
		BVHNode* rightChild;
		unsigned firstIndex, primitiveCount; //instead of primitive vector we will store the index of the first primitive in the array 
											 //and the count of primitives after it that are inside the box

		bool isLeaf() const {
			return !leftChild && !rightChild;
		}
	};

	std::vector<Intersectable*> allPrimitives;

	void addPrimitive(Intersectable *prim) override {
		allPrimitives.push_back(prim);
	}

	void clear() override {}

	void build(Purpose purpose) override {
	
		if (purpose == Purpose::Instances) {

			//properties when building a tree for instances

		} else if (purpose == Purpose::Mesh){

			//properties when building a tree for meshes
		}

		//take the center points of the primitives and store them in array (we have the number of our primitives so dont use a dynamic array)
		//convert the centers to coordinates in range [0;1] (use Offset function)
		//scale the floating point number by 2^10 = 1024 which would make it a 10 bit number
		//store the 3 coordinates in a 32bit variable (10 bits for each coordinate and 2 bits leftover)
		//radix sort the array with morton codes

		//build LBVHTreelets
		//take the sorted mortion codes and traverse them linearly by taking the start and end index of primitives inside the same cluster
		//for each cluster create a structure that holds the index of the first primitive in it, the number of primitives and the allocated treenodes
		//there is upper limit for the count of treenodes based on the primitives count - 2*N - 1
		//do multithreaded build of the LBVHTreelets
		//after getting 16x16x16 grid with treelets finish the BVH with top to bottom SAH build
		//finally flatten the tree into array for faster traversal

		//gives offset of a point in the range [0;1] in a bounding box (put it inside the BBox class
		/*Vector3f Offset(Point3f p) const {
			Vector3f o = p - pMin;
			if (pMax.x > pMin.x) o.x /= pMax.x - pMin.x;
			if (pMax.y > pMin.y) o.y /= pMax.y - pMin.y;
			if (pMax.z > pMin.z) o.z /= pMax.z - pMin.z;
			return o;
		}*/
	}
	bool isBuilt() const override { 
		return false; //later return if the root is not nullptr
	}

	bool intersect(const Ray &ray, float tMin, float tMax, Intersection &intersection) override { 
		return false; //later check 7.3.5 in https://www.pbr-book.org/4ed/Primitives_and_Intersection_Acceleration/Bounding_Volume_Hierarchies#fragment-CreateleafmonoBVHBuildNode-0
	}
};

AcceleratorPtr makeDefaultAccelerator() {
	// TODO: uncomment or add the acceleration structure you have implemented
	//return AcceleratorPtr(new KDTree());
	//return AcceleratorPtr(new BVHTree());
	return AcceleratorPtr(new OctTree());
}

