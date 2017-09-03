/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once

#include <nori/mesh.h>
#include <tbb/concurrent_vector.h>
#include <tbb/task.h>

NORI_NAMESPACE_BEGIN

/**
 * \brief Acceleration data structure for ray intersection queries
 *
 * The current implementation falls back to a brute force loop
 * through the geometry.
 */
class Accel {
public:
    /**
     * \brief Register a triangle mesh for inclusion in the acceleration
     * data structure
     *
     * This function can only be used before \ref build() is called
     */
    void addMesh(Mesh *mesh);

    /// Build the acceleration data structure
	void build();

    /// Return an axis-aligned box that bounds the scene
    const BoundingBox3f &getBoundingBox() const { return m_bbox; }

    /**
     * \brief Intersect a ray against all triangles stored in the scene and
     * return detailed intersection information
     *
     * \param ray
     *    A 3-dimensional ray data structure with minimum/maximum extent
     *    information
     *
     * \param its
     *    A detailed intersection record, which will be filled by the
     *    intersection query
     *
     * \param shadowRay
     *    \c true if this is a shadow ray query, i.e. a query that only aims to
     *    find out whether the ray is blocked or not without returning detailed
     *    intersection information.
     *
     * \return \c true if an intersection was found
     */
    bool rayIntersect(const Ray3f &ray, Intersection &its, bool shadowRay) const;

	size_t interiors() const { return m_interior; }

	size_t leaves() const { return m_leaf; }

	double averageOnLeaves() const { return 1.0 * m_total / leaves(); }

private:
	struct _Node;
    Mesh         *m_mesh = nullptr; ///< Mesh (only a single one for now)
    BoundingBox3f m_bbox;           ///< Bounding box of the entire scene
	_Node		 *m_root = nullptr; ///< Root of octree

	tbb::atomic<size_t>	  m_interior;   ///< Number of interior nodes
	tbb::atomic<size_t>	  m_leaf;		///< Number of leaf nodes
	tbb::atomic<size_t>	  m_total;		///< Number of total triangles on leaf nodes


	/// Build a node of octree
	_Node* buildTree(BoundingBox3f& box, tbb::concurrent_vector<uint32_t>& triangles, uint32_t depth);
	// Serial version
	_Node* buildTreeSerial(BoundingBox3f& box, tbb::concurrent_vector<uint32_t>& triangles, uint32_t depth);

	/// Return a 1/8 bounding box by index
	static BoundingBox3f calcBoundingBox(const BoundingBox3f& box, int index);

	/// RayIntersect variants for internal use
	bool rayIntersectInternal(const _Node* root, const BoundingBox3f& box, Ray3f &ray, Intersection &its, uint32_t &idx, bool shadowRay) const;

	static const uint32_t MAX_DEPTH = 9;
	static const size_t LEAF_SIZE = 10;
	static const size_t CUTOFF_SIZE = 80;
	static const size_t BLOCK_SIZE = 30;

	typedef struct _Node {
		union {
			struct {
				struct _Node* child[8];
				BoundingBox3f subbox[8];
			} inte;
			tbb::concurrent_vector<uint32_t>* triangles;
		};
		bool leaf;

		_Node() {

		}

		~_Node() {
			if (leaf) {
				delete triangles;
			}
			else {
				for (int i = 0; i < 8; ++i) {
					delete inte.child[i];
				}
			}
		}
	} Node;

	class BuildTask : public tbb::task {
	public:
		BuildTask(Accel& parent, Node*& root, BoundingBox3f& box, tbb::concurrent_vector<uint32_t>& triangles, uint32_t depth);
		tbb::task* execute();
	private:
		Accel& parent;
		Node*& root;
		BoundingBox3f box;
		tbb::concurrent_vector<uint32_t> triangles;
		uint32_t depth;
	};

};

NORI_NAMESPACE_END
