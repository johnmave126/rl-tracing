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

#include <nori/accel.h>
#include <Eigen/Geometry>
#include <algorithm>

NORI_NAMESPACE_BEGIN

void Accel::addMesh(Mesh *mesh) {
    if (m_mesh)
        throw NoriException("Accel: only a single mesh is supported!");
    m_mesh = mesh;
    m_bbox = m_mesh->getBoundingBox();
}

void Accel::build() {
	std::vector<uint32_t> all;
	for (uint32_t idx = 0; idx < m_mesh->getTriangleCount(); ++idx) {
		all.push_back(idx);
	}
	m_root = buildTree(m_bbox, all, 0);
}

bool Accel::rayIntersect(const Ray3f &ray_, Intersection &its, bool shadowRay) const {
    bool foundIntersection = false;  // Was an intersection found so far?
    uint32_t f = (uint32_t) -1;      // Triangle index of the closest intersection

    Ray3f ray(ray_); /// Make a copy of the ray (we will need to update its '.maxt' value)

    /* Delegate to octree to find an intersection */
	foundIntersection = rayIntersectInternal(m_root, m_bbox, ray, its, f, shadowRay);

    if (foundIntersection) {
        /* At this point, we now know that there is an intersection,
           and we know the triangle index of the closest such intersection.

           The following computes a number of additional properties which
           characterize the intersection (normals, texture coordinates, etc..)
        */

        /* Find the barycentric coordinates */
        Vector3f bary;
        bary << 1-its.uv.sum(), its.uv;

        /* References to all relevant mesh buffers */
        const Mesh *mesh   = its.mesh;
        const MatrixXf &V  = mesh->getVertexPositions();
        const MatrixXf &N  = mesh->getVertexNormals();
        const MatrixXf &UV = mesh->getVertexTexCoords();
        const MatrixXu &F  = mesh->getIndices();

        /* Vertex indices of the triangle */
        uint32_t idx0 = F(0, f), idx1 = F(1, f), idx2 = F(2, f);

        Point3f p0 = V.col(idx0), p1 = V.col(idx1), p2 = V.col(idx2);

        /* Compute the intersection positon accurately
           using barycentric coordinates */
        its.p = bary.x() * p0 + bary.y() * p1 + bary.z() * p2;

        /* Compute proper texture coordinates if provided by the mesh */
        if (UV.size() > 0)
            its.uv = bary.x() * UV.col(idx0) +
                bary.y() * UV.col(idx1) +
                bary.z() * UV.col(idx2);

        /* Compute the geometry frame */
        its.geoFrame = Frame((p1-p0).cross(p2-p0).normalized());

        if (N.size() > 0) {
            /* Compute the shading frame. Note that for simplicity,
               the current implementation doesn't attempt to provide
               tangents that are continuous across the surface. That
               means that this code will need to be modified to be able
               use anisotropic BRDFs, which need tangent continuity */

            its.shFrame = Frame(
                (bary.x() * N.col(idx0) +
                 bary.y() * N.col(idx1) +
                 bary.z() * N.col(idx2)).normalized());
        } else {
            its.shFrame = its.geoFrame;
        }
    }

    return foundIntersection;
}

Accel::Node* Accel::buildTree(BoundingBox3f& box, std::vector<uint32_t>& triangles, uint32_t depth) {
	Node *node = nullptr;
	if (triangles.size() > 0) {
		if (triangles.size() <= LEAF_SIZE || depth > MAX_DEPTH) {
			m_leaf++;
			m_total += triangles.size();
			node = new Node;
			node->leaf = true;
			node->triangles = new std::vector<uint32_t>(std::move(triangles));
		}
		else {
			m_interior++;
			node = new Node;
			node->leaf = false;
			std::vector<uint32_t> candidates[8];
			for (int i = 0; i < 8; ++i) {
				node->inte.subbox[i] = calcBoundingBox(box, i);
				for (std::vector<uint32_t>::iterator it = triangles.begin(); it != triangles.end(); ++it) {
					if (node->inte.subbox[i].overlaps(m_mesh->getBoundingBox(*it))) {
						candidates[i].push_back(*it);
					}
				}
			}
			for (int i = 0; i < 8; ++i) {
				node->inte.child[i] = buildTree(node->inte.subbox[i], candidates[i], depth + 1);
			}
		}
	}
	return node;
}

BoundingBox3f Accel::calcBoundingBox(const BoundingBox3f& box, int index) {
	Point3f min, max;
	for (int i = 0; i < 3; ++i) {
		if (index & (1<<i)) {
			min[i] = box.min[i];
			max[i] = (box.min[i] + box.max[i]) / 2;
		}
		else {
			min[i] = (box.min[i] + box.max[i]) / 2;
			max[i] = box.max[i];
		}
	}
	return BoundingBox3f(min, max);
}

/* Any order traversal version */
//bool Accel::rayIntersectInternal(const Accel::Node* root, const BoundingBox3f& box, Ray3f &ray, Intersection &its, uint32_t &idx, bool shadowRay) const {
//	if (!root) {
//		return false;
//	}
//	bool foundIntersection = false;  // Was an intersection found so far?
//	if (root->leaf) {
//		//Test all triangles in the leaf node
//		for (std::vector<uint32_t>::iterator it = root->triangles->begin(); it != root->triangles->end(); ++it) {
//			float u, v, t;
//			if (m_mesh->rayIntersect(*it, ray, u, v, t)) {
//				/* An intersection was found! Can terminate
//				immediately if this is a shadow ray query */
//				if (shadowRay)
//					return true;
//				ray.maxt = its.t = t;
//				its.uv = Point2f(u, v);
//				its.mesh = m_mesh;
//				idx = *it;
//				foundIntersection = true;
//			}
//		}
//	}
//	else {
//		for (int i = 0; i < 8; ++i) {
//			if (root->inte.subbox[i].rayIntersect(ray)) {
//				foundIntersection = rayIntersectInternal(root->inte.child[i], root->inte.subbox[i], ray, its, idx, shadowRay) || foundIntersection;
//				if (shadowRay && foundIntersection) {
//					break;
//				}
//			}
//		}
//	}
//	return foundIntersection;
//}

/* Near to far order traversal version */
bool Accel::rayIntersectInternal(const Accel::Node* root, const BoundingBox3f& box, Ray3f &ray, Intersection &its, uint32_t &idx, bool shadowRay) const {
	bool foundIntersection = false;  // Was an intersection found so far?
	if (root->leaf) {
		//Test all triangles in the leaf node
		for (std::vector<uint32_t>::iterator it = root->triangles->begin(); it != root->triangles->end(); ++it) {
			float u, v, t;
			if (m_mesh->rayIntersect(*it, ray, u, v, t)) {
				/* An intersection was found! Can terminate
				immediately if this is a shadow ray query */
				if (shadowRay)
					return true;
				ray.maxt = its.t = t;
				its.uv = Point2f(u, v);
				its.mesh = m_mesh;
				idx = *it;
				foundIntersection = true;
			}
		}
	}
	else {
		std::vector<std::pair<float, int> > candidate;
		for (int i = 0; i < 8; ++i) {
			if (root->inte.child[i]) {
				float nt, ft;
				if (root->inte.subbox[i].rayIntersect(ray, nt, ft)) {
					candidate.push_back(std::make_pair(nt, i));
					if (candidate.size() > 3) {
						break;
					}
				}
			}
		}
		if (candidate.size() > 0) {
			std::sort(candidate.begin(), candidate.end());
			for (int i = 0; i < candidate.size() && !foundIntersection; ++i) {
				foundIntersection = rayIntersectInternal(root->inte.child[candidate[i].second], root->inte.subbox[candidate[i].second], ray, its, idx, shadowRay);
			}
		}
	}
	return foundIntersection;
}

NORI_NAMESPACE_END

