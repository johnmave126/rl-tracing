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

#include <tracer/accel.h>

NORI_NAMESPACE_BEGIN

/**
 * \brief Main scene data structure
 *
 * This class holds information on scene objects and is responsible for
 * coordinating rendering jobs. It also provides useful query routines that
 * are mostly used by the \ref Integrator implementations.
 */
class Scene : public NoriObject {
public:
    /// Construct a new scene object
    Scene(const PropertyList & props);

    /// Release all memory
    virtual ~Scene();

    /// Return a pointer to the scene's kd-tree
    const Accel *getAccel() const { return m_accel; }

    /// Return a pointer to the scene's integrator
    const Integrator *getIntegrator() const { return m_integrator; }

    /// Return a pointer to the scene's integrator
    Integrator *getIntegrator() { return m_integrator; }

    /// Return a pointer to the scene's camera
    const Camera *getCamera() const { return m_camera; }

    /// Return a pointer to the scene's sample generator (const version)
    const Sampler *getSampler() const { return m_sampler; }

    /// Return a pointer to the scene's sample generator
    Sampler *getSampler() { return m_sampler; }

    /// Return a reference to an array containing all meshes
    const std::vector<Mesh *> &getMeshes() const { return m_meshes; }

    /// Return a reference to an array containing all emitting meshes
    const std::vector<Emitter *> &getEmitters() const { return m_emitters; }

    /// Return whether the rendering should be progressive
    bool isProgressive() const { return m_isprogressive; }

    /**
     * \brief Intersect a ray against all triangles stored in the scene
     * and return detailed intersection information
     *
     * \param ray
     *    A 3-dimensional ray data structure with minimum/maximum
     *    extent information
     *
     * \param its
     *    A detailed intersection record, which will be filled by the
     *    intersection query
     *
     * \return \c true if an intersection was found
     */
    bool rayIntersect(const Ray3f &ray, Intersection &its) const {
        return m_accel->rayIntersect(ray, its, false);
    }

    /**
     * \brief Intersect a ray against all triangles stored in the scene
     * and \a only determine whether or not there is an intersection.
     *
     * This method much faster than the other ray tracing function,
     * but the performance comes at the cost of not providing any
     * additional information about the detected intersection
     * (not even its position).
     *
     * \param ray
     *    A 3-dimensional ray data structure with minimum/maximum
     *    extent information
     *
     * \return \c true if an intersection was found
     */
    bool rayIntersect(const Ray3f &ray) const {
        Intersection its; /* Unused */
        return m_accel->rayIntersect(ray, its, true);
    }

    /// \brief Return an axis-aligned box that bounds the scene
    const BoundingBox3f &getBoundingBox() const {
        return m_accel->getBoundingBox();
    }

    /// Take a [0,1]^2 uniform sample and return a random sample over all emitters
    const Emitter* sampleEmitter(float& sample, float &pdf) const;

	const Emitter* sampleEmitter(const float& sample, float &pdf) const;

    /**
     * \brief Inherited from \ref NoriObject::activate()
     *
     * Initializes the internal data structures (kd-tree,
     * emitter sampling data structures, etc.)
     */
    void activate();

    /// Add a child object to the scene (meshes, integrators etc.)
    void addChild(NoriObject *obj);

    /// Return a string summary of the scene (for debugging purposes)
    std::string toString() const;

    EClassType getClassType() const { return EScene; }
private:
    std::vector<Mesh *> m_meshes;
    std::vector<Emitter *> m_emitters;
    Integrator *m_integrator = nullptr;
    Sampler *m_sampler = nullptr;
    Camera *m_camera = nullptr;
    Accel *m_accel = nullptr;
    DiscretePDF m_emitterpdf;
    bool m_isprogressive = false;
};

NORI_NAMESPACE_END
