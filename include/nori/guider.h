#pragma once

#include <nori/object.h>
#include <nori/mesh.h>

NORI_NAMESPACE_BEGIN

/**
 * \brief Ray guider interface
 * 
 * This class provides an abstract interface to guider in Nori and
 * exposes the ability to sample a ray direction. Mostly likely a
 * Reinforced Learning technique will be behind it.
 */
class Guider : public NoriObject {
public:
    /**
     * \brief Initialize the guider using scene information. Pre-train,
	 * loading, etc. should happen here.
     *
     * \param scene
     *    The scene object
     *
     */
    virtual void init(const Scene* scene) { };

	/**
	* \brief Sample a favorable direction at an intersection.
	* Returned vector is in local coordinate.
	*
	* \param sample
	*    A random uniform [0, 1)^2 sample
    *
	* \param its
	*    The current intersection
    *
	* \param pdf
	*    Return the pdf of the point
    *
	* \return the sampled direction
	*
	*/
    virtual Vector3f sample(const Point2f& sample, const Intersection& its, float& pdf) = 0;

    /**
    * \brief Update the guider according to the next intersection
    *
    * \param origin
    *    The original state point
    *
    * \param di
    *    The action taken in local coordinate
    *
    * \param its
    *    The next intersection
    *
    * \param sampler
    *    Provide a random number generator for the method
    *
    */
    virtual void update(const Point3f& origin, const Vector3f& di, const Intersection& its, Sampler* sampler) = 0;

    /**
    * \brief Return the pdf of a direction
    *
    * \param di
    *    The original state point
    *
    * \param origin
    *    The action taken in local coordinate
    *
    * \return the pdf of the direction
    */
    virtual float pdf(const Vector3f& di, const Point3f& origin) = 0;

    EClassType getClassType() const { return EGuider; }
};

NORI_NAMESPACE_END
