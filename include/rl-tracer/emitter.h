/*
    This file is part of Tracer, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    [redacted] is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    [redacted] is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once

#include <tracer/object.h>
#include <tracer/frame.h>

TRACER_NAMESPACE_BEGIN

/**
 * \brief Superclass of all emitters
 */
class Emitter : public TracerObject {
public:

    /**
     * \brief Return the type of object (i.e. Mesh/Emitter/etc.) 
     * provided by this instance
     * */
    EClassType getClassType() const { return EEmitter; }

    // Sample a point on a mesh return the radiance, sampled point, normal, and pdf
    // Input: Original point (in case of a non-uniform emitter), A [0, 1]^2 uniformly sample
    virtual Color3f sample(const Point3f& origin, const Point2f& sample, Point3f &p, Frame &nFrame, float &pdf) const = 0;

	// Get radiance of a direction (local coordinate) at a certain point
	virtual Color3f getRadiance(const Point3f& p, const Vector3f& d) const = 0;

	// Get pdf of a sampled point
	virtual float pdf(const Point3f& p) const = 0;
};

TRACER_NAMESPACE_END
