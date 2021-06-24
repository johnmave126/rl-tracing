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

#include <tracer/bsdf.h>

TRACER_NAMESPACE_BEGIN

/**
 * \brief Diffuse / Lambertian BRDF model
 */
class Probe : public BSDF {
public:
    Probe(const PropertyList &) { }

    Color3f eval(const BSDFQueryRecord &bRec) const {
        return Color3f(0.0f);
    }

    /// Compute the density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord &bRec) const {
        return 0.0f;
    }

    /// Draw a a sample from the BRDF model
    Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const {
        return Color3f(0.0f);
    }

    bool isDiffuse() const {
        return true;
    }

    bool isProbe() const {
        return true;
    }

    /// Return a human-readable summary
    std::string toString() const {
        return "Probe[]";
    }

    EClassType getClassType() const { return EBSDF; }
private:
};

TRACER_REGISTER_CLASS(Probe, "probe");
TRACER_NAMESPACE_END
