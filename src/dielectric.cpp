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
#include <tracer/frame.h>

inline float power5(float x) {
	float y = x * x;
	return y * y * x;
}

TRACER_NAMESPACE_BEGIN

/// Ideal dielectric BSDF
class Dielectric : public BSDF {
public:
    Dielectric(const PropertyList &propList) {
        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);

    }

    Color3f eval(const BSDFQueryRecord &) const {
        /* Discrete BRDFs always evaluate to zero in [redacted] */
        return Color3f(0.0f);
    }

    float pdf(const BSDFQueryRecord &) const {
        /* Discrete BRDFs always evaluate to zero in [redacted] */
        return 0.0f;
    }

    Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const {
		float cosd = Frame::cosTheta(bRec.wi);
		float eta = cosd <= 0.0f ? m_intIOR / m_extIOR : m_extIOR / m_intIOR;
		float sint2 = (1 - cosd * cosd) * eta * eta;
		float cost = cosd <= 0.0f ? sqrt(1.0f - sint2) : -sqrt(1.0f - sint2);
		float Fr = fresnel(cosd, m_extIOR, m_intIOR);
		bRec.measure = EDiscrete;
		if (sample.x() < Fr) {
			// Reflection in local coordinates
			bRec.wo = Vector3f(
				-bRec.wi.x(),
				-bRec.wi.y(),
				bRec.wi.z()
			);

			/* Relative index of refraction: no change */
			bRec.eta = 1.0f;
			return Color3f(1.0f);
		}
		else {
			bRec.eta = eta;
			bRec.wo = (-bRec.wi * bRec.eta + Vector3f(0.0f, 0.0f, cosd * bRec.eta + cost)).normalized();
			return Color3f(1.0f);
		}
    }

    std::string toString() const {
        return tfm::format(
            "Dielectric[\n"
            "  intIOR = %f,\n"
            "  extIOR = %f\n"
            "]",
            m_intIOR, m_extIOR);
    }
private:
    float m_intIOR, m_extIOR;
};

TRACER_REGISTER_CLASS(Dielectric, "dielectric");
TRACER_NAMESPACE_END
