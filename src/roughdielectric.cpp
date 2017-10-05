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

#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/warp.h>

inline float power5(float x) {
	float y = x * x;
	return y * y * x;
}

NORI_NAMESPACE_BEGIN

/// Ideal dielectric BSDF
class RoughDielectric : public BSDF {
public:
	RoughDielectric(const PropertyList &propList) {
		/* RMS surface roughness */
		m_alpha = propList.getFloat("alpha", 0.1f);

		/* Interior IOR (default: BK7 borosilicate optical glass) */
		m_intIOR = propList.getFloat("intIOR", 1.5046f);

		/* Exterior IOR (default: air) */
		m_extIOR = propList.getFloat("extIOR", 1.000277f);

    }

    Color3f eval(const BSDFQueryRecord &) const {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return Color3f(0.0f);
    }

    float pdf(const BSDFQueryRecord &) const {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return 0.0f;
    }

    Color3f sample(BSDFQueryRecord &bRec, const Point2f &_sample) const {
		Point2f sample = _sample;
		Vector3f wn = Warp::squareToBeckmann(sample, m_alpha);
		bRec.wo = (2.0f * wn.dot(bRec.wi) * wn - bRec.wi).normalized();
		// Note: Once you have implemented the part that computes the scattered
		// direction, the last part of this function should simply return the
		// BRDF value divided by the solid angle density and multiplied by the
		// cosine factor from the reflection equation, i.e.
		return eval(bRec) * Frame::cosTheta(bRec.wo) / pdf(bRec);
    }

    std::string toString() const {
        return tfm::format(
            "RoughDielectric[\n"
            "  alpha = %f,\n"
            "  intIOR = %f,\n"
            "  extIOR = %f\n"
            "]",
            m_alpha, m_intIOR, m_extIOR);
    }
private:
	float m_intIOR, m_extIOR, m_alpha;

	float G1(const Vector3f& wv, const Vector3f& wh) const {
		if (wv.dot(wh) / wv.z() <= -FLT_EPSILON)
			return 0;
		float b = wv.z() / m_alpha / sqrt(1.0f - wv.z() * wv.z());
		if (b < 1.6f) {
			return (3.535f * b + 2.181f * b * b) / (1.0f + 2.276f * b + 2.577f * b * b);
		}
		else {
			return 1.0f;
		}
	}
};

NORI_REGISTER_CLASS(RoughDielectric, "roughdielectric");
NORI_NAMESPACE_END
