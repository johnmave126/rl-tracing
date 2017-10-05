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

NORI_NAMESPACE_BEGIN

class Microfacet : public BSDF {
public:
    Microfacet(const PropertyList &propList) {
        /* RMS surface roughness */
        m_alpha = propList.getFloat("alpha", 0.1f);

        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);

        /* Albedo of the diffuse base material (a.k.a "kd") */
        m_kd = propList.getColor("kd", Color3f(0.5f));

        /* To ensure energy conservation, we must scale the 
           specular component by 1-kd. 

           While that is not a particularly realistic model of what 
           happens in reality, this will greatly simplify the 
           implementation. Please see the course staff if you're 
           interested in implementing a more realistic version 
           of this BRDF. */
        m_ks = 1 - m_kd.maxCoeff();
    }

    /// Evaluate the BRDF for the given pair of directions
    Color3f eval(const BSDFQueryRecord &bRec) const {
		Vector3f wh = (bRec.wi + bRec.wo).normalized();
		return m_kd * INV_PI + 
			m_ks * Warp::squareToBeckmannPdf(wh, m_alpha) * fresnel(wh.dot(bRec.wi), m_extIOR, m_intIOR)
				* G1(bRec.wi, wh) * G1(bRec.wo, wh) / (4 * Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) * Frame::cosTheta(wh));
    }

    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord &bRec) const {
		if (Frame::cosTheta(bRec.wo) <= 0.0f || Frame::cosTheta(bRec.wi) <= 0.0f)
			return 0.0f;
		Vector3f wh = (bRec.wi + bRec.wo).normalized();
		return m_ks * Warp::squareToBeckmannPdf(wh, m_alpha) / (4.0f * wh.dot(bRec.wo)) + (1 - m_ks) * std::max(Frame::cosTheta(bRec.wo), 0.0f)  * INV_PI;
    }

    /// Sample the BRDF
    Color3f sample(BSDFQueryRecord &bRec, const Point2f &_sample) const {
		Point2f sample = _sample;
		if (Frame::cosTheta(bRec.wi) <= 0)
			return Color3f(0.0f);
		if (sample.x() < m_ks) {
			//Reflect
			sample.x() /= m_ks;
			Vector3f wn = Warp::squareToBeckmann(sample, m_alpha);
			bRec.wo = (2.0f * wn.dot(bRec.wi) * wn - bRec.wi).normalized();
		}
		else {
			//Diffuse
			sample.x() = (sample.x() - m_ks) / (1.0f - m_ks);
			bRec.wo = Warp::squareToCosineHemisphere(sample);
		}
		bRec.measure = ESolidAngle;
		bRec.eta = 1.0f;
		if (Frame::cosTheta(bRec.wo) <= 0)
			return Color3f(0.0f);
        // Note: Once you have implemented the part that computes the scattered
        // direction, the last part of this function should simply return the
        // BRDF value divided by the solid angle density and multiplied by the
        // cosine factor from the reflection equation, i.e.
        return eval(bRec) * Frame::cosTheta(bRec.wo) / pdf(bRec);
    }

    bool isDiffuse() const {
        /* While microfacet BRDFs are not perfectly diffuse, they can be
           handled by sampling techniques for diffuse/non-specular materials,
           hence we return true here */
        return true;
    }

    std::string toString() const {
        return tfm::format(
            "Microfacet[\n"
            "  alpha = %f,\n"
            "  intIOR = %f,\n"
            "  extIOR = %f,\n"
            "  kd = %s,\n"
            "  ks = %f\n"
            "]",
            m_alpha,
            m_intIOR,
            m_extIOR,
            m_kd.toString(),
            m_ks
        );
    }
private:
    float m_alpha;
    float m_intIOR, m_extIOR;
    float m_ks;
    Color3f m_kd;

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

NORI_REGISTER_CLASS(Microfacet, "microfacet");
NORI_NAMESPACE_END
