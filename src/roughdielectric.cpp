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

    Color3f eval(const BSDFQueryRecord & bRec) const {
		if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) > 0) {
			//Reflect
			Vector3f wh = (bRec.wi + bRec.wo).normalized();
			return Warp::squareToBeckmannPdf(wh, m_alpha) * fresnel(wh.dot(bRec.wi), m_extIOR, m_intIOR)
				* G1(bRec.wi, wh) * G1(bRec.wo, wh) / abs(4 * Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) * Frame::cosTheta(wh));
		}
		else {
			//Transmission
			Vector3f wh;
			float ii, io;
			if (Frame::cosTheta(bRec.wi) > 0) //Ext->Int
				ii = m_extIOR, io = m_intIOR;
			else
				io = m_extIOR, ii = m_intIOR;
			wh = (-(ii * bRec.wi + io * bRec.wo)).normalized();
			float iioht = ii * bRec.wi.dot(wh) + io * bRec.wo.dot(wh);
			return Warp::squareToBeckmannPdf(wh, m_alpha) * (1.0f - fresnel(wh.dot(bRec.wi), m_extIOR, m_intIOR))
				* G1(bRec.wi, wh) * G1(bRec.wo, wh) * io * io / iioht / iioht * abs(bRec.wi.dot(wh) * bRec.wo.dot(wh) / Frame::cosTheta(bRec.wi) / Frame::cosTheta(bRec.wo));
		}
    }

    float pdf(const BSDFQueryRecord &bRec) const {
		if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) > 0) {
			//Reflect
			Vector3f wh = (bRec.wi + bRec.wo).normalized();
			if (Frame::cosTheta(wh) < 0)
				wh = -wh;
			return Warp::squareToBeckmannPdf(wh, m_alpha) * fresnel(wh.dot(bRec.wi), m_extIOR, m_intIOR) / abs(4.0f * wh.dot(bRec.wo));
		}
		else {
			//Transmission
			Vector3f wh;
			float ii, io;
			if (Frame::cosTheta(bRec.wi) > 0) //Ext->Int
				ii = m_extIOR, io = m_intIOR;
			else
				io = m_extIOR, ii = m_intIOR;
			wh = (-(ii * bRec.wi + io * bRec.wo)).normalized();
			if (Frame::cosTheta(wh) < 0)
				wh = -wh;
			float iioht = ii * bRec.wi.dot(wh) + io * bRec.wo.dot(wh);
			return Warp::squareToBeckmannPdf(wh, m_alpha) * (1 - fresnel(wh.dot(bRec.wi), m_extIOR, m_intIOR)) / abs(4.0f * wh.dot(bRec.wo));
		}
    }

    Color3f sample(BSDFQueryRecord &bRec, const Point2f &_sample) const {
		// We need to reuse the sample
		// PCG32 provides 23 random bits for each float, so we have 46 random bits in total
		// We will provide Beckmann 18+14 random bits and 14 random bits for Fresnel
		float x = ldexp(_sample.x(), 9), y = ldexp(_sample.y(), 5);
		float partialx, partialy;
		Point2f sample = Point2f(modf(x, &partialx), modf(y, &partialy));
		float sample1d = ldexp(ldexp(partialx, 5) + partialy, -14);

		Vector3f wn = Warp::squareToBeckmann(sample, m_alpha);
		float cosd = wn.dot(bRec.wi);
		float Fr = fresnel(cosd, m_extIOR, m_intIOR);
		if (sample1d < Fr) {
			// Reflect
			bRec.wo = (2.0f * wn.dot(bRec.wi) * wn - bRec.wi).normalized();
			bRec.eta = 1.0f;
			if (Frame::cosTheta(bRec.wo) * Frame::cosTheta(bRec.wi) < 0)
				return Color3f(0.0f);
		}
		else {
			// Transmission
			bRec.eta = cosd <= 0.0f ? m_intIOR / m_extIOR : m_extIOR / m_intIOR;
			float sint2 = (1 - cosd * cosd) * bRec.eta * bRec.eta;
			float cost = cosd <= 0.0f ? sqrt(1.0f - sint2) : -sqrt(1.0f - sint2);
			bRec.wo = (-bRec.wi * bRec.eta + (cosd * bRec.eta + cost) * wn).normalized();
			if (Frame::cosTheta(bRec.wo) * Frame::cosTheta(bRec.wi) > 0)
				return Color3f(0.0f);
		}
		bRec.measure = ESolidAngle;
		// Note: Once you have implemented the part that computes the scattered
		// direction, the last part of this function should simply return the
		// BRDF value divided by the solid angle density and multiplied by the
		// cosine factor from the reflection equation, i.e.
		return Color3f(1.0f);
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
		if (wv.dot(wh) / wv.z() <= -Epsilon)
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
