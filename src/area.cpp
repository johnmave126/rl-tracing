#include <tracer/emitter.h>
#include <tracer/mesh.h>

NORI_NAMESPACE_BEGIN

class AreaLight : public Emitter {
public:
    AreaLight(const PropertyList &props) {
        m_radiance = props.getColor("radiance");
    }

    Color3f sample(const Point3f& origin, const Point2f& sample, Point3f &p, Frame &nFrame, float &pdf) const {
        if (!m_mesh)
            throw NoriException("AreaLight must be a children of some mesh");
        m_mesh->samplePosition(sample, p, nFrame, pdf);
		return getRadiance(p, nFrame.toLocal(origin - p).normalized());
    }

	Color3f getRadiance(const Point3f& p, const Vector3f& d) const {
		if (Frame::cosTheta(d) < 0.0f)
			return Color3f(0.0f);
		return m_radiance;
	}

	float pdf(const Point3f& p) const {
		return 1.0f / m_mesh->getSurfaceArea();
	}

    std::string toString() const {
        return tfm::format(
			"AreaLight[\n"
			"  radiance = %s\n"
			"]", m_radiance.toString());
    }

    void setParent(NoriObject *parent) {
        m_mesh = dynamic_cast<const Mesh*>(parent);
    }

protected:
    const Mesh* m_mesh;
    Color3f m_radiance;
};

NORI_REGISTER_CLASS(AreaLight, "area");
NORI_NAMESPACE_END
