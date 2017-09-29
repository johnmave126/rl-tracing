#include <nori/emitter.h>
#include <nori/mesh.h>

NORI_NAMESPACE_BEGIN

class AreaLight : public Emitter {
public:
    AreaLight(const PropertyList &props) {
        m_radiance = props.getColor("radiance");
    }

    Color3f sample(const Point3f& origin, const Point2f& sample, Point3f &p, Normal3f &n, float &pdf) const {
        if (!m_mesh)
            throw NoriException("AreaLight must be a children of some mesh");
        m_mesh->samplePosition(sample, p, n, pdf);
        if ((p - origin).dot(n) > 0.0f) {
            return m_radiance;
        }
        return Color3f(0.0f, 0.0f, 0.0f);
    }

    std::string toString() const {
        return tfm::format("AreaLight[radiance = %s]", m_radiance.toString());
    }

    void setParent(NoriObject *parent) {
        m_mesh = dynamic_cast<const Mesh*>(parent);
    }

protected:
    const Mesh* m_mesh;
    Color3f m_radiance;
};

NORI_REGISTER_CLASS(AreaLight, "AreaLight");
NORI_NAMESPACE_END
