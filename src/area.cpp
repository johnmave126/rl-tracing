#include <emitter.h>

NORI_NAMESPACE_BEGIN

class AreaLight : public Emitter {
public:
    AreaLight(const PropertyList &props) {
        m_radiance = props.getColor('radiance');
    }

    std::string toString() const {
        return tfm::format("AreaLight[radiance = %s]", m_radiance.toString());
    }

protected:
    Color3f m_radiance;
};

NORI_REGISTER_CLASS(AreaLight, "AreaLight");
NORI_NAMESPACE_END
