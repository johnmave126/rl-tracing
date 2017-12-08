#include <nori/integrator.h>
#include <nori/scene.h>

NORI_NAMESPACE_BEGIN

class SimpleIntegrator : public Integrator {
public:
	SimpleIntegrator(const PropertyList &props) {
		m_position = props.getPoint("position");
		m_energy = props.getColor("energy");
	}

	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
		/* Find the surface that is visible in the requested direction */
		Intersection its;
		if (!scene->rayIntersect(ray, its))
			return Color3f(0.0f);
		Vector3f direction = m_position - its.p;
		/* Check shadow ray*/
		if (scene->rayIntersect(Ray3f(its.p, direction, Epsilon, 1.0f - Epsilon)))
			return Color3f(0.0f);
		return m_energy * std::max(0.0f, its.shFrame.n.dot(direction) / direction.norm()) / (4 * M_PI * M_PI) / direction.squaredNorm();
	}

	std::string toString() const {
		return tfm::format(
			"SimpleIntegrator[\n"
			"  position = %s,\n"
			"  energy = %s\n"
			"]",
			m_position.toString(),
			m_energy.toString()
		);
	}
private:
	Point3f m_position;
	Color3f m_energy;
};

NORI_REGISTER_CLASS(SimpleIntegrator, "simple");
NORI_NAMESPACE_END