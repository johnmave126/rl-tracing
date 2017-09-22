#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class AOIntegrator : public Integrator {
public:
	AOIntegrator(const PropertyList &props) { }

	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
		/* Find the surface that is visible in the requested direction */
		Intersection its;
		if (!scene->rayIntersect(ray, its))
			return Color3f(0.0f);
		Vector3f w = Warp::squareToCosineHemisphere(sampler->next2D());
		Vector3f shadow = its.shFrame.toWorld(w);
		return Color3f(!scene->rayIntersect(Ray3f(its.p, shadow)));
	}

	std::string toString() const {
		return tfm::format("AmbientOcclusionIntegrator[]");
	}
};

NORI_REGISTER_CLASS(AOIntegrator, "ao");
NORI_NAMESPACE_END