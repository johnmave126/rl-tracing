#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/emitter.h>
#include <nori/sampler.h>
#include <nori/bsdf.h>

NORI_NAMESPACE_BEGIN

class WhittedIntegrator : public Integrator {
public:
    WhittedIntegrator(const PropertyList &props) { }

	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
		/* Find the surface that is visible in the requested direction */
		Intersection its;
		if (!scene->rayIntersect(ray, its))
			return Color3f(0.0f);
		const BSDF* bsdf = its.mesh->getBSDF();
		if (bsdf->isDiffuse()) {
			float emitter_pdf, surface_pdf;
			const Emitter* emitter = scene->sampleEmitter(sampler->next1D(), emitter_pdf);
			if (!emitter)
				return Color3f(0.0f);
			Point3f source;
			Frame enFrame;
			Color3f radiance = emitter->sample(its.p, sampler->next2D(), source, enFrame, surface_pdf);
			Vector3f inc_ray = source - its.p;
			float inc_norm = inc_ray.squaredNorm();
			if (scene->rayIntersect(Ray3f(its.p, inc_ray, Epsilon, 1.0f - Epsilon)))
				return Color3f(0.0f);
			inc_ray.normalize();
			BSDFQueryRecord brec = BSDFQueryRecord(its.shFrame.toLocal(-ray.d.normalized()), its.shFrame.toLocal(inc_ray), ESolidAngle);
			return bsdf->eval(brec) * radiance * (abs(its.shFrame.n.dot(inc_ray) * enFrame.n.dot(inc_ray)) / inc_norm / surface_pdf / emitter_pdf);
		}
		else if(sampler->next1D() < 0.95f) {
			//Go ahead
			BSDFQueryRecord brec = BSDFQueryRecord(its.shFrame.toLocal(-ray.d.normalized()));
			Color3f c = bsdf->sample(brec, sampler->next2D());
			return c * Li(scene, sampler, Ray3f(its.p, its.shFrame.toWorld(brec.wo))) / 0.95f;
		}
		return Color3f(0.0f);
	}

	std::string toString() const {
        return "WhittedIntegrator[]";
	}
};

NORI_REGISTER_CLASS(WhittedIntegrator, "whitted");
NORI_NAMESPACE_END