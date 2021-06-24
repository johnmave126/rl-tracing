#include <tracer/integrator.h>
#include <tracer/scene.h>
#include <tracer/bsdf.h>
#include <tracer/emitter.h>
#include <tracer/mesh.h>
#include <tracer/sampler.h>

NORI_NAMESPACE_BEGIN

class PathSimpleIntegrator : public Integrator {
public:
	PathSimpleIntegrator(const PropertyList &props) { }

	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
		/* Find the surface that is visible in the requested direction */
		Intersection its;
		Ray3f ray_ = ray;
		if (!scene->rayIntersect(ray_, its))
			return Color3f(0.0f);
		Color3f result = Color3f(0.0f);
		Color3f alpha = Color3f(1.0f);
		int k = 0;
		bool last_specular = false;
		while (true) {
			bool need_shading = true;
			const Vector3f wi = its.shFrame.toLocal(-ray_.d.normalized());
			if (its.mesh->isEmitter()) {
				if (last_specular || k == 0) {
					//Last hop specular or primary ray
					result += alpha * its.mesh->getEmitter()->getRadiance(its.p, wi);
					need_shading = false;
				}
			}
			const BSDF* bsdf = its.mesh->getBSDF();
			if (!bsdf) {
				break;
			}
			last_specular = !bsdf->isDiffuse();
			if (bsdf->isDiffuse() && need_shading && Frame::cosTheta(wi) > 0) {
				//Shade it
				float emitter_pdf, surface_pdf;
				const Emitter* emitter = scene->sampleEmitter(sampler->next1D(), emitter_pdf);
				do {
					if (!emitter)
						break;
					Point3f source;
					Frame enFrame;
					Color3f radiance = emitter->sample(its.p, sampler->next2D(), source, enFrame, surface_pdf);
					Vector3f inc_ray = source - its.p;
					if (its.shFrame.n.dot(inc_ray) <= 0)
						break;
					float inc_norm = inc_ray.squaredNorm();
					if (scene->rayIntersect(Ray3f(its.p, inc_ray, Epsilon, 1.0f - Epsilon)))
						break;
					inc_ray.normalize();
					BSDFQueryRecord brec = BSDFQueryRecord(wi, its.shFrame.toLocal(inc_ray), ESolidAngle);
					result += alpha * bsdf->eval(brec) * radiance * (its.shFrame.n.dot(inc_ray) * enFrame.n.dot(-inc_ray) / inc_norm / surface_pdf / emitter_pdf);
				} while (false);
			}
			if (k <= 2 || sampler->next1D() < 0.95f) {
				BSDFQueryRecord brec = BSDFQueryRecord(wi);
				alpha *= bsdf->sample(brec, sampler->next2D())  / (k <= 2 ? 1.0f : 0.95f);
				ray_ = Ray3f(its.p, its.shFrame.toWorld(brec.wo));
				if (!scene->rayIntersect(ray_, its))
					break;
			}
			else {
				break;
			}
			k++;
		}
		return result;
	}

	std::string toString() const {
		return "PathSimpleIntegrator[]";
	}
};

NORI_REGISTER_CLASS(PathSimpleIntegrator, "path_simple");
NORI_NAMESPACE_END