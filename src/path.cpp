#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/mesh.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN

class PathIntegrator : public Integrator {
public:
	PathIntegrator(const PropertyList &props) { }

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
					break;
				}
			}
			const BSDF* bsdf = its.mesh->getBSDF();
			if (!bsdf) {
				break;
			}
			last_specular = !bsdf->isDiffuse();
			if (bsdf->isDiffuse() && need_shading && Frame::cosTheta(wi) > 0) {
				float emitter_shading_pdf = 0.0f, hemisphere_shading_pdf;
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
					emitter_shading_pdf = surface_pdf * emitter_pdf / enFrame.n.dot(-inc_ray) * inc_norm;
					hemisphere_shading_pdf = bsdf->pdf(brec);
					result += alpha * bsdf->eval(brec) * radiance  / (emitter_shading_pdf + hemisphere_shading_pdf) * its.shFrame.n.dot(inc_ray);
				} while (false);
			}
			if (k <= 2 || sampler->next1D() < 0.95f) {
				BSDFQueryRecord brec = BSDFQueryRecord(wi);
				Color3f backup = alpha;
				alpha *= bsdf->sample(brec, sampler->next2D()) / (k <= 2 ? 1.0f : 0.95f);
				ray_ = Ray3f(its.p, its.shFrame.toWorld(brec.wo));
				if (!scene->rayIntersect(ray_, its))
					break;
				if (its.mesh->isEmitter() && bsdf->isDiffuse() && need_shading && Frame::cosTheta(wi) > 0) {
					float emitter_shading_pdf = 0.0f, hemisphere_shading_pdf = bsdf->pdf(brec);
					float emitter_pdf, surface_pdf;
					emitter_pdf = 1.0f / scene->getEmitters().size();
					surface_pdf = its.mesh->getEmitter()->pdf(its.p);
					emitter_shading_pdf = surface_pdf * emitter_pdf / its.shFrame.n.dot(-ray_.d) * (its.p - ray_.o).squaredNorm();
					Color3f radiance = its.mesh->getEmitter()->getRadiance(its.p, its.shFrame.toLocal(-ray_.d));
					if(radiance.maxCoeff() > 0)
						result += backup * bsdf->eval(brec) * radiance * Frame::cosTheta(brec.wo) / (emitter_shading_pdf + hemisphere_shading_pdf);
				}
			}
			else {
				break;
			}
			k++;
		}
		return result;
	}

	std::string toString() const {
		return "PathIntegrator[]";
	}
};

NORI_REGISTER_CLASS(PathIntegrator, "path");
NORI_NAMESPACE_END