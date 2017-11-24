#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/mesh.h>
#include <nori/sampler.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class PathNaiveIntegrator : public Integrator {
public:
    PathNaiveIntegrator(const PropertyList &props) { }

	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
		/* Find the surface that is visible in the requested direction */
		Intersection its;
		Ray3f ray_ = ray;
		if (!scene->rayIntersect(ray_, its))
			return Color3f(0.0f);
		Color3f alpha = Color3f(1.0f);
        Vector3f last_ray;
		while (true) {
			const Vector3f wi = its.shFrame.toLocal(-ray_.d.normalized());
			if (its.mesh->isEmitter()) {
                //cout << tfm::format("%d %s %s\n", k, alpha, its.mesh->toString());
				return alpha * its.mesh->getEmitter()->getRadiance(its.p, wi);
			}
			const BSDF* bsdf = its.mesh->getBSDF();
			if (!bsdf) {
				break;
			}
			if (bsdf->isDiffuse()) {
                last_ray = Warp::squareToCosineHemisphere(sampler->next2D());
                BSDFQueryRecord brec = BSDFQueryRecord(wi, last_ray, ESolidAngle);
                alpha *= bsdf->eval(brec);
			}
            else if (!bsdf->isDiffuse()) {
                BSDFQueryRecord brec = BSDFQueryRecord(wi);
                alpha *= bsdf->sample(brec, sampler->next2D());
                last_ray = brec.wo;
            }
            else {
                break;
            }
            ray_ = Ray3f(its.p, its.shFrame.toWorld(last_ray));
            if (!scene->rayIntersect(ray_, its))
                break;
		}
        return Color3f(0.0f);
	}

	std::string toString() const {
        return "PathNaiveIntegrator[]";
	}
};

NORI_REGISTER_CLASS(PathNaiveIntegrator, "path_naive");
NORI_NAMESPACE_END