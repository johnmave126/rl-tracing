#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/mesh.h>
#include <nori/sampler.h>
#include <nori/guider.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class PathGuidedIntegrator : public Integrator {
public:
    PathGuidedIntegrator(const PropertyList &props) { }

    void addChild(NoriObject *obj) {
        switch (obj->getClassType()) {
            case EGuider:
                if (m_guider)
                    throw NoriException("There can only be one guider per integrator!");
                m_guider = static_cast<Guider*>(obj);
                break;
            default:
                throw NoriException("PathGuidedIntegrator::addChild(<%s>) is not supported!",
                    classTypeName(obj->getClassType()));
        }
    }

    void activate() {
        if (!m_guider)
            throw NoriException("No guider was specified!");
    }

    void preprocess(const Scene *scene) {
        m_guider->init(scene);
    }

	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
		/* Find the surface that is visible in the requested direction */
		Intersection its;
		Ray3f ray_ = ray;
		if (!scene->rayIntersect(ray_, its))
			return Color3f(0.0f);
		Color3f alpha = Color3f(1.0f);
		int k = 0;
        Vector3f last_ray;
		while (true) {
			const Vector3f wi = its.shFrame.toLocal(-ray_.d.normalized());
            if (k > 0) {
                m_guider->update(ray_.o, last_ray, its, sampler);
            }
			if (its.mesh->isEmitter()) {
                //cout << tfm::format("%d %s %s\n", k, alpha, its.mesh->toString());
				return alpha * its.mesh->getEmitter()->getRadiance(its.p, wi);
			}
			const BSDF* bsdf = its.mesh->getBSDF();
			if (!bsdf) {
				break;
			}
			if (bsdf->isDiffuse()) {
                float pdf = INV_TWOPI;
                last_ray = /*Warp::squareToUniformHemisphere(sampler->next2D());*/m_guider->sample(sampler->next2D(), its, pdf);
                BSDFQueryRecord brec = BSDFQueryRecord(wi, last_ray, ESolidAngle);
                alpha *= bsdf->eval(brec) * Frame::cosTheta(last_ray) / pdf;
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
			k++;
		}
        return Color3f(0.0f);
	}

	std::string toString() const {
        return tfm::format(
            "PathGuidedIntegrator[\n"
            "  guider = %s\n"
            "]",
            indent(m_guider->toString())
        );
	}
protected:
    Guider* m_guider = nullptr;
};

NORI_REGISTER_CLASS(PathGuidedIntegrator, "path_guided");
NORI_NAMESPACE_END