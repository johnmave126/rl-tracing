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
		Intersection its, last_its;
		Ray3f ray_ = ray;
		if (!scene->rayIntersect(ray_, its))
			return Color3f(0.0f);
		Color3f alpha = Color3f(1.0f);
		int k = 0;
		while (true) {
			const Vector3f wi = its.shFrame.toLocal(-ray_.d.normalized());
            if (k > 0) {
                m_guider->update(last_its, its, sampler);
            }
			if (its.mesh->isEmitter()) {
				return alpha * its.mesh->getEmitter()->getRadiance(its.p, wi);
			}
			const BSDF* bsdf = its.mesh->getBSDF();
			if (!bsdf) {
				break;
			}
            BSDFQueryRecord brec = BSDFQueryRecord(wi);
			if (bsdf->isDiffuse()) {
                float pdf;
                brec.wo = m_guider->sample(sampler->next2D(), its, pdf);
                brec.measure = ESolidAngle;
                alpha *= bsdf->eval(brec) * Frame::cosTheta(brec.wo) / pdf;
			}
            else {
                alpha *= bsdf->sample(brec, sampler->next2D());
            }
            ray_ = Ray3f(its.p, its.shFrame.toWorld(brec.wo));
            last_its = its;
            if (!scene->rayIntersect(ray_, its))
                break;
			k++;
		}
        return Color3f(0.0f);
	}

    void done() {
        m_guider->done();
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