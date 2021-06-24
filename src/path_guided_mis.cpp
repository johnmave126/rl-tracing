#include <tracer/integrator.h>
#include <tracer/scene.h>
#include <tracer/bsdf.h>
#include <tracer/emitter.h>
#include <tracer/mesh.h>
#include <tracer/sampler.h>
#include <tracer/guider.h>

TRACER_NAMESPACE_BEGIN

#define CHECK_VALID(x) (isnan(x) || isinf(x) || x < 0.0f)

class PathGuidedMISIntegrator : public Integrator {
public:
    PathGuidedMISIntegrator(const PropertyList &props) { }

    void addChild(TracerObject *obj) {
        switch (obj->getClassType()) {
        case EGuider:
            if (m_guider)
                throw TracerException("There can only be one guider per integrator!");
            m_guider = static_cast<Guider*>(obj);
            break;
        default:
            throw TracerException("PathGuidedIntegrator::addChild(<%s>) is not supported!",
                classTypeName(obj->getClassType()));
        }
    }

    void activate() {
        if (!m_guider)
            throw TracerException("No guider was specified!");
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
				}
                else {
                    float emitter_pdf, surface_pdf;
                    emitter_pdf = 1.0f / scene->getEmitters().size();
                    surface_pdf = its.mesh->getEmitter()->pdf(its.p);
                    float geom = (its.p - ray_.o).squaredNorm() / abs(Frame::cosTheta(wi));
                    float emitter_shading_pdf = emitter_pdf * surface_pdf * geom;
                    float hemisphere_shading_pdf = m_guider->pdf(last_its.shFrame.toLocal((its.p - last_its.p).normalized()), last_its);
                    bool isresult_nan = CHECK_VALID(result.r());
                    result += alpha * its.mesh->getEmitter()->getRadiance(its.p, wi) * hemisphere_shading_pdf / (emitter_shading_pdf + hemisphere_shading_pdf);
                    if (!isresult_nan && CHECK_VALID(result.r())) {
                        cout << tfm::format("65: alpha: %s\nh_pdf: %f, e_pdf: %f, geom: %f\nits.p: %s, last_its.p: %s, its.n: %s, last_its.n: %s\n",
                            alpha.toString(), hemisphere_shading_pdf, emitter_shading_pdf, geom, its.p.toString(), last_its.p.toString(), its.shFrame.n.toString(), last_its.shFrame.n.toString());
                    }
                }
			}
            if (k > 0) {
                //Update Guider
                m_guider->update(last_its, its, sampler);
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
					if (its.shFrame.n.dot(inc_ray) <= 0 || enFrame.n.dot(-inc_ray) <= 0 || radiance.sum() < Epsilon)
						break;
					float inc_norm = inc_ray.squaredNorm();
                    if (scene->rayIntersect(Ray3f(its.p, inc_ray, Epsilon, 1 - Epsilon)))
                        break;
                    //Intersection e_its;
                    //if (!scene->rayIntersect(Ray3f(its.p, inc_ray), e_its) || (source - e_its.p).norm() > Epsilon)
                    //    break;
                    inc_ray.normalize();
                    Vector3f local_inc_ray = its.shFrame.toLocal(inc_ray);

					BSDFQueryRecord brec = BSDFQueryRecord(wi, local_inc_ray, ESolidAngle);
					emitter_shading_pdf = surface_pdf * emitter_pdf / abs(enFrame.n.dot(inc_ray)) * inc_norm;
                    hemisphere_shading_pdf = m_guider->pdf(local_inc_ray, its);
                    bool isresult_nan = CHECK_VALID(result.r());
					result += alpha * bsdf->eval(brec) * radiance  / (emitter_shading_pdf + hemisphere_shading_pdf) * Frame::cosTheta(local_inc_ray);
                    if (!isresult_nan && CHECK_VALID(result.r())) {
                        cout << tfm::format("112: alpha: %s\nh_pdf: %f, e_pdf: %f, radiance: %s\nits.p: %s, source: %s, its.n: %s, enFrame.n: %s\n",
                            alpha.toString(), hemisphere_shading_pdf, emitter_shading_pdf, radiance.toString(), its.p.toString(), source.toString(), its.shFrame.n.toString(), enFrame.n.toString());
                    }
                    //m_guider->update(its, e_its, sampler);
				} while (false);
			}
			if (k <= 2 || sampler->next1D() < 0.95f) {
                BSDFQueryRecord brec = BSDFQueryRecord(wi);
                if (last_specular) {
                    alpha *= bsdf->sample(brec, sampler->next2D()) / (k <= 2 ? 1.0f : 0.95f);
                }
                else {
                    //Use guider to decide next direction
                    float pdf;
                    brec.wo = m_guider->sample(sampler->next2D(), its, pdf);
                    brec.measure = ESolidAngle;
                    pdf *= k <= 2 ? 1.0f : 0.95f;
                    alpha *= bsdf->eval(brec) * Frame::cosTheta(brec.wo) / pdf;
                }
                last_its = its;
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

    void done() {
        m_guider->done();
    }

	std::string toString() const {
		return "PathGuidedMISIntegrator[]";
	}
protected:
    Guider* m_guider = nullptr;
};

TRACER_REGISTER_CLASS(PathGuidedMISIntegrator, "path_guided_mis");
TRACER_NAMESPACE_END