/*
    This file is part of Tracer, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    [redacted] is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    [redacted] is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <tracer/scene.h>
#include <tracer/bitmap.h>
#include <tracer/integrator.h>
#include <tracer/sampler.h>
#include <tracer/camera.h>
#include <tracer/emitter.h>
#include <tracer/timer.h>

TRACER_NAMESPACE_BEGIN

Scene::Scene(const PropertyList &props) {
    m_accel = new Accel();
    m_isprogressive = props.getBoolean("progressive", false);
}

Scene::~Scene() {
    delete m_accel;
    delete m_sampler;
    delete m_camera;
    delete m_integrator;
}

void Scene::activate() {
    m_accel->build();

    if (!m_integrator)
        throw TracerException("No integrator was specified!");
    if (!m_camera)
        throw TracerException("No camera was specified!");
    
    if (!m_sampler) {
        /* Create a default (independent) sampler */
        m_sampler = static_cast<Sampler*>(
            TracerObjectFactory::createInstance("independent", PropertyList()));
    }
    for (auto it = m_meshes.begin(); it != m_meshes.end(); ++it) {
        if ((*it)->isEmitter()) {
            m_emitters.push_back((*it)->getEmitter());
            m_emitterpdf.append(1.0f);
        }
    }
    m_emitterpdf.normalize();

    cout << endl;
    cout << "Configuration: " << toString() << endl;
    cout << endl;
}

const Emitter* Scene::sampleEmitter(float &sample, float &pdf) const {
    if (m_emitters.size() == 0) {
        return nullptr;
    }
    size_t idx = m_emitterpdf.sampleReuse(sample);
    pdf = 1.0f / m_emitters.size();
    return m_emitters[idx];
}

const Emitter* Scene::sampleEmitter(const float& sample, float &pdf) const {
	float s = sample;
	return sampleEmitter(s, pdf);
}

void Scene::addChild(TracerObject *obj) {
    switch (obj->getClassType()) {
        case EMesh: {
                Mesh *mesh = static_cast<Mesh *>(obj);
                m_accel->addMesh(mesh);
                m_meshes.push_back(mesh);
            }
            break;
        
        case EEmitter: {
                //Emitter *emitter = static_cast<Emitter *>(obj);
                /* TBD */
                throw TracerException("Scene::addChild(): You need to implement this for emitters");
            }
            break;

        case ESampler:
            if (m_sampler)
                throw TracerException("There can only be one sampler per scene!");
            m_sampler = static_cast<Sampler *>(obj);
            break;

        case ECamera:
            if (m_camera)
                throw TracerException("There can only be one camera per scene!");
            m_camera = static_cast<Camera *>(obj);
            break;
        
        case EIntegrator:
            if (m_integrator)
                throw TracerException("There can only be one integrator per scene!");
            m_integrator = static_cast<Integrator *>(obj);
            break;

        default:
            throw TracerException("Scene::addChild(<%s>) is not supported!",
                classTypeName(obj->getClassType()));
    }
}

std::string Scene::toString() const {
    std::string meshes;
    for (size_t i=0; i<m_meshes.size(); ++i) {
        meshes += std::string("  ") + indent(m_meshes[i]->toString(), 2);
        if (i + 1 < m_meshes.size())
            meshes += ",";
        meshes += "\n";
    }

    return tfm::format(
        "Scene[\n"
        "  integrator = %s,\n"
        "  sampler = %s\n"
        "  camera = %s,\n"
        "  meshes = {\n"
        "  %s  }\n"
        "]",
        indent(m_integrator->toString()),
        indent(m_sampler->toString()),
        indent(m_camera->toString()),
        indent(meshes, 2)
    );
}

TRACER_REGISTER_CLASS(Scene, "scene");
TRACER_NAMESPACE_END
