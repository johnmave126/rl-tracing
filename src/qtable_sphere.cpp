#include <functional>
#include <nori/guider.h>
#include <nori/scene.h>
#include <nori/sampler.h>
#include <nori/emitter.h>
#include <nori/integrator.h>
#include <nori/bsdf.h>
#include <nori/warp.h>
#include <tbb/concurrent_vector.h>
#include <tbb/concurrent_hash_map.h>
#include <algorithm>
#include <fstream>
#include <iostream>

#undef NDEBUG
#include <assert.h>

NORI_NAMESPACE_BEGIN

class QTableSphereGuider : public Guider {
protected:
    struct Wrapper {
        float* map;
        int* visit;

        ~Wrapper() {
            if (map)
                delete[] map;
            if (visit)
                delete[] visit;
        }

        void init(int width, int height) {
            map = new float[width * height];
            visit = new int[width * height];
            memset(visit, 0, width * height * sizeof(int));
            for (int i = 0; i < width * height; i++) {
                map[i] = 1.0f;
            }
        }
    };
    typedef tbb::concurrent_hash_map<int, Wrapper> WrapperMap;
public:
    QTableSphereGuider(const PropertyList &props): m_storage(10000) {
        m_sceneResolution = props.getInteger("sceneResolution", 50);
        m_angleResolution = props.getInteger("angleResolution", 8);
        try {
            m_alpha = props.getFloat("alpha");
            m_useVisit = false;
        }
        catch (NoriException e) {
            //alpha doesn't exist
        }
        //Store m_angleResolution x m_angleResolution array for each element in (2 x m_angleResolution) x m_angleResolution array
        m_hemishpereMap = new int[m_angleResolution * m_angleResolution * m_angleResolution * m_angleResolution * 2];
        m_importFilename = props.getString("import", "");
        m_exportFilename = props.getString("export", "");
    }

    ~QTableSphereGuider() {
        delete[] m_hemishpereMap;
    }

    /* Integrator need to call this in preprocess() */
    void init(const Scene *scene) {
        m_sceneBox = scene->getBoundingBox();
        Point3f orig_max = m_sceneBox.max;
        m_sceneBox.expandBy(orig_max + Vector3f(Epsilon));
        m_sceneBlockSize = (m_sceneBox.max - m_sceneBox.min) / m_sceneResolution;

        int width = (m_angleResolution << 1), height = m_angleResolution;

        //Pre-compute hemisphere map
        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                Vector3f center_di = Warp::squareToUniformSphere(Point2f((i + 0.5f) / width, (j + 0.5f) / height));
                Frame center_frame = Frame(center_di);
                for (int k = 0; k < m_angleResolution; k++) {
                    for (int l = 0; l < m_angleResolution; l++) {
                        Vector3f local_di = Warp::squareToUniformHemisphere(Point2f((k + 0.5f) / m_angleResolution, (l + 0.5f) / m_angleResolution));
                        int mapped_idx = locateDirection(center_frame.toWorld(local_di));
                        getHemisphereMap(i, j, k, l) = mapped_idx;
                    }
                }
            }
        }
        if (m_importFilename.length() > 0) {
            // Import from file
            std::ifstream file(m_importFilename, std::ios::in | std::ios::binary);
            if (!file.is_open()) {
                throw NoriException("Cannot open file %s", m_importFilename.c_str());
            }
	        cout << "Importing QTable ...";
	        cout.flush();
            float bounds[6];
            file.read((char*)bounds, 6 * sizeof(float));
            m_sceneBox = BoundingBox3f(Point3f(bounds[0], bounds[1], bounds[2]), Point3f(bounds[3], bounds[4], bounds[5]));
            m_sceneBlockSize = (m_sceneBox.max - m_sceneBox.min) / m_sceneResolution;
            while (!file.eof()) {
                int block_idx;
                file.read((char*)&block_idx, sizeof(int));
                WrapperMap::accessor access;
                m_storage.insert(access, block_idx);
                access->second.init(width, height);
                file.read((char*)access->second.map, width * height * sizeof(float));
                file.read((char*)access->second.visit, width * height * sizeof(int));
                access.release();
            }
            file.close();
	        cout << " done." << endl;
        }

        //cout << "Testing locateDirection" << endl;
        //for (int i = 0; i < width; i++) {
        //    for (int j = 0; j < height; j++) {
        //        int x, y;
        //        Point2f s = Point2f((i + 0.5f) / width, (j + 0.5f) / height);
        //        locateDirection(Warp::squareToUniformSphere(s), x, y);
        //        if (i != x || j != y) {
        //            cout << "Oops..." << i << ' ' << x << ' ' << j << ' ' << y << endl;
        //        }
        //    }
        //}
    }
    Vector3f sample(const Point2f& sample, const Intersection& its, float& pdf) {
        int nx, ny;
        int block_idx = locateBlock(its.p);
        locateDirection(its.shFrame.n, nx, ny);
        assert(nx < 2 * m_angleResolution);
        assert(ny < m_angleResolution);
        Point2f result;
        WrapperMap::const_accessor const_access;
        WrapperMap::accessor access;
        auto do_sample = [&](auto& accessor) -> Point2f {
            float *weights = new float[m_angleResolution + 1];
            float total_weight;
            int x, y;
            float px, py;
            float t;
            weights[0] = 0.0f;
            for (int i = 1; i <= m_angleResolution; i++) {
                weights[i] = weights[i - 1];
                for (int j = 0; j < m_angleResolution; j++) {
                    int mapped_idx = this->getHemisphereMap(nx, ny, i - 1, j);
                    assert(mapped_idx < 2 * m_angleResolution * m_angleResolution);
                    weights[i] += accessor->second.map[mapped_idx];
                }
            }
            total_weight = weights[m_angleResolution];
            t = sample.x() * weights[m_angleResolution];
            x = std::upper_bound(weights, weights + m_angleResolution + 1, t) - weights - 1;
            assert(x < m_angleResolution);
            px = x + (t - weights[x]) / (weights[x + 1] - weights[x]);
            for (int i = 1; i <= m_angleResolution; i++) {
                weights[i] = weights[i - 1];
                int mapped_idx = this->getHemisphereMap(nx, ny, x, i - 1);
                assert(mapped_idx < 2 * m_angleResolution * m_angleResolution);
                weights[i] += accessor->second.map[mapped_idx];
            }
            t = sample.y() * weights[m_angleResolution];
            y = std::upper_bound(weights, weights + m_angleResolution + 1, t) - weights - 1;
            py = y + (t - weights[y]) / (weights[y + 1] - weights[y]);
            assert(y < m_angleResolution);
            int idx = this->getHemisphereMap(nx, ny, x, y);
            assert(idx < 2 * m_angleResolution * m_angleResolution);
            pdf = accessor->second.map[idx] / total_weight * m_angleResolution * m_angleResolution * INV_TWOPI;
            delete[] weights;
            return Point2f(px / m_angleResolution, py / m_angleResolution);
        };
        if (m_storage.find(const_access, block_idx)) {
            result = do_sample(const_access);
        }
        else {
            if (m_storage.insert(access, block_idx)) {
                access->second.init(2 * m_angleResolution, m_angleResolution);
            }
            result = do_sample(access);
        }
        return Warp::squareToUniformHemisphere(result);
    }

    void update(const Intersection& origin, const Intersection& dest, Sampler* sampler) {
        const Vector3f& ray = (dest.p - origin.p).normalized(),
            dest_wi = dest.shFrame.toLocal(-ray),
            orig_wi = origin.shFrame.toLocal(ray);
        int nx, ny, ox, oy;
        locateDirection(dest.shFrame.n, nx, ny);
        assert(nx < 2 * m_angleResolution);
        assert(ny < m_angleResolution);
        int block_orig_idx = locateBlock(origin.p);
        locateDirectionHemisphere(orig_wi, ox, oy);
        assert(ox < m_angleResolution);
        assert(oy < m_angleResolution);
        int block_dest_idx = locateBlock(dest.p),
            angle_orig_idx = getHemisphereMap(nx, ny, ox, oy);

        assert(angle_orig_idx < 2 * m_angleResolution * m_angleResolution);

        WrapperMap::const_accessor const_access_dest;
        WrapperMap::accessor access_dest, access_orig;

        float integral_term = 0.0f;
        const BSDF *bsdf = dest.mesh->getBSDF();
        BSDFQueryRecord brec = BSDFQueryRecord(dest_wi);
        auto integral_job = [&](auto& accessor) {
            if (bsdf->isDiffuse()) {
                brec.measure = ESolidAngle;
                for (int i = 0; i < m_angleResolution; i++) {
                    for (int j = 0; j < m_angleResolution; j++) {
                        Point2f sample = (sampler->next2D() + Point2f(i, j)) / m_angleResolution;
                        brec.wo = Warp::squareToUniformHemisphere(sample);
                        float eval = bsdf->eval(brec).maxCoeff();
                        float normal_q = accessor->second.map[getHemisphereMap(nx, ny, i, j)];
                        float term = normal_q * Frame::cosTheta(brec.wo) * eval;
                        integral_term += term;
                    }
                }
            }
            else {
                //We need to sample the incident ray because brdf is always 0
                for (int i = 0; i < m_angleResolution * m_angleResolution; i++) {
                    bsdf->sample(brec, sampler->next2D());
                    int idx = locateDirection(dest.shFrame.toWorld(brec.wo));
                    integral_term += accessor->second.map[idx];
                }
            }
            accessor.release();
        };
        if (m_storage.find(const_access_dest, block_dest_idx)) {
            integral_job(const_access_dest);
        }
        else {
            if (m_storage.insert(access_dest, block_dest_idx)) {
                access_dest->second.init(2 * m_angleResolution, m_angleResolution);
            }
            integral_job(access_dest);
        }
        integral_term *= 2.0f * M_PI / m_angleResolution / m_angleResolution;
        if (dest.mesh->isEmitter()) {
            integral_term += dest.mesh->getEmitter()->getRadiance(dest.p, dest_wi).sum();
        }


        if (m_storage.insert(access_orig, block_orig_idx)) {
            access_orig->second.init(2 * m_angleResolution, m_angleResolution);
        }
        float alpha = m_useVisit ? 1.0f / (1 + access_orig->second.visit[angle_orig_idx]) : m_alpha;
        float oldval = access_orig->second.map[angle_orig_idx];
        float newval = (1.0f - alpha) * oldval + alpha * integral_term;
        if (newval < UPDATE_THREASHOLD)
            newval = UPDATE_THREASHOLD;
        access_orig->second.map[angle_orig_idx] = newval;
        access_orig->second.visit[angle_orig_idx]++;
    }

    float pdf(const Vector3f& di, const Intersection& origin) {
        int nx, ny, ox, oy;
        locateDirection(origin.shFrame.n, nx, ny);
        locateDirectionHemisphere(di, ox, oy);
        assert(nx < 2 * m_angleResolution);
        assert(ny < m_angleResolution);
        assert(ox < m_angleResolution);
        assert(oy < m_angleResolution);
        int block_idx = locateBlock(origin.p), angle_idx = getHemisphereMap(nx, ny, ox, oy);
        assert(angle_idx < 2 * m_angleResolution * m_angleResolution);

        WrapperMap::const_accessor const_access;
        WrapperMap::accessor access;

        auto calc_pdf = [&](auto& accessor) -> float {
            float total_weight = 0.0f;
            for (int i = 0; i < m_angleResolution; i++) {
                for (int j = 0; j < m_angleResolution; j++) {
                    int mapped_idx = this->getHemisphereMap(nx, ny, i, j);
                    total_weight += accessor->second.map[mapped_idx];
                }
            }
            return accessor->second.map[angle_idx] / total_weight * m_angleResolution * m_angleResolution * INV_TWOPI;
        };
        if (m_storage.find(const_access, block_idx)) {
            return calc_pdf(const_access);
        }
        else {
            if (m_storage.insert(access, block_idx)) {
                access->second.init(2 * m_angleResolution, m_angleResolution);
            }
            return calc_pdf(access);
        }
        return 0.0f;
    }

    void done() {
        if (m_exportFilename.length() > 0) {
            std::ofstream file(m_exportFilename, std::ios::binary | std::ios::out);
            if (!file.is_open()) {
                throw NoriException("Cannot open file %s to write", m_importFilename.c_str());
            }
            // Export to the file
            std::cout << "Exporting QTableSphere to " << m_exportFilename << " ... ";
            std::cout.flush();
            float bounds[6] = { m_sceneBox.min.x(), m_sceneBox.min.y(), m_sceneBox.min.z(), m_sceneBox.max.x(), m_sceneBox.max.y(), m_sceneBox.max.z() };
            file.write((char*)bounds, 6 * sizeof(float));
            for (WrapperMap::const_iterator it = m_storage.begin(); it != m_storage.end(); ++it) {
                file.write((char*)&it->first, sizeof(int));
                file.write((char*)it->second.map, 2 * m_angleResolution * m_angleResolution * sizeof(float));
                file.write((char*)it->second.visit, 2 * m_angleResolution * m_angleResolution * sizeof(int));
            }
            file.close();
            std::cout << "done." << std::endl;
        }
    }

	std::string toString() const {
        return tfm::format(
            "QTableSphereGuider[\n"
            "  alpha = %s,\n"
            "  sceneResolution = %d,\n"
            "  angleResolution = %d\n"
            "]",
            m_useVisit? "1/(1 + visit)" : tfm::format("%f", m_alpha), m_sceneResolution, m_angleResolution);
	}

protected:

    int locateBlock(const Point3f& pos) const {
        Vector3f offset = pos - m_sceneBox.min;
        int x = offset.x() / m_sceneBlockSize.x(),
            y = offset.y() / m_sceneBlockSize.y(),
            z = offset.z() / m_sceneBlockSize.z();
        return (x * m_sceneResolution + y) * m_sceneResolution + z;
    }

    inline int locateDirection(const Vector3f& di, int &ix, int& iy) const {
        float x = std::min(1 - 1e-6f, di.z());
        float y = 0.0f;
        if (di.z() > -1 + 1e-6f && di.z() < 1 - 1e-6f) {
            y = sphericalCoordinates(di).y();
        }
        x = (x + 1.0f) / 2.0f;
        y /= 2.0f * M_PI;
        ix = x * 2 * m_angleResolution;
        iy = y * m_angleResolution;
        return ix * m_angleResolution + iy;
    }

    int locateDirection(const Vector3f& di) const {
        int dummyx, dummyy;
        return locateDirection(di, dummyx, dummyy);
    }

    inline int locateDirectionHemisphere(const Vector3f& di, int &ix, int& iy) const {
        float x = std::min(1 - 1e-6f, di.z());
        float y = 0.0f;
        if (di.z() < 1 - 1e-6f) {
            y = sphericalCoordinates(di).y();
        }
        y /= 2 * M_PI;
        ix = x * m_angleResolution;
        iy = y * m_angleResolution;
        return ix * m_angleResolution + iy;
    }

    int locateDirectionHemisphere(const Vector3f& di) const {
        int dummyx, dummyy;
        return locateDirectionHemisphere(di, dummyx, dummyy);
    }

    inline int& getHemisphereMap(int nx, int ny, int x, int y) {
        return m_hemishpereMap[(((nx * m_angleResolution) + ny) * m_angleResolution + x) * m_angleResolution + y];
    }

    friend class QTableVisualizationIntegrator;

    int m_sceneResolution;
    int m_angleResolution;
    float m_alpha;
    bool m_useVisit = true;
    Vector3f m_sceneBlockSize;
    BoundingBox3f m_sceneBox;
    WrapperMap m_storage;
    int* m_hemishpereMap;
    const float UPDATE_THREASHOLD = 0.1f;
    std::string m_importFilename;
    std::string m_exportFilename;
};


class QTableVisualizationIntegrator : public Integrator {
public:
    QTableVisualizationIntegrator(const PropertyList &props) {
        /* No parameters this time */
    }

    void addChild(NoriObject *obj) {
        switch (obj->getClassType()) {
        case EGuider:
            if (m_guider)
                throw NoriException("There can only be one guider per integrator!");
            m_guider = dynamic_cast<QTableSphereGuider*>(obj);
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
        Intersection its, its_;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);

        if (!its.mesh->getBSDF()->isProbe())
            return Color3f(0.0f);

        if (!scene->rayIntersect(Ray3f(its.p, -its.shFrame.n), its_) || its_.mesh->getBSDF()->isProbe())
            return Color3f(0.0f);
        QTableSphereGuider::WrapperMap::const_accessor const_access;
        int block_idx = m_guider->locateBlock(its_.p);
        int nx, ny, ox, oy;
        m_guider->locateDirection(its_.shFrame.n, nx, ny);
        m_guider->locateDirectionHemisphere(its_.shFrame.toLocal(its.shFrame.n), ox, oy);
        int angle_idx = m_guider->getHemisphereMap(nx, ny, ox, oy);

        if (m_guider->m_storage.find(const_access, block_idx)) {
            float maxq = 0.0f;

            for (int i = 0; i < m_guider->m_angleResolution; i++) {
                for (int j = 0; j < m_guider->m_angleResolution; j++) {
                    maxq = std::max(maxq, const_access->second.map[m_guider->getHemisphereMap(nx, ny, i, j)]);
                }
            }
            return Color3f(const_access->second.map[angle_idx] / maxq, 1.0f - std::min(1.0f, const_access->second.map[angle_idx] / maxq), 0.0f);
        }
        return Color3f(0.0f);
    }

    std::string toString() const {
        return "QTableVisualizationIntegrator[]";
    }

protected:
    QTableSphereGuider* m_guider = nullptr;
};

NORI_REGISTER_CLASS(QTableSphereGuider, "qtable_sphere");
NORI_REGISTER_CLASS(QTableVisualizationIntegrator, "qtable_visualization");
NORI_NAMESPACE_END
