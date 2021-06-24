#include <functional>
#include <tracer/guider.h>
#include <tracer/scene.h>
#include <tracer/sampler.h>
#include <tracer/emitter.h>
#include <tracer/bsdf.h>
#include <tracer/warp.h>
#include <tbb/spin_rw_mutex.h>
#include <tbb/concurrent_vector.h>
#include <tbb/concurrent_hash_map.h>

TRACER_NAMESPACE_BEGIN

class QTableGuider : public Guider {
protected:
    template<typename Scalar>
    struct RangeTree;
    struct Wrapper {
        RangeTree<float>* tree;
        int* visit;

        ~Wrapper() {
            if (tree)
                delete tree;
            if (visit)
                delete[] visit;
        }

        void init(int width, int height) {
            tree = new RangeTree<float>(width, height);
            visit = new int[width * height];
            memset(visit, 0, width * height * sizeof(int));
        }
    };
    typedef tbb::concurrent_hash_map<int, Wrapper> WrapperMap;
public:
    QTableGuider(const PropertyList &props): m_storage(10000) {
        m_sceneResolution = props.getInteger("sceneResolution", 50);
        m_angleResolution = props.getInteger("angleResolution", 8);
        try {
            m_alpha = props.getFloat("alpha");
            m_useVisit = false;
        }
        catch (TracerException e) {
            //alpha doesn't exist
        }
    }

    /* Integrator need to call this in preprocess() */
    void init(const Scene *scene) {
        m_sceneBox = scene->getBoundingBox();
        Point3f orig_max = m_sceneBox.max;
        m_sceneBox.expandBy(orig_max + Vector3f(Epsilon));
        m_sceneBlockSize = (m_sceneBox.max - m_sceneBox.min) / m_sceneResolution;
    }

    Vector3f sample(const Point2f& sample, const Intersection& its, float& pdf) {
        int block_idx = locateBlock(its.p);
        Point2f result;
        WrapperMap::const_accessor const_access;
        WrapperMap::accessor access;
        if (m_storage.find(const_access, block_idx)) {
            result = const_access->second.tree->warp(sample, pdf);
        }
        else {
            if (m_storage.insert(access, block_idx)) {
                access->second.init(m_angleResolution, m_angleResolution);
            }
            result = access->second.tree->warp(sample, pdf);
        }
        pdf *= INV_TWOPI;
        return Warp::squareToUniformHemisphere(result);
    }

    void update(const Intersection& origin, const Intersection& dest, Sampler* sampler) {
        const Vector3f& ray = (dest.p - origin.p).normalized(),
                origin_wo = origin.shFrame.toLocal(ray),
                dest_wi = dest.shFrame.toLocal(-ray);
        int ox, oy;
        int block_orig_idx = locateBlock(origin.p), angle_orig_idx = locateDirection(origin_wo, ox, oy);
        int block_dest_idx = locateBlock(dest.p);

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
                        float normal_q = accessor->second.tree->get(i, j);
                        float term = normal_q * Frame::cosTheta(brec.wo) * eval;
                        integral_term += term;
                    }
                }
            }
            else {
                //We need to sample the incident ray because brdf is always 0
                for (int i = 0; i < m_angleResolution * m_angleResolution; i++) {
                    bsdf->sample(brec, sampler->next2D());
                    int tx, ty;
                    locateDirection(brec.wo, tx, ty);
                    integral_term += accessor->second.tree->get(tx, ty);
                }
            }
            accessor.release();
        };
        if (m_storage.find(const_access_dest, block_dest_idx)) {
            integral_job(const_access_dest);
        }
        else {
            if (m_storage.insert(access_dest, block_dest_idx)) {
                access_dest->second.init(m_angleResolution, m_angleResolution);
            }
            integral_job(access_dest);
        }
        integral_term *= 2.0f * M_PI / m_angleResolution / m_angleResolution;
        if (dest.mesh->isEmitter()) {
            integral_term += dest.mesh->getEmitter()->getRadiance(dest.p, dest_wi).sum();
        }


        if (m_storage.insert(access_orig, block_orig_idx)) {
            access_orig->second.init(m_angleResolution, m_angleResolution);
        }
        float alpha = m_useVisit ? 1.0f / (1 + access_orig->second.visit[angle_orig_idx]) : m_alpha;
        float oldval = access_orig->second.tree->get(ox, oy);
        float newval = (1.0f - alpha) * oldval + alpha * integral_term;
        access_orig->second.tree->update(ox, oy, newval);
        access_orig->second.visit[angle_orig_idx]++;
    }

    float pdf(const Vector3f& di, const Intersection& origin) {
        int ox, oy;
        int block_idx = locateBlock(origin.p), angle_idx = locateDirection(di, ox, oy);
        WrapperMap::const_accessor const_access;
        WrapperMap::accessor access;
        if (m_storage.find(const_access, block_idx)) {
            return const_access->second.tree->getPdf(ox, oy);
        }
        else {
            if (m_storage.insert(access, block_idx)) {
                access->second.init(m_angleResolution, m_angleResolution);
            }
            return access->second.tree->getPdf(ox, oy);
        }
        return 0.0f;
    }

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
        if (di.z() < 1 - 1e-6f) {
            y = sphericalCoordinates(di).y();
        }
        y /= 2 * M_PI;
        ix = x * m_angleResolution;
        iy = y * m_angleResolution;
        return ix * m_angleResolution + iy;
    }

    int locateDirection(const Vector3f& di) const {
        int dummyx, dummyy;
        return locateDirection(di, dummyx, dummyy);
    }

	std::string toString() const {
        return tfm::format(
            "QTableGuider[\n"
            "  alpha = %s,\n"
            "  sceneResolution = %d,\n"
            "  angleResolution = %d\n"
            "]",
            m_useVisit? "1/(1 + visit)" : tfm::format("%f", m_alpha), m_sceneResolution, m_angleResolution);
	}

protected:
    /* A simple 2D Range Tree with for log(n) update and sample. */
    /* This is not thread safe so synchronization should be done on top */
    template<typename Scalar>
    struct RangeTree {
        struct Node {
            struct Node* left;
            struct Node* right;
            Scalar sum;
        };

        RangeTree(int width, int height, Scalar initv = 1.0f) : RangeTree(width, height, [initv](int, int) -> Scalar { return initv; }) { }

        RangeTree(int width, int height, const std::function<Scalar(int, int)>& initializer) : m_width(width), m_height(height), m_data(new Scalar[width * height]), m_yroots(new Node*[width]) {
            //Build all the columns
            Scalar *xsum = new Scalar[width];
            for (int i = 0; i < m_width; i++) {
                xsum[i] = build1D(&m_yroots[i], 0, height, [=, &initializer](int y) -> Scalar { Scalar v = initializer(i, y); m_data[i * width + y] = v; return v; });
            }
            build1D(&m_xroot, 0, width, [xsum](int x) -> Scalar { return xsum[x]; });
            delete[] xsum;
        }

        ~RangeTree() {
            delete[] m_data;
            for (int i = 0; i < m_width; i++) {
                free1D(&m_yroots[i]);
            }
            delete[] m_yroots;
            free1D(&m_xroot);
        }

        Scalar build1D(Node **root, int l, int r, const std::function<Scalar(int)>& initializer) {
            if (r == l + 1) {
                //Leaf
                *root = new Node{ nullptr, nullptr, initializer(l) };
            }
            else {
                *root = new Node;
                int m = (l + r) >> 1;
                (*root)->sum = build1D(&(*root)->left, l, m, initializer) + build1D(&(*root)->right, m, r, initializer);
            }
            return (*root)->sum;
        }

        TPoint<Scalar, 2> warp(const TPoint<Scalar, 2>& _sample, float& pdf) {
            TPoint<Scalar, 2> sample = _sample;
            int l = 0, r = m_width, yl = 0, yr = m_height;
            Node *root = m_xroot;
            Scalar total_sum = root->sum;
            while (l + 1 < r) {
                Scalar breakdown = root->left->sum / root->sum;
                if (sample.x() < breakdown) {
                    //left
                    sample.x() /= breakdown;
                    root = root->left;
                    r = (l + r) >> 1;
                }
                else {
                    sample.x() = (sample.x() - breakdown) / ((Scalar)1 - breakdown);
                    root = root->right;
                    l = (l + r) >> 1;
                }
            }
            root = m_yroots[l];
            while (yl + 1 < yr) {
                Scalar breakdown = root->left->sum / root->sum;
                if (sample.y() < breakdown) {
                    //left
                    sample.y() /= breakdown;
                    root = root->left;
                    yr = (yl + yr) >> 1;
                }
                else {
                    sample.y() = (sample.y() - breakdown) / ((Scalar)1 - breakdown);
                    root = root->right;
                    yl = (yl + yr) >> 1;
                }
            }
            pdf = root->sum / total_sum * m_width * m_height;
            return TPoint<Scalar, 2>(((Scalar)l + sample.x()) / m_width, ((Scalar)yl + sample.y()) / m_height);
        }

        Scalar rebuildX(Node *root, int l, int r) {
            if (r == l + 1) {
                root->sum = rebuildY(m_yroots[l], l * m_width, 0, m_height);
            }
            else {
                int m = (l + r) >> 1;
                root->sum = rebuildX(root->left, l, m) + rebuildX(root->right, m, r);
            }
            return root->sum;
        }

        Scalar rebuildY(Node *root, int x, int l, int r) {
            if (r == l + 1) {
                root->sum = m_data[x + l];
            }
            else {
                int m = (l + r) >> 1;
                root->sum = rebuildY(root->left, x, l, m) + rebuildY(root->right, x, m, r);
            }
            return root->sum;
        }

        void update(int i, int j, Scalar newval) {
            if(newval < WEIGHT_THREASHOLD)
                newval = WEIGHT_THREASHOLD;
            Scalar diff = newval - get(i, j);
            int l = 0, r = m_width;
            Node *root = m_xroot;
            m_data[i * m_width + j] = newval;
            while (root) {
                root->sum += diff;
                int m = (l + r) >> 1;
                if (i < m) {
                    root = root->left;
                    r = m;
                }
                else {
                    root = root->right;
                    l = m;
                }
            }
            root = m_yroots[l];
            l = 0; r = m_height;
            while (root) {
                root->sum += diff;
                int m = (l + r) >> 1;
                if (j < m) {
                    root = root->left;
                    r = m;
                }
                else {
                    root = root->right;
                    l = m;
                }
            }
        }

        inline Scalar get(int i, int j) const {
            return m_data[i * m_width + j];
        }

        inline Scalar getPdf(int i, int j) const {
            return m_data[i * m_width + j] / m_xroot->sum * m_width * m_height * INV_TWOPI;
        }

        void free1D(Node **root) {
            if (*root) {
                free1D(&(*root)->left);
                free1D(&(*root)->right);
                delete root;
            }
        }
        
        Scalar *m_data = nullptr;
        Node *m_xroot = nullptr;
        Node **m_yroots = nullptr;
        int m_width, m_height;
        const float WEIGHT_THREASHOLD = 0.1f;
    };
    int m_sceneResolution;
    int m_angleResolution;
    float m_alpha;
    bool m_useVisit = true;
    Vector3f m_sceneBlockSize;
    BoundingBox3f m_sceneBox;
    WrapperMap m_storage;

};

TRACER_REGISTER_CLASS(QTableGuider, "qtable");
TRACER_NAMESPACE_END