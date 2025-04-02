#include <stdlib.h>
#include <vector>
#include "VapourSynth4.h"
#include "VSHelper4.h"

typedef struct {
    VSNode *node;
    int min_area;
} AreaFilterData;

class DisjointSet {
private:
    std::vector<int> parent;
    std::vector<int> size;

public:
    DisjointSet(int max_elements) : parent(max_elements), size(max_elements, 1) {
        for (auto i = 0; i < max_elements; i++)
            parent[i] = i;
    }
    
    auto find(int x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    }
    
    auto merge(int x, int y) {
        auto root_x = find(x);
        auto root_y = find(y);
        if (root_x == root_y) 
            return;
    
        if (size[root_x] < size[root_y]) {
            parent[root_x] = root_y;
            size[root_y] += size[root_x];
        } else {
            parent[root_y] = root_x;
            size[root_x] += size[root_y];
        }
    }
    auto getSize(int x) { return size[find(x)]; }
};
    
template<typename T>
static auto processPlane(const T* VS_RESTRICT srcp, T* VS_RESTRICT dstp, int width, int height, ptrdiff_t src_stride, ptrdiff_t dst_stride, int min_area, T fg_value) {
    auto src_stride_elements = src_stride / sizeof(T);
    auto dst_stride_elements = dst_stride / sizeof(T);
        
    std::vector<int> labels(width * height, 0);
    DisjointSet ds(width * height + 2);

    struct NeighborOffset { int dy, dx; };
    constexpr NeighborOffset neighbors[] = {{-1, -1}, {-1, 0}, {-1, 1}, {0, -1}};
        
    auto next_label = 1;
    
    for (auto y = 0; y < height; y++) {
        for (auto x = 0; x < width; x++) {
            if (srcp[y * src_stride_elements + x] != fg_value)
                continue;
                
            auto min_label = 0;
            int valid_neighbors[4] = {0};
            auto num_valid = 0;
    
            for (auto& off : neighbors) {
                auto ny = y + off.dy;
                auto nx = x + off.dx;
                if (ny >= 0 && nx >= 0 && ny < height && nx < width && srcp[ny * src_stride_elements + nx] == fg_value) {
                    auto neighbor_label = labels[ny * width + nx];
                    if (neighbor_label > 0) {
                        valid_neighbors[num_valid++] = neighbor_label;
                        if (min_label == 0 || neighbor_label < min_label)
                            min_label = neighbor_label;
                    }
                }
            }
    
            if (min_label == 0) {
                labels[y * width + x] = next_label++;
            } 
            else {
                labels[y * width + x] = min_label;
                for (auto i = 0; i < num_valid; i++) {
                    if (valid_neighbors[i] != min_label)
                        ds.merge(min_label, valid_neighbors[i]);
                }
            }
        }
    }
    
    auto max_label = next_label - 1;
    std::vector<int> component_sizes(max_label + 1, 0);
    for (auto i = 0; i < width * height; i++) {
        if (labels[i] > 0) {
            auto root = ds.find(labels[i]);
            if (root <= max_label)
                component_sizes[root]++;
        }
    }
    
    for (auto y = 0; y < height; y++) {
        auto row = reinterpret_cast<T*>(reinterpret_cast<uint8_t*>(dstp) + y * dst_stride);
        memset(row, 0, width * sizeof(T));
    }
    
    for (auto y = 0; y < height; y++) {
        for (auto x = 0; x < width; x++) {
            auto label = labels[y * width + x];
            if (label > 0 && component_sizes[ds.find(label)] >= min_area)
                dstp[y * dst_stride_elements + x] = fg_value;
        }
    }
}

static const VSFrame *VS_CC areaFilterGetFrame(int n, int activationReason, void *instanceData, [[maybe_unused]] void **frameData, VSFrameContext *frameCtx, VSCore *core, const VSAPI *vsapi) {
    auto d = (AreaFilterData *)instanceData;

    if (activationReason == arInitial) {
        vsapi->requestFrameFilter(n, d->node, frameCtx);
    } else if (activationReason == arAllFramesReady) {
        auto src = vsapi->getFrameFilter(n, d->node, frameCtx);
        auto fi = vsapi->getVideoFrameFormat(src);
        auto height = vsapi->getFrameHeight(src, 0);
        auto width = vsapi->getFrameWidth(src, 0);

        auto dst = vsapi->newVideoFrame(fi, width, height, src, core);

        for (auto plane = 0; plane < fi->numPlanes; plane++) {
            const void *srcp = vsapi->getReadPtr(src, plane);
            auto src_stride = vsapi->getStride(src, plane);
            void *dstp = vsapi->getWritePtr(dst, plane);
            auto dst_stride = vsapi->getStride(dst, plane);
            
            auto plane_width = vsapi->getFrameWidth(src, plane);
            auto plane_height = vsapi->getFrameHeight(src, plane);

            if (fi->sampleType == stInteger) {
                if (fi->bitsPerSample == 8) {
                    processPlane<uint8_t>(
                        static_cast<const uint8_t*>(srcp),
                        static_cast<uint8_t*>(dstp),
                        plane_width, plane_height, src_stride, dst_stride,
                        d->min_area, static_cast<uint8_t>(255)
                    );
                } else{
                    auto max_value = static_cast<uint16_t>((1 << fi->bitsPerSample) - 1);
                    processPlane<uint16_t>(
                        static_cast<const uint16_t*>(srcp),
                        static_cast<uint16_t*>(dstp),
                        plane_width, plane_height, src_stride, dst_stride,
                        d->min_area, max_value
                    );
                }
            } else if (fi->sampleType == stFloat) {
                processPlane<float>(
                    static_cast<const float*>(srcp),
                    static_cast<float*>(dstp),
                    plane_width, plane_height, src_stride, dst_stride,
                    d->min_area, 1.0f
                );
            }
        }

        vsapi->freeFrame(src);
        
        return dst;
    }
    return NULL;
}

static auto VS_CC areaFilterFree(void *instanceData, [[maybe_unused]] VSCore *core, const VSAPI *vsapi) {
    auto d = (AreaFilterData *)instanceData;
    vsapi->freeNode(d->node);
    free(d);
}

static auto VS_CC areaFilterCreate(const VSMap *in, VSMap *out, [[maybe_unused]] void *userData, VSCore *core, const VSAPI *vsapi) {
    AreaFilterData d;
    AreaFilterData *data;
    auto err = 0;

    d.node = vsapi->mapGetNode(in, "clip", 0, 0);
    auto vi = vsapi->getVideoInfo(d.node);

    if (!vsh::isConstantVideoFormat(vi)) {
        vsapi->mapSetError(out, "AreaFilter: only clips with constant format are accepted");
        vsapi->freeNode(d.node);
        return;
    }

    if (!(((vi->format.bitsPerSample == 8 || vi->format.bitsPerSample == 16) && vi->format.sampleType == stInteger) || (vi->format.bitsPerSample==32 && vi->format.sampleType == stFloat))){
        vsapi->mapSetError(out, "AreaFilter: only 8-16 bit integer or 32 bit float input are accepted");
        vsapi->freeNode(d.node);
        return;
    }

    d.min_area = vsapi->mapGetInt(in, "min_area", 0, &err);
    if (err) {
        vsapi->mapSetError(out, "AreaFilter: min_area must be set");
        vsapi->freeNode(d.node);
        return;
    }
    
    if (d.min_area <= 0) {
        vsapi->mapSetError(out, "AreaFilter: min_area must be greater than 0");
        vsapi->freeNode(d.node);
        return;
    }

    data = (AreaFilterData *)malloc(sizeof(d));
    *data = d;
    
    VSFilterDependency deps[] = {{d.node, rpStrictSpatial}};
    vsapi->createVideoFilter(out, "AreaFilter", vi, areaFilterGetFrame, areaFilterFree, fmParallel, deps, 1, data, core);
}

VS_EXTERNAL_API(void) VapourSynthPluginInit2(VSPlugin *plugin, const VSPLUGINAPI *vspapi) {
    vspapi->configPlugin("com.yuygfgg.areafilter", "areafilter", "VapourSynth Area Filter Plugin", VS_MAKE_VERSION(1, 0), VAPOURSYNTH_API_VERSION, 0, plugin);
    vspapi->registerFunction("Filter", "clip:vnode;min_area:int;", "clip:vnode;", areaFilterCreate, NULL, plugin);
}