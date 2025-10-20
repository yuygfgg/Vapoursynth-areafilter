#include "VSHelper4.h"
#include "VapourSynth4.h"
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <format>
#include <vector>

typedef struct {
    int component_count;
    std::vector<int> size_percentiles;
    std::vector<int> component_sizes;
} ComponentStats;

typedef ComponentStats (*ProcessPlaneFn)(const void*, void*, int, int,
                                         ptrdiff_t, ptrdiff_t, int, float,
                                         float);

typedef struct {
    VSNode* node;
    union {
        int min_area;
        float percentage;
    };
    VSSampleType sample_type;
    int bits_per_sample;
    uint16_t max_value;
    float fg_value;
    bool write_props;
    ProcessPlaneFn process_plane_fn;
} FilterData;

class DisjointSet {
  private:
    std::vector<int> parent;
    std::vector<int> size;

  public:
    DisjointSet(int max_elements)
        : parent(max_elements), size(max_elements, 1) {
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

struct NeighborOffset {
    int dy, dx;
};

// clang-format off
constexpr NeighborOffset EIGHT_NEIGHBORS[] = {
    {-1, -1}, {-1, 0}, {-1, 1}, 
    { 0, -1},                      { 0, 1},
    { 1, -1}, { 1, 0}, { 1, 1}
};
// clang-format on

constexpr auto EIGHT_NEIGHBORS_COUNT = 8;

// clang-format off
constexpr NeighborOffset FOUR_NEIGHBORS[] = {
                        {-1, 0}, 
    {0, -1},                    {0, 1}, 
                        {1, 0}
};
// clang-format on

constexpr auto FOUR_NEIGHBORS_COUNT = 4;

template <auto use_8_neighbors> struct NeighborhoodTraits;

template <> struct NeighborhoodTraits<true> {
    static constexpr auto neighbors = EIGHT_NEIGHBORS;
    static constexpr auto count = EIGHT_NEIGHBORS_COUNT;
};

template <> struct NeighborhoodTraits<false> {
    static constexpr auto neighbors = FOUR_NEIGHBORS;
    static constexpr auto count = FOUR_NEIGHBORS_COUNT;
};

template <bool use_8_neighbors, bool use_percentage, typename T>
static inline auto processPlane(const T* VS_RESTRICT srcp, T* VS_RESTRICT dstp,
                                auto width, auto height, auto src_stride,
                                auto dst_stride, auto min_area, auto fg_value,
                                auto percentage = 0.0f) noexcept {
    auto src_stride_elements = src_stride / sizeof(T);
    auto dst_stride_elements = dst_stride / sizeof(T);

    std::vector<int> labels(width * height, 0);
    DisjointSet ds(width * height + 2);

    constexpr auto neighbors = NeighborhoodTraits<use_8_neighbors>::neighbors;
    constexpr auto num_neighbors = NeighborhoodTraits<use_8_neighbors>::count;

    auto next_label = 1;

    for (auto y = 0; y < height; y++) {
        for (auto x = 0; x < width; x++) {
            if (srcp[y * src_stride_elements + x] != fg_value)
                continue;

            auto min_label = 0;
            std::vector<int> valid_neighbors;

            for (auto i = 0; i < num_neighbors; i++) {
                auto dy = neighbors[i].dy;
                auto dx = neighbors[i].dx;
                auto ny = y + dy;
                auto nx = x + dx;

                if (ny >= 0 && nx >= 0 && ny < height && nx < width &&
                    srcp[ny * src_stride_elements + nx] == fg_value) {
                    auto neighbor_label = labels[ny * width + nx];
                    if (neighbor_label > 0) {
                        valid_neighbors.push_back(neighbor_label);
                        if (min_label == 0 || neighbor_label < min_label)
                            min_label = neighbor_label;
                    }
                }
            }

            if (min_label == 0) {
                labels[y * width + x] = next_label++;
            } else {
                labels[y * width + x] = min_label;

                for (auto nl : valid_neighbors) {
                    if (nl != min_label)
                        ds.merge(min_label, nl);
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

    ComponentStats stats;
    std::vector<int> non_zero_sizes;
    for (auto i = 1; i <= max_label; i++) {
        if (component_sizes[i] > 0) {
            non_zero_sizes.push_back(component_sizes[i]);
        }
    }

    stats.component_count = non_zero_sizes.size();

    stats.component_sizes = non_zero_sizes;

    stats.size_percentiles.resize(21);

    if (!non_zero_sizes.empty()) {
        std::sort(non_zero_sizes.begin(), non_zero_sizes.end());

        for (auto i = 0; i <= 20; i++) {
            auto percentile = i * 5.0f;
            auto idx = static_cast<int>(
                (percentile / 100.0f) * (non_zero_sizes.size() - 1) + 0.5f);
            idx = std::min(std::max(0, idx),
                           static_cast<int>(non_zero_sizes.size() - 1));
            stats.size_percentiles[i] = non_zero_sizes[idx];
        }
    } else {
        for (auto i = 0; i <= 20; i++) {
            stats.size_percentiles[i] = 0;
        }
    }

    for (auto y = 0; y < height; y++) {
        auto row = reinterpret_cast<T*>(reinterpret_cast<uint8_t*>(dstp) +
                                        y * dst_stride);
        std::memset(row, 0, width * sizeof(T));
    }

    auto size_threshold = 0;

    if constexpr (use_percentage) {
        if (!non_zero_sizes.empty()) {
            std::sort(non_zero_sizes.begin(), non_zero_sizes.end(),
                      std::greater<int>());

            auto total_area = 0;
            for (auto size : non_zero_sizes) {
                total_area += size;
            }

            auto area_to_keep =
                static_cast<int>(total_area * percentage / 100.0f + 0.5f);
            auto current_area = 0;

            for (auto size : non_zero_sizes) {
                current_area += size;
                size_threshold = size;
                if (current_area >= area_to_keep) {
                    break;
                }
            }
        }
    }

    for (auto y = 0; y < height; y++) {
        for (auto x = 0; x < width; x++) {
            auto label = labels[y * width + x];
            if (label > 0) {
                auto component_size = component_sizes[ds.find(label)];
                auto keep = false;

                if constexpr (use_percentage) {
                    keep = (component_size >= size_threshold);
                } else {
                    keep = (component_size >= min_area);
                }

                if (keep) {
                    dstp[y * dst_stride_elements + x] = fg_value;
                }
            }
        }
    }

    return stats;
}

template <bool use_8_neighbors, bool use_percentage, typename T>
static inline auto
processPlaneWrapper(const void* srcp, void* dstp, int width, int height,
                    ptrdiff_t src_stride, ptrdiff_t dst_stride, int min_area,
                    float fg_value, float percentage) noexcept {
    return processPlane<use_8_neighbors, use_percentage, T>(
        static_cast<const T*>(srcp), static_cast<T*>(dstp), width, height,
        src_stride, dst_stride, min_area, static_cast<T>(fg_value),
        use_percentage ? percentage : 0.0f);
}

static inline auto setFrameProperties(auto dst, auto stats, auto vsapi) {
    vsapi->mapSetInt(vsapi->getFramePropertiesRW(dst), "ComponentCount",
                     stats.component_count, maReplace);

    for (auto i = 0; i <= 20; i++) {
        auto propName = std::format("SizePercentile{}", i * 5);
        vsapi->mapSetInt(vsapi->getFramePropertiesRW(dst), propName.c_str(),
                         stats.size_percentiles[i], maReplace);
    }
}

static inline const VSFrame* VS_CC
areaFilterGetFrame(auto n, auto activationReason, auto instanceData,
                   [[maybe_unused]] auto frameData, auto frameCtx, auto core,
                   auto vsapi) noexcept {
    auto d = static_cast<FilterData*>(instanceData);

    if (activationReason == arInitial) {
        vsapi->requestFrameFilter(n, d->node, frameCtx);
    } else if (activationReason == arAllFramesReady) {
        auto src = vsapi->getFrameFilter(n, d->node, frameCtx);
        auto fi = vsapi->getVideoFrameFormat(src);
        auto height = vsapi->getFrameHeight(src, 0);
        auto width = vsapi->getFrameWidth(src, 0);

        auto dst = vsapi->newVideoFrame(fi, width, height, src, core);

        std::vector<ComponentStats> plane_stats;

        for (auto plane = 0; plane < fi->numPlanes; plane++) {
            const void* srcp = vsapi->getReadPtr(src, plane);
            auto src_stride = vsapi->getStride(src, plane);
            void* dstp = vsapi->getWritePtr(dst, plane);
            auto dst_stride = vsapi->getStride(dst, plane);

            auto plane_width = vsapi->getFrameWidth(src, plane);
            auto plane_height = vsapi->getFrameHeight(src, plane);

            auto stats = d->process_plane_fn(
                srcp, dstp, plane_width, plane_height, src_stride, dst_stride,
                d->min_area, d->fg_value, 0.0f);

            plane_stats.push_back(stats);
        }

        if (d->write_props && !plane_stats.empty()) {
            setFrameProperties(dst, plane_stats[0], vsapi);
        }

        vsapi->freeFrame(src);

        return dst;
    }
    return nullptr;
}

static inline const VSFrame* VS_CC
relFilterGetFrame(auto n, auto activationReason, auto instanceData,
                  [[maybe_unused]] auto frameData, auto frameCtx, auto core,
                  auto vsapi) noexcept {
    auto d = static_cast<FilterData*>(instanceData);

    if (activationReason == arInitial) {
        vsapi->requestFrameFilter(n, d->node, frameCtx);
    } else if (activationReason == arAllFramesReady) {
        auto src = vsapi->getFrameFilter(n, d->node, frameCtx);
        auto fi = vsapi->getVideoFrameFormat(src);
        auto height = vsapi->getFrameHeight(src, 0);
        auto width = vsapi->getFrameWidth(src, 0);

        auto dst = vsapi->newVideoFrame(fi, width, height, src, core);

        std::vector<ComponentStats> plane_stats;

        for (auto plane = 0; plane < fi->numPlanes; plane++) {
            const void* srcp = vsapi->getReadPtr(src, plane);
            auto src_stride = vsapi->getStride(src, plane);
            void* dstp = vsapi->getWritePtr(dst, plane);
            auto dst_stride = vsapi->getStride(dst, plane);

            auto plane_width = vsapi->getFrameWidth(src, plane);
            auto plane_height = vsapi->getFrameHeight(src, plane);

            auto stats = d->process_plane_fn(
                srcp, dstp, plane_width, plane_height, src_stride, dst_stride,
                0, d->fg_value, d->percentage);

            plane_stats.push_back(stats);
        }

        if (d->write_props && !plane_stats.empty()) {
            setFrameProperties(dst, plane_stats[0], vsapi);
        }

        vsapi->freeFrame(src);

        return dst;
    }
    return nullptr;
}

static inline auto VS_CC filterFree(auto instanceData,
                                    [[maybe_unused]] auto core,
                                    auto vsapi) noexcept {
    auto d = static_cast<FilterData*>(instanceData);
    vsapi->freeNode(d->node);
    free(d);
}

static inline auto validateInput(auto in, auto out, auto vsapi, auto& d,
                                 auto filter_name) noexcept {
    d.node = vsapi->mapGetNode(in, "clip", 0, 0);
    auto vi = vsapi->getVideoInfo(d.node);

    if (!vsh::isConstantVideoFormat(vi)) {
        vsapi->mapSetError(
            out, std::format("{}: only clips with constant format are accepted",
                             filter_name)
                     .c_str());
        vsapi->freeNode(d.node);
        return false;
    }

    if (!((vi->format.bitsPerSample >= 8 && vi->format.bitsPerSample <= 16 &&
           vi->format.sampleType == stInteger) ||
          (vi->format.bitsPerSample == 32 &&
           vi->format.sampleType == stFloat))) {
        vsapi->mapSetError(
            out, std::format("{}: only 8-16 bit integer or 32 bit float input "
                             "are accepted, got {} bit {}",
                             filter_name, vi->format.bitsPerSample,
                             vi->format.sampleType == stInteger ? "integer"
                                                                : "float")
                     .c_str());
        vsapi->freeNode(d.node);
        return false;
    }

    return true;
}

static inline void setupCommonFilterData(auto& d, auto vsapi) {
    auto vi = vsapi->getVideoInfo(d.node);
    d.sample_type = static_cast<VSSampleType>(vi->format.sampleType);
    d.bits_per_sample = vi->format.bitsPerSample;

    if (d.sample_type == stInteger) {
        d.max_value = static_cast<uint16_t>((1 << d.bits_per_sample) - 1);
    } else {
        d.max_value = 0;
    }

    if (d.sample_type == stInteger) {
        if (d.bits_per_sample == 8) {
            d.fg_value = 255.0f;
        } else {
            d.fg_value = static_cast<float>(d.max_value);
        }
    } else {
        d.fg_value = 1.0f;
    }
}

static inline void selectProcessFunction(auto& d, auto use_8_neighbors,
                                         auto use_percentage) {
    if (d.sample_type == stInteger) {
        if (d.bits_per_sample == 8) {
            if (use_8_neighbors) {
                d.process_plane_fn =
                    use_percentage ? processPlaneWrapper<true, true, uint8_t>
                                   : processPlaneWrapper<true, false, uint8_t>;
            } else {
                d.process_plane_fn =
                    use_percentage ? processPlaneWrapper<false, true, uint8_t>
                                   : processPlaneWrapper<false, false, uint8_t>;
            }
        } else {
            if (use_8_neighbors) {
                d.process_plane_fn =
                    use_percentage ? processPlaneWrapper<true, true, uint16_t>
                                   : processPlaneWrapper<true, false, uint16_t>;
            } else {
                d.process_plane_fn =
                    use_percentage
                        ? processPlaneWrapper<false, true, uint16_t>
                        : processPlaneWrapper<false, false, uint16_t>;
            }
        }
    } else {
        if (use_8_neighbors) {
            d.process_plane_fn = use_percentage
                                     ? processPlaneWrapper<true, true, float>
                                     : processPlaneWrapper<true, false, float>;
        } else {
            d.process_plane_fn = use_percentage
                                     ? processPlaneWrapper<false, true, float>
                                     : processPlaneWrapper<false, false, float>;
        }
    }
}

static inline auto VS_CC areaFilterCreate(const VSMap* in, VSMap* out,
                                          [[maybe_unused]] void* userData,
                                          VSCore* core,
                                          const VSAPI* vsapi) noexcept {
    FilterData d;
    FilterData* data;
    auto err = 0;

    constexpr auto filter_name = "AreaFilter";

    if (!validateInput(in, out, vsapi, d, filter_name)) {
        return;
    }

    setupCommonFilterData(d, vsapi);

    d.min_area = vsapi->mapGetInt(in, "min_area", 0, &err);
    if (err) {
        vsapi->mapSetError(
            out, std::format("{}: min_area must be set", filter_name).c_str());
        vsapi->freeNode(d.node);
        return;
    }

    if (d.min_area <= 0) {
        vsapi->mapSetError(
            out, std::format("{}: min_area must be greater than 0, got {}",
                             filter_name, d.min_area)
                     .c_str());
        vsapi->freeNode(d.node);
        return;
    }

    auto use_8_neighbors = !!vsapi->mapGetInt(in, "neighbors8", 0, &err);
    if (err)
        use_8_neighbors = false;

    d.write_props = !!vsapi->mapGetInt(in, "write_props", 0, &err);
    if (err)
        d.write_props = true;

    selectProcessFunction(d, use_8_neighbors, false);

    data = static_cast<FilterData*>(malloc(sizeof(d)));
    *data = d;

    VSFilterDependency deps[] = {{d.node, rpStrictSpatial}};
    vsapi->createVideoFilter(out, filter_name, vsapi->getVideoInfo(d.node),
                             areaFilterGetFrame, filterFree, fmParallel, deps,
                             1, data, core);
}

static inline auto VS_CC relFilterCreate(auto in, auto out,
                                         [[maybe_unused]] auto userData,
                                         auto core, auto vsapi) noexcept {
    FilterData d;
    FilterData* data;
    auto err = 0;

    constexpr auto filter_name = "RelFilter";

    if (!validateInput(in, out, vsapi, d, filter_name)) {
        return;
    }

    setupCommonFilterData(d, vsapi);

    d.percentage =
        static_cast<float>(vsapi->mapGetFloat(in, "percentage", 0, &err));
    if (err) {
        vsapi->mapSetError(out,
                           std::format("{}: percentage must be set, got {}",
                                       filter_name, d.percentage)
                               .c_str());
        vsapi->freeNode(d.node);
        return;
    }

    if (d.percentage <= 0.0f || d.percentage > 100.0f) {
        vsapi->mapSetError(
            out, std::format("{}: percentage must be in the range (0, 100], "
                             "got {}",
                             filter_name, d.percentage)
                     .c_str());
        vsapi->freeNode(d.node);
        return;
    }

    auto use_8_neighbors = !!vsapi->mapGetInt(in, "neighbors8", 0, &err);
    if (err)
        use_8_neighbors = false;

    d.write_props = !!vsapi->mapGetInt(in, "write_props", 0, &err);
    if (err)
        d.write_props = true;

    selectProcessFunction(d, use_8_neighbors, true);

    data = static_cast<FilterData*>(malloc(sizeof(d)));
    *data = d;

    VSFilterDependency deps[] = {{d.node, rpStrictSpatial}};
    vsapi->createVideoFilter(out, filter_name, vsapi->getVideoInfo(d.node),
                             relFilterGetFrame, filterFree, fmParallel, deps, 1,
                             data, core);
}

VS_EXTERNAL_API(void)
VapourSynthPluginInit2(VSPlugin* plugin, const VSPLUGINAPI* vspapi) {
    vspapi->configPlugin("com.yuygfgg.areafilter", "areafilter",
                         "VapourSynth Area Filter Plugin",
                         VS_MAKE_VERSION(3, 0), VAPOURSYNTH_API_VERSION, 0,
                         plugin);
    vspapi->registerFunction(
        "AreaFilter",
        "clip:vnode;min_area:int;neighbors8:int:opt;write_props:int:opt;",
        "clip:vnode;", areaFilterCreate, NULL, plugin);
    vspapi->registerFunction(
        "RelFilter",
        "clip:vnode;percentage:float;neighbors8:int:opt;write_props:int:opt;",
        "clip:vnode;", relFilterCreate, NULL, plugin);
}