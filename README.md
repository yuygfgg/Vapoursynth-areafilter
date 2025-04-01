# Vapoursynth areafilter Plugin

VapourSynth plugin that filters binarized images by removing connected components with an area below a specified threshold.

## Usage

```python
core.areafilter.Filter(clip clip, int min_area)
```

## Building

```bash
meson build
ninja -C build install
```