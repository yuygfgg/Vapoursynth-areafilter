# VapourSynth AreaFilter Plugin

VapourSynth plugin for processing binary images by filtering connected components based on area size.

## Usage

### AreaFilter

Removes connected components with an area below a specified threshold:

```python
core.areafilter.Filter(clip clip, int min_area, neighbors8=False)
```

Parameters:
- `clip`: Input clip (binary image)
- `min_area`: Minimum area threshold in pixels
- `neighbors8`: Use 8-neighborhood connectivity when True, 4-neighborhood when False (default: False)

### RelFilter

Keeps only the largest connected components up to a specified percentage of the total area:

```python
core.areafilter.RelFilter(clip clip, int percentage, neighbors8=False)
```

Parameters:
- `clip`: Input clip (binary image)
- `percentage`: Percentage of largest components to keep (1-100)
- `neighbors8`: Use 8-neighborhood connectivity when True, 4-neighborhood when False (default: False)

## Frame Properties

Both filters output the following frame properties:
- `ComponentCount`: Number of connected components
- `SizePercentile0`, `SizePercentile5`, ..., `SizePercentile100`: Component size percentiles

## Building

```bash
meson build
ninja -C build install
```
