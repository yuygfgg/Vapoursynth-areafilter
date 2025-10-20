# VapourSynth AreaFilter Plugin

VapourSynth plugin for processing binary images by filtering connected components based on area size.

## Usage

### AreaFilter

Removes connected components with an area below a specified threshold:

```python
core.areafilter.AreaFilter(clip clip, int min_area, neighbors8=False, write_props=True)
```

Parameters:
- `clip`: Input clip (binary image)
- `min_area`: Minimum area threshold in pixels
- `neighbors8`: Use 8-neighborhood connectivity when True, 4-neighborhood when False (default: False)
- `write_props`: Write frame properties when True (default: True)

### RelFilter

Keeps only the largest connected components up to a specified percentage of the total area:

```python
core.areafilter.RelFilter(clip clip, int percentage, neighbors8=False, write_props=True)
```

Parameters:
- `clip`: Input clip (binary image)
- `percentage`: Percentage of largest components to keep (1-100)
- `neighbors8`: Use 8-neighborhood connectivity when True, 4-neighborhood when False (default: False)
- `write_props`: Write frame properties when True (default: True)

## Frame Properties

Both filters output the following frame properties (when `write_props=True`):
- `ComponentCount`: Number of connected components
- `SizePercentile0`, `SizePercentile5`, ..., `SizePercentile100`: Component size percentiles

## Building

```bash
meson build
ninja -C build install
```
