project('areafilter', 'cpp',
  default_options: ['buildtype=release', 'warning_level=2', 'b_lto=true', 'b_ndebug=if-release', 'cpp_std=c++20'],
  license: 'GPL-3.0-or-later',
  meson_version: '>=0.51.0',
)

cxx = meson.get_compiler('cpp')

gcc_syntax = cxx.get_argument_syntax() == 'gcc'

if gcc_syntax
  vapoursynth_dep = dependency('vapoursynth', version: '>=55').partial_dependency(compile_args: true, includes: true)
  install_dir = vapoursynth_dep.get_variable(pkgconfig: 'libdir') / 'vapoursynth'
else
  vapoursynth_dep = dependency('', required: false)
  install_dir = get_option('libdir') / 'vapoursynth'
endif


sources = [
  'areafilter/areafilter.cpp'
]

shared_module('areafilter', sources,
  dependencies: [vapoursynth_dep],
  install: true,
  install_dir: install_dir,
  gnu_symbol_visibility: 'hidden'
)