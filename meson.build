# This file is part of tblite.
# SPDX-Identifier: LGPL-3.0-or-later
#
# tblite is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# tblite is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with tblite.  If not, see <https://www.gnu.org/licenses/>.

project(
  'qvSZP',
  'fortran',
  version: '2.0',
  license: 'LGPL-3.0-or-later',
  meson_version: '>=0.57.2',
  default_options: [
    'buildtype=debugoptimized',
    'default_library=both',
  ],
)
install = not (meson.is_subproject() and get_option('default_library') == 'static')
has_cc = add_languages('c', required: get_option('api') or get_option('python'), native: false)

# General configuration information
lib_deps = []
inc_dirs = []
subdir('config')

# Collect source of the project
srcs = []
subdir('src')

# Library target
qvSZP_lib = library(
  meson.project_name(),
  sources: srcs,
  version: meson.project_version(),
  dependencies: lib_deps,
  include_directories: inc_dirs,
  install: install,
)

# Export dependency for other projects and test suite
qvSZP_inc = [qvSZP_lib.private_dir_include()]
qvSZP_dep = declare_dependency(
  link_with: qvSZP_lib,
  include_directories: qvSZP_inc,
  dependencies: lib_deps,
  variables: {'includedir': meson.current_source_dir() / 'include'},
)

# Add executable targets
subdir('app')

if install
  install_subdir(
    'include/qvSZP',
    install_dir: get_option('includedir'),
  )

  module_id = meson.project_name() / fc_id + '-' + fc.version()
  meson.add_install_script(
    find_program(files('config'/'install-mod.py')),
    get_option('includedir') / module_id,
  )

  pkg = import('pkgconfig')
  pkg.generate(
    qvSZP_lib,
    description: 'qvSZP basis',
    subdirs: ['', module_id],
  )

endif
