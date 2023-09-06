exe_args = []
if fc_id == 'gcc'
  exe_args += '-fall-intrinsics'
endif

qvSZP_exe = executable(
  meson.project_name(),
  sources: files(
    'main.f90',
  ),
  dependencies: qvSZP_dep,
  install: install,
  fortran_args: exe_args,
  link_language: 'fortran',
)
