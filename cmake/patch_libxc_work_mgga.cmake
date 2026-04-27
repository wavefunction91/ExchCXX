# patch_libxc_work_mgga.cmake
# Portable replacement for:
#   sed -i -e "s/p->info->family != XC_KINETIC/p->info->kind != XC_KINETIC/g" src/work_mgga_inc.c
# Run from the libxc source directory via PATCH_COMMAND.

set(_file "src/work_mgga_inc.c")
if(NOT EXISTS "${_file}")
  message(FATAL_ERROR "patch_libxc_work_mgga.cmake: ${_file} not found (cwd=${CMAKE_CURRENT_SOURCE_DIR})")
endif()

file(READ "${_file}" _content)
string(REPLACE "p->info->family != XC_KINETIC" "p->info->kind != XC_KINETIC" _content "${_content}")
file(WRITE "${_file}" "${_content}")
