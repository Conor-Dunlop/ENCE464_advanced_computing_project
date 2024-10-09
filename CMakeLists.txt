# CMakeList.txt : CMake project for assignment_2, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.13)

# Enable Hot Reload for MSVC compilers if supported.
if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

add_compile_options(-mavx512f -march=znver4)

project ("poisson")

# Add source to this project's executable.
add_executable (poisson "main.cpp")
target_include_directories(poisson PRIVATE ./include)
target_sources(poisson PRIVATE "poisson_mt.cpp")
set_property(TARGET poisson PROPERTY CXX_STANDARD 20)


# TODO: Add tests and install targets if needed.
