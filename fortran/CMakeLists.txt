cmake_minimum_required(VERSION 3.18)

project(gwr_parameteridentifikation_fortran LANGUAGES Fortran)

set(CMAKE_Fortran_STANDARD 95)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)
set(CMAKE_Fortran_EXTENSIONS OFF)

set(GWR_BLAS_LAPACK_LIBRARY "" CACHE FILEPATH
    "Optional explicit combined BLAS/LAPACK library (for reproducible benchmarks)")

if(GWR_BLAS_LAPACK_LIBRARY)
    if(NOT EXISTS "${GWR_BLAS_LAPACK_LIBRARY}")
        message(FATAL_ERROR
            "GWR_BLAS_LAPACK_LIBRARY does not exist: ${GWR_BLAS_LAPACK_LIBRARY}")
    endif()
    set(GWR_LINEAR_ALGEBRA_LIBRARIES "${GWR_BLAS_LAPACK_LIBRARY}")
    message(STATUS
        "Using explicit BLAS/LAPACK library: ${GWR_BLAS_LAPACK_LIBRARY}")
else()
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
    set(GWR_LINEAR_ALGEBRA_LIBRARIES
        ${LAPACK_LIBRARIES}
        ${BLAS_LIBRARIES}
    )
endif()

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    add_compile_options(
        $<$<COMPILE_LANGUAGE:Fortran>:-std=f95>
        $<$<COMPILE_LANGUAGE:Fortran>:-Wall>
        $<$<COMPILE_LANGUAGE:Fortran>:-Wextra>
        $<$<COMPILE_LANGUAGE:Fortran>:-Wimplicit-interface>
        $<$<COMPILE_LANGUAGE:Fortran>:-Wconversion-extra>
    )
    # A non-default integer kind gives SYSTEM_CLOCK nanosecond resolution in
    # gfortran.  This Fortran-2003 timing extension is isolated here; all
    # numerical reference sources remain checked with strict Fortran 95.
    set_source_files_properties(
        src/performance.f90
        PROPERTIES COMPILE_OPTIONS -std=gnu
    )
endif()

add_library(gwr_parameteridentifikation
    src/precision.f90
    src/model.f90
    src/linear_algebra.f90
    src/evaluation.f90
)
target_link_libraries(gwr_parameteridentifikation PUBLIC
    ${GWR_LINEAR_ALGEBRA_LIBRARIES}
)

# Keep the high-resolution timing extension outside the strict-F95 numerical
# library.  Only benchmark executables need this support module.
add_library(gwr_benchmark_support src/performance.f90)
target_link_libraries(gwr_benchmark_support PUBLIC gwr_parameteridentifikation)

add_executable(gwr_reference src/run_reference.f90)
target_link_libraries(gwr_reference PRIVATE gwr_parameteridentifikation)

add_executable(gwr_benchmark src/run_benchmark.f90)
target_link_libraries(gwr_benchmark PRIVATE gwr_benchmark_support)

function(add_gwr_thread_benchmark thread_count cpu_mask)
    set(GWR_THREAD_COUNT "${thread_count}")
    set(GWR_CPU_MASK "${cpu_mask}")
    set(generated_source
        "${CMAKE_CURRENT_BINARY_DIR}/run_thread_benchmark_${thread_count}.f90")
    configure_file(
        src/run_thread_benchmark.f90.in
        "${generated_source}"
        @ONLY
    )
    add_executable(gwr_thread_benchmark_${thread_count} "${generated_source}")
    target_link_libraries(
        gwr_thread_benchmark_${thread_count}
        PRIVATE gwr_benchmark_support
    )
endfunction()

add_gwr_thread_benchmark(1 "0")
add_gwr_thread_benchmark(2 "0-1")
add_gwr_thread_benchmark(4 "0-3")
add_gwr_thread_benchmark(8 "0-7")
add_gwr_thread_benchmark(16 "0-15")
add_gwr_thread_benchmark(20 "0-19")
add_gwr_thread_benchmark(40 "0-39")

include(CTest)
if(BUILD_TESTING)
    add_executable(gwr_fortran_tests tests/test_reference.f90)
    target_link_libraries(gwr_fortran_tests PRIVATE gwr_parameteridentifikation)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        target_compile_options(gwr_fortran_tests PRIVATE -fcheck=all)
    endif()

    add_test(NAME fortran_reference_tests COMMAND gwr_fortran_tests)

    find_package(Python3 COMPONENTS Interpreter)
    if(Python3_Interpreter_FOUND)
        add_test(
            NAME compare_with_python
            COMMAND
                ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/tools/compare_with_python.py
                --fortran-executable $<TARGET_FILE:gwr_reference>
                --python-source
                ${CMAKE_CURRENT_SOURCE_DIR}/../python/src
                --output
                ${CMAKE_CURRENT_BINARY_DIR}/reference_fortran.csv
        )
    endif()
endif()

