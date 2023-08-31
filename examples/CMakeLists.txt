add_executable(collinear collinear.f90)
target_link_libraries(collinear PUBLIC fasten blas lapack)
target_include_directories(collinear PUBLIC $<TARGET_PROPERTY:fasten,Fortran_MODULE_DIRECTORY>)

