add_library(solver STATIC solver.cpp)
target_include_directories(solver PUBLIC "${PROJECT_SOURCE_DIR}/include")
