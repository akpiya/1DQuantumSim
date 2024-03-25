include(FetchContent)
FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG v3.4.0
)

FetchContent_MakeAvailable(Catch2)
list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)

enable_testing()
add_executable(unit_tests testing.cpp ${PROJECT_SOURCE_DIR}/include/solver.h)

target_link_libraries(unit_tests PRIVATE Catch2::Catch2WithMain solver)

target_include_directories(unit_tests PUBLIC
        "${PROJECT_SOURCE_DIR}/include"
        "${PROJECT_SOURCE_DIR}/src"
)