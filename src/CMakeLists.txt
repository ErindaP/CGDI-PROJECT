add_executable(testPolyscope testPolyscope.cpp)
target_link_libraries(testPolyscope polyscope geometry-central)
add_executable(testSphere testSphere.cpp)
target_link_libraries(testSphere polyscope geometry-central)
