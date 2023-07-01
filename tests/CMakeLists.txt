project(test)

# root
set(TEST "${PROJECT_NAME}-root")
add_executable(${TEST} "test-root.cc")
target_link_libraries(${TEST} libnm)
add_dependencies(${TEST} libnm)

# integral
set(TEST "${PROJECT_NAME}-integral")
add_executable(${TEST} "test-integral.cc")
target_link_libraries(${TEST} libnm)
add_dependencies(${TEST} libnm)

# abscissa
set(TEST "${PROJECT_NAME}-abscissa")
add_executable(${TEST} "test-abscissa.cc")
target_link_libraries(${TEST} libnm)
add_dependencies(${TEST} libnm)
