project(test)

set(LIBNM "libnm")

# integral
set(TEST "${PROJECT_NAME}-integral")
add_executable(${TEST} "test-integral.cc")
target_link_libraries(${TEST} ${LIBNM})
add_dependencies(${TEST} ${LIBNM})
