project(test)

set(LIBNM "libnm")

# root
set(TEST "${PROJECT_NAME}-root")
add_executable(${TEST} "test-root.cc")
target_link_libraries(${TEST} ${LIBNM})
add_dependencies(${TEST} ${LIBNM})

# integral
set(TEST "${PROJECT_NAME}-integral")
add_executable(${TEST} "test-integral.cc")
target_link_libraries(${TEST} ${LIBNM})
add_dependencies(${TEST} ${LIBNM})

# abscissa
set(TEST "${PROJECT_NAME}-abscissa")
add_executable(${TEST} "test-abscissa.cc")
target_link_libraries(${TEST} ${LIBNM})
add_dependencies(${TEST} ${LIBNM})
