project(libnm)

file(GLOB SRCS "*.h" "*.hh" "*.hpp" "*.c" "*.cc" "*.cpp")
add_library(${PROJECT_NAME} ${SRCS})

target_include_directories(${PROJECT_NAME} PUBLIC .)
