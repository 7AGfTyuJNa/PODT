cmake_minimum_required(VERSION 3.28)
project(PODT-code VERSION 0.1.0)

include(CTest)
enable_testing()

include_directories(/usr/local/include)
link_directories(/usr/local/lib)

add_executable(PODT-code main.cpp ODTE-tester.cpp TDSC20.cpp ODTE.cpp HSS.cpp helper.cpp)
target_link_libraries(PODT-code ntl m gmp pthread relic)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)