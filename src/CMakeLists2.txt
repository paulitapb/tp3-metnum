cmake_minimum_required(VERSION 2.8.12)
project(tp2)

set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I ../")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I /usr/include/python3.8 ")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I /home/paula/Documents/COMPUTACION/Metodos-2c2022/tp2-metnum/src/pybind11")

IF(CMAKE_BUILD_TYPE MATCHES Debug)
    message("Debug mode")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g3 -O0")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")
ENDIF()

IF(CMAKE_BUILD_TYPE MATCHES Release)
    message("Release mode")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
ENDIF()

# si se usa un virtual environment setearlo aca
#set(PYTHON_EXECUTABLE "/path/to/bin/python")
#set(PYTHON_LIBRARY "/path/to/lib/libpytho.so")

# Prender y apagar este flag para soporte OpenMP
#set(USE_OpenMP ON)
# Script cmake ``multi-plataforma`` (testeado en MacOS y Ubuntu)
#include(cmake/OpenMP.cmake)

add_subdirectory(pybind11)

# Metnum seria el nombre del modulo a exportar
# Metnum.cpp contendria la definicion de los metodos a exportar
pybind11_add_module(metnum
                    ${PROJECT_SOURCE_DIR}/metnum.cpp
                    ${PROJECT_SOURCE_DIR}/utils.cpp
                    ${PROJECT_SOURCE_DIR}/eigen.cpp)

# Por cada módulo que use OpenMP ponerlo como dependencia de esta forma
#target_link_libraries(metnum LINK_PUBLIC OpenMP::OpenMP_CXX)

# Esta variable se usa para fijar el directorio de instalación
set(CMAKE_INSTALL_PREFIX
        ${PROJECT_SOURCE_DIR}
        CACHE PATH "Python custom modules directory" FORCE)

message(CMAKE_INSTALL_PREFIX = "${CMAKE_INSTALL_PREFIX}")

install(TARGETS metnum DESTINATION ${CMAKE_INSTALL_PREFIX}/notebooks)

# Crear un binario para testing manteniendo las dependencias de python por separado
add_executable(tp2
        ${PROJECT_SOURCE_DIR}/main.cpp
        ${PROJECT_SOURCE_DIR}/utils.cpp
        ${PROJECT_SOURCE_DIR}/eigen.cpp)


# si se quiere hacer un ejecutable que incluya pybind11, utilizar las
# siguientes 3 instrucciones de cmake.

# include_directories(pybind11/include)

# add_executable(tp2_pybind
#         src/main_pybind.cpp
#         src/eigen.cpp
#         src/utils.cpp)

# target_link_libraries(tp2_pybind PRIVATE pybind11::embed)
