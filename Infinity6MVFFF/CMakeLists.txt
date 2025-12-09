cmake_minimum_required(VERSION 3.16...3.27)
project(Simulation LANGUAGES CXX) 

# Encuentra Geant4
find_package(Geant4 REQUIRED ui_all vis_all)

# Localiza fuentes y headers
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

# Crea el ejecutable con main en la raíz
add_executable(sim main.cc ${sources} ${headers})

# Incluye directorio de headers locales
target_include_directories(sim PRIVATE ${PROJECT_SOURCE_DIR}/include)

# Vincula las bibliotecas de Geant4
target_link_libraries(sim PRIVATE ${Geant4_LIBRARIES})

# Estándar C++
set_target_properties(sim PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

# (Opcional) Copiar archivos adicionales si los tienes
# set(SIM_SCRIPTS run.mac ...)
# foreach(_script ${SIM_SCRIPTS})
#   configure_file(${PROJECT_SOURCE_DIR}/${_script} ${PROJECT_BINARY_DIR}/${_script} COPYONLY)
# endforeach()
