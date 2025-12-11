cmake_minimum_required(VERSION 3.20)

project(Simulation LANGUAGES CXX)

find_package(Geant4 REQUIRED ui_all vis_all)

# --- Configuración TetGen (Windows) ---
set(TETGEN_DIR "C:/tetgen")                  
set(TETGEN_LIB "${TETGEN_DIR}/build/Release/tet.lib") 

file(GLOB SOURCES ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB HEADERS ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(sim main.cc ${SOURCES} ${HEADERS})

target_include_directories(sim PRIVATE
    ${PROJECT_SOURCE_DIR}/include   # Aquí debe estar CADMesh.hh
    ${TETGEN_DIR}                   # Aquí debe estar tetgen.h
)

target_link_libraries(sim PRIVATE
    ${Geant4_LIBRARIES}
    ${TETGEN_LIB}
)

# --- CORRECCIÓN CRÍTICA ---
# Habilita los bloques #ifdef CADMESH_USE_TETGEN en CADMesh.hh
target_compile_definitions(sim PRIVATE CADMESH_USE_TETGEN) 

set_target_properties(sim PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)