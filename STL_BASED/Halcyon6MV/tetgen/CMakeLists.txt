cmake_minimum_required(VERSION 3.20)
project(TetGenLib LANGUAGES CXX)

# Archivos fuente
set(TETGEN_SOURCES
    tetgen.cxx
    predicates.cxx
)

# Crear librería estática
add_library(tet STATIC ${TETGEN_SOURCES})

# Definir que compila como librería, no ejecutable
target_compile_definitions(tet PRIVATE TETLIBRARY)

# Incluir headers
target_include_directories(tet PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
