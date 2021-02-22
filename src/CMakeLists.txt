
add_subdirectory(app)
add_subdirectory(model)
add_subdirectory(util)

set(pifo2_source_files main.cpp)
add_executable(pifo2 ${pifo2_source_files})

target_link_libraries( pifo2 PUBLIC atlas eckit eccodes app model util)
