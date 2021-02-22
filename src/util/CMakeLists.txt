set(util_source_files GribFile.cpp Regridding.cpp WGribFormat.cpp)
add_library(util ${util_source_files})
target_link_libraries(util atlas eckit eccodes)