set(app_source_files Application.cpp DataProcessor.cpp ModelRun.cpp)
add_library(app ${app_source_files})
target_link_libraries(app PUBLIC atlas eckit eccodes model util)