
set(model_source_files Model.cpp BarotropicDynamics.cpp 
    fdm/ConformalProjectionFiniteDifferenceMethod.cpp
    fdm/AGridBarotropicDynamics.cpp)
add_library(model ${model_source_files})
target_link_libraries(model PUBLIC atlas eckit eccodes util)