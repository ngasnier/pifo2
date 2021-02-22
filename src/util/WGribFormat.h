#pragma once

#include <string>
#include "atlas/grid.h"
#include "atlas/field.h"

namespace pifo {
    class WGribFormat {
        public:
            static void readField(const std::string& file, const atlas::RegularGrid& grid, atlas::Field& field);
            static void writeField(const std::string& file, const atlas::RegularGrid& grid, const atlas::Field& field);
            static void writeLonLat(const std::string& lonfile, const std::string& latfile, const atlas::RegularGrid& grid);
    };
}