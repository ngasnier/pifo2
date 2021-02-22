#include "WGribFormat.h"
#include "atlas/array/ArrayShape.h"
#include "atlas/array/ArrayView.h"
#include <stdexcept>
#include <fstream>

namespace pifo {

    void WGribFormat::readField(const std::string& file, const atlas::RegularGrid& grid, atlas::Field& field)
    {
        auto field_view = atlas::array::make_view<double, 1>(field);
        std::ifstream infile;
        atlas::idx_t nx, ny;
        double val;
        atlas::idx_t k=0;
        infile.open(file);
        infile >> nx;
        infile >> ny;
        if (nx!=grid.nx() || ny!=grid.ny()) throw std::runtime_error("file has not the same dimension as the grid to load. "+std::to_string(nx)+" "+std::to_string(ny));
        for (atlas::idx_t j=0;j<grid.ny();j++)
        {
            for (atlas::idx_t i=0;i<grid.nx();i++)
            {
                if (!infile.eof())
                {
                    infile >> val;
                    field_view(k) = val;
                }
                k++;
            }
        }
        infile.close();        
    }

    void WGribFormat::writeField(const std::string& file, const atlas::RegularGrid& grid, const atlas::Field& field)
    {
        auto field_view = atlas::array::make_view<double, 1>(field);
        std::ofstream outfile;
        atlas::idx_t k=0;
        outfile.open(file, std::ofstream::trunc);
        outfile << grid.nx() << " " << grid.ny() << std::endl;
        for (atlas::idx_t j=0;j<grid.ny();j++)
        {
            for (atlas::idx_t i=0;i<grid.nx();i++)
            {
                outfile << field_view(k) << std::endl;
                k++;
            }
        }
        outfile.close();        
    }

    void WGribFormat::writeLonLat(const std::string& lonfile, const std::string& latfile, const atlas::RegularGrid& grid)
    {
        std::ofstream outlonfile;
        std::ofstream outlatfile;
        outlonfile.open(lonfile, std::ofstream::trunc);
        outlatfile.open(latfile, std::ofstream::trunc);
        outlonfile << grid.nx() << " " << grid.ny() << std::endl;
        outlatfile << grid.nx() << " " << grid.ny() << std::endl;
        for (atlas::idx_t j=0;j<grid.ny();j++)
        {
            for (atlas::idx_t i=0;i<grid.nx();i++)
            {
                atlas::PointLonLat lonlat = grid.lonlat(i, j);
                outlonfile << lonlat.lon() << std::endl;
                outlatfile << lonlat.lat() << std::endl;
            }
        }
        outlonfile.close();
        outlatfile.close();
    }
}
