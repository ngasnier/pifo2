#include "atlas/library/Library.h"
#include "atlas/runtime/Log.h"

#include "atlas/util/Config.h"
#include "atlas/grid.h"
#include "atlas/grid/Distribution.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/meshgenerator.h"
#include "atlas/meshgenerator/MeshGenerator.h"
#include "atlas/mesh.h"
#include "atlas/field.h"
#include "atlas/array/DataType.h"
#include "atlas/array/ArrayShape.h"
#include "atlas/array/ArrayView.h"
#include "atlas/interpolation/Interpolation.h"
#include "atlas/interpolation/method/Method.h"
#include "atlas/interpolation/method/MethodFactory.h"

#include <iostream>
#include <fstream>
#include <cmath>

#include "../util/GribFile.h"
#include "../util/WGribFormat.h"
#include "../util/Regridding.h"

#include "DataProcessor.h"

namespace pifo {
    namespace app  {
        using namespace std;
        using namespace atlas;
        using namespace atlas::util;
 
        void dataToField(const double* data, const RegularGrid& grid, Field& field)
        {
            auto field_view = atlas::array::make_view<double, 1>(field);
            idx_t k = 0;
            for (idx_t j=0;j<grid.ny();j++)
            {
                for (idx_t i=0;i<grid.nx();i++)
                {
                    field_view(k) = data[k];
                    k++;
                }
            }
        }

        void multiplyConst(const RegularGrid& grid, Field& field, double value)
        {
            auto field_view = atlas::array::make_view<double, 1>(field);
            idx_t k = 0;
            for (idx_t j=0;j<grid.ny();j++)
            {
                for (idx_t i=0;i<grid.nx();i++)
                {
                    field_view(k) = field_view(k) * value;
                    k++;
                }
            }
        }

        void addConst(const RegularGrid& grid, Field& field, double value)
        {
            auto field_view = atlas::array::make_view<double, 1>(field);
            idx_t k = 0;
            for (idx_t j=0;j<grid.ny();j++)
            {
                for (idx_t i=0;i<grid.nx();i++)
                {
                    field_view(k) = field_view(k) + value;
                    k++;
                }
            }
        }

        void scaleField(const RegularGrid& grid, Field& field, const Field& m)
        {
            auto field_view = atlas::array::make_view<double, 1>(field);
            auto m_view = atlas::array::make_view<double, 1>(m);
            idx_t k = 0;
            for (idx_t j=0;j<grid.ny();j++)
            {
                for (idx_t i=0;i<grid.nx();i++)
                {
                    field_view(k) = field_view(k)/m_view(k);
                    k++;
                }
            }
        }

        void calcScalingFactor(const RegularGrid& grid, Field& field)
        {
            auto field_view = atlas::array::make_view<double, 1>(field);
            idx_t k = 0;
            for (idx_t j=0;j<grid.ny();j++)
            {
                for (idx_t i=0;i<grid.nx();i++)
                {
                    PointLonLat lonlat = grid.lonlat(i, j);
                    field_view(k) = 1/cos(lonlat.lat()*3.14159265/180.0);
                    k++;
                }
            }
        }
 
        void calcCoriolisFactor(const RegularGrid& grid, Field& field)
        {
            auto field_view = atlas::array::make_view<double, 1>(field);
            idx_t k = 0;
            for (idx_t j=0;j<grid.ny();j++)
            {
                for (idx_t i=0;i<grid.nx();i++)
                {
                    PointLonLat lonlat = grid.lonlat(i, j);
                    field_view(k) = 2 * 7.292115e-5 * sin(lonlat.lat()*3.14159265/180.0);
                    k++;
                }
            }
        }

        void DataProcessor::run()
        {
            atlas::Log::info () << "loading lonlat grid" << std::endl ;
            Config latlon_config("regular_lonlat.yml");
            RegularGrid latlon_grid(latlon_config);

            FieldSet latlon_fields("latlon_fields");
            auto latlon_prmsl = std::make_unique<Field>("prmsl", atlas::array::make_datatype<double>(), atlas::array::make_shape(latlon_grid.size()));
            auto latlon_z500 = std::make_unique<Field>("z500", atlas::array::make_datatype<double>(), atlas::array::make_shape(latlon_grid.size()));
            auto latlon_u = std::make_unique<Field>("u", atlas::array::make_datatype<double>(), atlas::array::make_shape(latlon_grid.size()));
            auto latlon_v = std::make_unique<Field>("v", atlas::array::make_datatype<double>(), atlas::array::make_shape(latlon_grid.size()));
            latlon_fields.add(*latlon_prmsl);
            latlon_fields.add(*latlon_z500);
            latlon_fields.add(*latlon_u);
            latlon_fields.add(*latlon_v);

            double* latlon_data = new double[latlon_grid.size()];
            double* latlon_lats = new double[latlon_grid.size()];
            double* latlon_lons = new double[latlon_grid.size()];

            auto latlon_pressure = atlas::array::make_view<double, 1>(*latlon_prmsl);
            idx_t nb_lats = 0;
            idx_t nb_lons = 0;
            double prev_lat = -1000000000;
            double prev_lon = -1000000000;
            idx_t k = 0;
            for (idx_t j=0;j<latlon_grid.ny();j++)
            {
                for (idx_t i=0;i<latlon_grid.nx();i++)
                {
                    latlon_pressure(k) = latlon_data[k];
                    k++;
                }
                PointLonLat lonlat = latlon_grid.lonlat(0, j);
                if (lonlat.lat()!=prev_lat) { prev_lat = latlon_lats[nb_lats] = lonlat.lat(); nb_lats++; }
            }
            for (idx_t i=0;i<latlon_grid.nx();i++)
            {
                PointLonLat lonlat = latlon_grid.lonlat(i, 0);
                if (lonlat.lon()!=prev_lon) { prev_lon = latlon_lons[nb_lons] = lonlat.lon(); nb_lons++; }
            }

            atlas::Log::info () << "loading mercator grid" << std::endl ;
            Config mercator_config("regional_mercator.yml");
            RegularGrid mercator_grid(mercator_config);
                       
            FieldSet mercator_fields("mercator_fields");
            auto mercator_prmsl = std::make_unique<Field>("prmsl", atlas::array::make_datatype<double>(), atlas::array::make_shape(mercator_grid.size()));
            auto mercator_u = std::make_unique<Field>("u", atlas::array::make_datatype<double>(), atlas::array::make_shape(mercator_grid.size()));
            auto mercator_v = std::make_unique<Field>("v", atlas::array::make_datatype<double>(), atlas::array::make_shape(mercator_grid.size()));
            auto mercator_phi = std::make_unique<Field>("phi", atlas::array::make_datatype<double>(), atlas::array::make_shape(mercator_grid.size()));
            auto mercator_f = std::make_unique<Field>("f", atlas::array::make_datatype<double>(), atlas::array::make_shape(mercator_grid.size()));
            auto mercator_m = std::make_unique<Field>("m", atlas::array::make_datatype<double>(), atlas::array::make_shape(mercator_grid.size()));
            mercator_fields.add(*mercator_prmsl);
            mercator_fields.add(*mercator_u);
            mercator_fields.add(*mercator_v);
            mercator_fields.add(*mercator_phi);

            double* mercator_data = new double[mercator_grid.size()];
            double* mercator_lats = new double[mercator_grid.size()];
            double* mercator_lons = new double[mercator_grid.size()];
            k=0;
            for (idx_t j=0;j<mercator_grid.ny();j++)
            {
                for (idx_t i=0;i<mercator_grid.nx();i++)
                {
                    PointLonLat lonlat = mercator_grid.lonlat(i, j);
                    mercator_lats[k] = lonlat.lat();
                    mercator_lons[k] = lonlat.lon();
                    k++;
                }
            }

            try  {
                GribFile grb("/home/nicolas/Meteo/Products/modeldata/global/gfs/2019112406/gfs.t06z.pgrb2.0p50.f000");

                atlas::Log::info () << "calculating m" << std::endl ;
                calcScalingFactor(mercator_grid, *mercator_m);
                WGribFormat::writeField("m.txt", mercator_grid, *mercator_m);

                atlas::Log::info () << "calculating f" << std::endl ;
                calcCoriolisFactor(mercator_grid, *mercator_f);
                WGribFormat::writeField("f.txt", mercator_grid, *mercator_f);

                atlas::Log::info () << "loading prmsl" << std::endl ;
                grb.getData("prmsl", "meanSea", 0, latlon_data);
                dataToField(latlon_data, latlon_grid, *latlon_prmsl);
                Log :: info () << "interpolating prmsl" << std::endl ;
                Regridding::bilinearRegrid(latlon_lons, nb_lons, latlon_lats, nb_lats, latlon_data, latlon_grid.size(), true, mercator_lons, mercator_lats, mercator_data, mercator_grid.size());
                dataToField(mercator_data, mercator_grid, *mercator_prmsl);
                WGribFormat::writeField("prmsl.txt", mercator_grid, *mercator_prmsl);

                atlas::Log::info () << "loading z500" << std::endl ;
                grb.getData("gh", "isobaricInhPa", 500, latlon_data);
                dataToField(latlon_data, latlon_grid, *latlon_z500);
                Log :: info () << "interpolating z500" << std::endl ;
                Regridding::bilinearRegrid(latlon_lons, nb_lons, latlon_lats, nb_lats, latlon_data, latlon_grid.size(), true, mercator_lons, mercator_lats, mercator_data, mercator_grid.size());
                dataToField(mercator_data, mercator_grid, *mercator_phi);
                multiplyConst(mercator_grid, *mercator_phi, 9.8066);
                addConst(mercator_grid, *mercator_phi, -40000);
                WGribFormat::writeField("phi.txt", mercator_grid, *mercator_phi);
               
                atlas::Log::info () << "loading u500" << std::endl ;
                grb.getData("u", "isobaricInhPa", 500, latlon_data);
                dataToField(latlon_data, latlon_grid, *latlon_u);
                Log :: info () << "interpolating u500" << std::endl ;
                Regridding::bilinearRegrid(latlon_lons, nb_lons, latlon_lats, nb_lats, latlon_data, latlon_grid.size(), true, mercator_lons, mercator_lats, mercator_data, mercator_grid.size());
                dataToField(mercator_data, mercator_grid, *mercator_u);
                scaleField(mercator_grid, *mercator_u, *mercator_m);
                WGribFormat::writeField("U.txt", mercator_grid, *mercator_u);

                atlas::Log::info () << "loading v500" << std::endl ;
                grb.getData("v", "isobaricInhPa", 500, latlon_data);
                dataToField(latlon_data, latlon_grid, *latlon_v);
                Log :: info () << "interpolating v500" << std::endl ;
                Regridding::bilinearRegrid(latlon_lons, nb_lons, latlon_lats, nb_lats, latlon_data, latlon_grid.size(), true, mercator_lons, mercator_lats, mercator_data, mercator_grid.size());
                dataToField(mercator_data, mercator_grid, *mercator_v);
                scaleField(mercator_grid, *mercator_v, *mercator_m);
                WGribFormat::writeField("V.txt", mercator_grid, *mercator_v);

                WGribFormat::writeLonLat("lons.txt", "lats.txt", mercator_grid);
            }
            catch (const std::exception& ex)
            {
                atlas::Log::error() << ex.what() << endl;
            }

            atlas::Log::info() << std::endl;

            delete latlon_data;
            delete latlon_lats;
            delete latlon_lons;

            delete mercator_data;
            delete mercator_lats;
            delete mercator_lons;            
        }
    }
}