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
#include "atlas/runtime/Trace.h"
#include "atlas/parallel/omp/omp.h"

#include <iostream>
#include <fstream>
#include <memory>
#include <cmath>

#include "../util/WGribFormat.h"

#include "ModelRun.h"
#include "../model/Model.h"

namespace pifo {
    namespace app  {
        void ModelRun::run()
        {
            atlas::Log::info() << "max threads : " << atlas_omp_get_max_threads() << std::endl;
            atlas_omp_set_num_threads(atlas_omp_get_max_threads()/2);

            atlas::Log::info () << "loading mercator grid" << std::endl ;
            atlas::util::Config mercator_config("regional_mercator.yml");
            atlas::RegularGrid mercator_grid(mercator_config);
            Model model(mercator_grid);
            
            auto fields = model.pronosticFieldSet().field_names();
            for (long unsigned int i=0;i<fields.size();i++)
            {
                atlas::Log::info () << "loading " << fields[i] << std::endl ;
                WGribFormat::readField(fields[i]+".txt", mercator_grid, model.pronosticFieldSet().field(fields[i]));
            }

            fields = model.parameterFieldSet().field_names();
            for (long unsigned int i=0;i<fields.size();i++)
            {
                atlas::Log::info () << "loading " << fields[i] << std::endl ;
                WGribFormat::readField(fields[i]+".txt", mercator_grid, model.parameterFieldSet().field(fields[i]));
            }

            atlas::Trace timer( Here(), "barotrope" );
            timer.start();
            double prevTime = 0.0;
            while (model.getTime()<3*3600)//  3*3600
            {
                model.step();
                timer.pause();
                atlas::Log::info () << "time " << model.getTime() << "s (elapsed : " << (timer.elapsed()-prevTime) << "s)" << std::endl ;
                prevTime = timer.elapsed();
                timer.resume();
            }
            timer.stop();
            atlas::Log::info () << "iteration finished. Total time : " << timer.elapsed() << "s)" << std::endl ;

            fields = model.pronosticFieldSet().field_names();
            for (long unsigned int i=0;i<fields.size();i++)
            {
                atlas::Log::info () << "writing " << fields[i] << std::endl ;
                WGribFormat::writeField(fields[i]+"_001.txt", mercator_grid, model.pronosticFieldSet().field(fields[i]));
            }
        }
    }
}