#pragma once

#include <memory>

#include "atlas/runtime/Log.h"
#include "atlas/grid.h"
#include "atlas/array/ArrayView.h"
#include "atlas/functionspace/StructuredColumns.h"
//#include "atlas/meshgenerator.h"
//#include "atlas/mesh.h"
#include "atlas/field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/parallel/omp/omp.h"
#include "model/fdm/ConformalProjectionFiniteDifferenceMethod.h"
#include "model/fdm/AGridBarotropicDynamics.h"

namespace pifo {
    class Model {
    public:
        Model(atlas::RegularGrid pgrid)
            : grid(atlas::RegularGrid(pgrid)),
            functionSpace(atlas::functionspace::StructuredColumns(pgrid)) // , atlas::option::halo(1)
        {
            pronosticFields = atlas::FieldSet("pronostics");
            pronosticFields.add(functionSpace.createField<double>(atlas::option::name("U")));
            pronosticFields.add(functionSpace.createField<double>(atlas::option::name("V")));
            pronosticFields.add(functionSpace.createField<double>(atlas::option::name("phi")));

            parameterFields = atlas::FieldSet("parameters");
            parameterFields.add(functionSpace.createField<double>(atlas::option::name("m")));
            parameterFields.add(functionSpace.createField<double>(atlas::option::name("f")));

            diagnosticFields = atlas::FieldSet("diagnostics");
            diagnosticFields.add(functionSpace.createField<double>(atlas::option::name("K")));
            diagnosticFields.add(functionSpace.createField<double>(atlas::option::name("zeta")));

            internalFields = atlas::FieldSet("internal");
            internalFields.add(functionSpace.createField<double>(atlas::option::name("U_t")));
            internalFields.add(functionSpace.createField<double>(atlas::option::name("V_t")));
            internalFields.add(functionSpace.createField<double>(atlas::option::name("phi_t")));
            internalFields.add(functionSpace.createField<double>(atlas::option::name("U_tdcy")));
            internalFields.add(functionSpace.createField<double>(atlas::option::name("V_tdcy")));
            internalFields.add(functionSpace.createField<double>(atlas::option::name("phi_tdcy")));
            
            method = std::unique_ptr<ConformalProjectionFiniteDifferenceMethod>(
                new ConformalProjectionFiniteDifferenceMethod(
                    functionSpace, 
                    parameterFields.field("f"), 
                    parameterFields.field("m")));
            dynamics = std::unique_ptr<AGridBarotropicDynamics>(new AGridBarotropicDynamics(*method));

            dt = 15;
            atlas::Log::info() << "init model for dx=" << method->getDx() 
                << " dy=" << method->getDy() 
                << " dt=" << dt 
                << " nx=" << grid.nx() 
                << " ny=" << grid.ny()
                << " halo=" << functionSpace.halo() 
                << " size=" << functionSpace.size() << "/" << pronosticFields.field("U").size()
                << " sizeHalo=" << functionSpace.sizeHalo()
                << std::endl;
            atlas::Log::info() << "iterate i : " 
                << "(" << functionSpace.i_begin(0) << "," << functionSpace.i_begin_halo(0) << "," << functionSpace.i_end_halo(0) << "," << functionSpace.i_end(0) << ")"
                << std::endl;
            atlas::Log::info() << "iterate j : " 
                << "(" << functionSpace.j_begin() << "," << functionSpace.j_begin_halo() << "," << functionSpace.j_end_halo() << "," << functionSpace.j_end() << ")"
                << std::endl;
            

            time = 0;
        }

        ~Model()
        {
        }

        const atlas::RegularGrid& getGrid()
        {
            return grid;        
        }

        atlas::FieldSet& pronosticFieldSet()
        {
            return pronosticFields;
        }

        atlas::FieldSet& parameterFieldSet()
        {
            return parameterFields;
        }

        atlas::FieldSet& diagnosticFieldSet()
        {
            return diagnosticFields;
        }

        atlas::FieldSet& internalFieldSet()
        {
            return internalFields;
        }

        double getDt()
        {
            return dt;
        }

        void setDt(double pdt)
        {
            dt = pdt;
        }

        double getTime()
        {
            return time;
        }

        void step()
        {
            calcK();
            calcZeta();
            calcU_tdcy();
            calcV_tdcy();
            calcphi_tdcy();

            if (time==0)
            {
                stepEuler();
            }
            else
            {
                stepLeapFrog();
            }

            finalizeStep();

            time += dt;
        }


    private:
        atlas::RegularGrid grid;
        atlas::functionspace::StructuredColumns functionSpace;
        atlas::FieldSet pronosticFields;
        atlas::FieldSet parameterFields;
        atlas::FieldSet diagnosticFields;
        atlas::FieldSet internalFields;
        
        std::unique_ptr<ConformalProjectionFiniteDifferenceMethod> method;
        std::unique_ptr<AGridBarotropicDynamics> dynamics;

        double dt;
        double time;

        void stepEuler()
        {    
            for (atlas::idx_t v=0;v<pronosticFields.size();v++)
            {
                auto vv = pronosticFields.field(v);
                auto var_t = internalFields.field(pronosticFields.field_names()[v]+"_t");
                auto var_tdcy = internalFields.field(pronosticFields.field_names()[v]+"_tdcy");
                a_bc(vv, var_tdcy, dt, var_t);
            }
        }

        void stepLeapFrog()
        {             
            for (atlas::idx_t v=0;v<pronosticFields.size();v++)
            {
                auto var_t = internalFields.field(pronosticFields.field_names()[v]+"_t");
                auto var_tdcy = internalFields.field(pronosticFields.field_names()[v]+"_tdcy");
                a_bc(var_t, var_tdcy, 2*dt, var_t);
            }
        }

        void a_bc(atlas::Field& a, atlas::Field& b, double c, atlas::Field& dest)
        {
#ifndef PIFO_FAST_MODE
            auto x = atlas::array::make_view<double, 1>(a);
            auto y = atlas::array::make_view<double, 1>(b);
            auto z = atlas::array::make_view<double, 1>(dest);
#else
            auto x = (double*)a.storage();
            auto y = (double*)b.storage();
            auto z = (double*)dest.storage();
#endif
            atlas::idx_t size = a.size();
            #pragma omp parallel for 
            for(atlas::idx_t i=0;i<size;i++)
            {
                z[i] = x[i]+c*y[i];
            }
        }

        void swap(atlas::Field& a, atlas::Field& b)
        {
#ifndef PIFO_FAST_MODE
            auto x = atlas::array::make_view<double, 1>(a);
            auto y = atlas::array::make_view<double, 1>(b);
#else
            auto x = (double*)a.storage();
            auto y = (double*)b.storage();
#endif
            atlas::idx_t size = a.size();
            #pragma omp parallel for 
            for(atlas::idx_t i=0;i<size;i++)
            {
                double tmp;
                tmp = x[i];
                x[i] = y[i];
                y[i] = tmp;
            }

        }

        void finalizeStep()
        {
            for (atlas::idx_t v=0;v<pronosticFields.size();++v)
            {
                swap(pronosticFields.field(pronosticFields.field_names()[v]), internalFields.field(pronosticFields.field_names()[v]+"_t"));
            }
        }

        void calcU_tdcy()
        {
            dynamics->calcU_tdcy(pronosticFields.field("V"),
                pronosticFields.field("phi"),
                diagnosticFields.field("zeta"),
                diagnosticFields.field("K"),
                internalFields.field("U_tdcy"));
        }

        void calcV_tdcy()
        {
            dynamics->calcV_tdcy(pronosticFields.field("U"),
                pronosticFields.field("phi"),
                diagnosticFields.field("zeta"),
                diagnosticFields.field("K"),
                internalFields.field("V_tdcy"));
        }

        void calcphi_tdcy()
        {
            dynamics->calcphi_tdcy(
                pronosticFields.field("U"),
                pronosticFields.field("V"),
                pronosticFields.field("phi"),
                internalFields.field("phi_tdcy"));
        }

        void calcK()
        {
            dynamics->calcK(pronosticFields.field("U"),
                pronosticFields.field("V"),
                diagnosticFields.field("K"));
        }

        void calcZeta()
        {
            dynamics->calcZeta(pronosticFields.field("U"),
                pronosticFields.field("V"),
                diagnosticFields.field("zeta"));
       }
    };
}