#pragma once

#include <string>

#include "atlas/numerics/Method.h"
#include "atlas/functionspace/StructuredColumns.h"

namespace pifo { 
    class ConformalProjectionFiniteDifferenceMethod : public atlas::numerics::Method {
    public:
        ConformalProjectionFiniteDifferenceMethod(atlas::functionspace::StructuredColumns& fs, atlas::Field& fField, atlas::Field& mField);
        
        ~ConformalProjectionFiniteDifferenceMethod();

        virtual const std::string& name() const;

        atlas::Field& getM() const
        {
            return m;
        }

        atlas::Field& getF() const
        {
            return f;
        }

        double getDx() const
        {
            return dx;
        }   

        double getDy() const
        {
            return dy;
        }     

        atlas::functionspace::StructuredColumns& getFunctionSpace() const
        {
            return functionSpace;
        }

    private:
        atlas::functionspace::StructuredColumns& functionSpace;
        atlas::Field& m;
        atlas::Field& f;
        double dx;
        double dy;
    };
}