#include "ConformalProjectionFiniteDifferenceMethod.h"

namespace pifo {
    ConformalProjectionFiniteDifferenceMethod::ConformalProjectionFiniteDifferenceMethod(
        atlas::functionspace::StructuredColumns& fs, 
        atlas::Field& fField, 
        atlas::Field& mField) : atlas::numerics::Method(),
            functionSpace(fs),
            f(fField),
            m(mField)
    {
        if (functionSpace.grid().nx(0)>0)
            dx = functionSpace.grid().x(1, 0) - functionSpace.grid().x(0, 0);
        else
            dx = 0.0;
        if (dx<0) dx = -dx;

        if (functionSpace.grid().ny()>0)
            dy = functionSpace.grid().y(1) - functionSpace.grid().y(0);
        else
            dy = 0;
        if (dy<0) dy = -dy;
    }

    ConformalProjectionFiniteDifferenceMethod::~ConformalProjectionFiniteDifferenceMethod()
    {

    }

    const std::string& ConformalProjectionFiniteDifferenceMethod::name() const
    {
        return "ConformalProjectionFiniteDifferenceMethod";
    }
}