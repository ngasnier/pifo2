#pragma once

#include "model/BarotropicDynamics.h"
#include "ConformalProjectionFiniteDifferenceMethod.h"

namespace pifo
{
    class AGridBarotropicDynamics : public BarotropicDynamicsImpl
    {
    public:
        AGridBarotropicDynamics(const atlas::numerics::Method&);
        AGridBarotropicDynamics(const atlas::numerics::Method&, const eckit::Parametrisation&);
        virtual ~AGridBarotropicDynamics();

        virtual void calcU_tdcy(
            const atlas::Field& V, 
            const atlas::Field& phi, 
            const atlas::Field& zeta,
            const atlas::Field& K,
            atlas::Field& U_tdcy) const;

        virtual void calcV_tdcy(
            const atlas::Field& U, 
            const atlas::Field& phi, 
            const atlas::Field& zeta,
            const atlas::Field& K,
            atlas::Field& V_tdcy) const;

        virtual void calcphi_tdcy(
            const atlas::Field& U,
            const atlas::Field& V,
            const atlas::Field& phi,
            atlas::Field& phi_tdcy) const;

        virtual void calcK(
            const atlas::Field& U,
            const atlas::Field& V,
            atlas::Field& K) const;

        virtual void calcZeta(
            const atlas::Field& U,
            const atlas::Field& V,
            atlas::Field& zeta) const;
    private:
        ConformalProjectionFiniteDifferenceMethod const* fdm;
        double dx;
        double dy;
    };

}