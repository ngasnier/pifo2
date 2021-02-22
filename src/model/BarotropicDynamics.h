#pragma once

#include "eckit/config/Parametrisation.h"
#include "atlas/numerics/Method.h"
#include "atlas/field.h"

namespace pifo
{
    class BarotropicDynamicsImpl
    {
    public:
        BarotropicDynamicsImpl(const atlas::numerics::Method&, const eckit::Parametrisation&);
        virtual ~BarotropicDynamicsImpl();

        virtual void calcU_tdcy(
            const atlas::Field& V, 
            const atlas::Field& phi, 
            const atlas::Field& zeta,
            const atlas::Field& K,
            atlas::Field& U_tdcy) const = 0;

        virtual void calcV_tdcy(
            const atlas::Field& U, 
            const atlas::Field& phi, 
            const atlas::Field& zeta,
            const atlas::Field& K,
            atlas::Field& V_tdcy) const = 0;

        virtual void calcphi_tdcy(
            const atlas::Field& U,
            const atlas::Field& V,
            const atlas::Field& phi,
            atlas::Field& phi_tdcy) const = 0;

        virtual void calcK(
            const atlas::Field& U,
            const atlas::Field& V,
            atlas::Field& K) const = 0;

        virtual void calcZeta(
            const atlas::Field& U,
            const atlas::Field& V,
            atlas::Field& zeta) const = 0;

    private:
        const atlas::numerics::Method& method;
    };

    class BarotropicDynamics
    {
    public:
        BarotropicDynamics() = default;

        void calcU_tdcy(
            const atlas::Field& V, 
            const atlas::Field& phi, 
            const atlas::Field& zeta,
            const atlas::Field& K,
            atlas::Field& U_tdcy) const;

        void calcV_tdcy(
            const atlas::Field& U, 
            const atlas::Field& phi, 
            const atlas::Field& zeta,
            const atlas::Field& K,
            atlas::Field& V_tdcy) const;

        void calcphi_tdcy(
            const atlas::Field& U,
            const atlas::Field& V,
            const atlas::Field& phi,
            atlas::Field& phi_tdcy) const;


        void calcK(
            const atlas::Field& U,
            const atlas::Field& V,
            atlas::Field& K) const;

        void calcZeta(
            const atlas::Field& U,
            const atlas::Field& V,
            atlas::Field& zeta) const;
    private:
        BarotropicDynamicsImpl* impl = nullptr;
    };
} // namespace pifo