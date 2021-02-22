#include "model/BarotropicDynamics.h"

namespace pifo {
    BarotropicDynamicsImpl::BarotropicDynamicsImpl(const atlas::numerics::Method& pMethod, const eckit::Parametrisation& pParam)
    :method(pMethod)
    {
    }
        
    
    BarotropicDynamicsImpl::~BarotropicDynamicsImpl() = default;


    void BarotropicDynamics::calcU_tdcy(
        const atlas::Field& V, 
        const atlas::Field& phi, 
        const atlas::Field& zeta,
        const atlas::Field& K,
        atlas::Field& U_tdcy) const
    {
        impl->calcU_tdcy(V, phi, zeta, K, U_tdcy);
    }
    
    void BarotropicDynamics::calcV_tdcy(
        const atlas::Field& U, 
        const atlas::Field& phi, 
        const atlas::Field& zeta,
        const atlas::Field& K,
        atlas::Field& V_tdcy) const
    {
        impl->calcV_tdcy(U, phi, zeta, K, V_tdcy);
    }

    void BarotropicDynamics::calcphi_tdcy(
        const atlas::Field& U,
        const atlas::Field& V,
        const atlas::Field& phi,
        atlas::Field& phi_tdcy) const
    {
        impl->calcphi_tdcy(U, V, phi, phi_tdcy);
    }

    void BarotropicDynamics::calcK(
        const atlas::Field& U,
        const atlas::Field& V,
        atlas::Field& K) const
    {
        impl->calcK(U, V, K);
    }

    void BarotropicDynamics::calcZeta(
        const atlas::Field& U,
        const atlas::Field& V,
        atlas::Field& zeta) const
    {
        impl->calcZeta(U, V, zeta);
    }
}