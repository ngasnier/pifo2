#include "model/fdm/AGridBarotropicDynamics.h"
#include "model/fdm/ConformalProjectionFiniteDifferenceMethod.h"
#include "atlas/array/ArrayView.h"
#include <stdexcept>

namespace pifo {
    AGridBarotropicDynamics::AGridBarotropicDynamics(const atlas::numerics::Method& pMethod)
        : AGridBarotropicDynamics(pMethod, atlas::util::NoConfig())
    {

    }

    AGridBarotropicDynamics::AGridBarotropicDynamics(const atlas::numerics::Method& pMethod, const eckit::Parametrisation& pParam)
        : BarotropicDynamicsImpl(pMethod, pParam)
    {
        fdm = dynamic_cast<const pifo::ConformalProjectionFiniteDifferenceMethod*>( &pMethod );
        if (!fdm) 
        {
            throw std::runtime_error( "pifo::AGridBarotropicDynamics needs a pifo::ConformalProjectionFiniteDifferenceMethod "+ Here());
        }
        
        dx = fdm->getDx();
        dy = fdm->getDy();
    }

    AGridBarotropicDynamics::~AGridBarotropicDynamics() = default;


    void AGridBarotropicDynamics::calcU_tdcy(
        const atlas::Field& pV, 
        const atlas::Field& pphi, 
        const atlas::Field& pzeta,
        const atlas::Field& pK,
        atlas::Field& pU_tdcy) const
    {
        auto V = atlas::array::make_view<double, 1>(pV);
        auto phi = atlas::array::make_view<double, 1>(pphi);
        auto tourbillon = atlas::array::make_view<double, 1>(pzeta);
        auto K = atlas::array::make_view<double, 1>(pK);
        auto f = atlas::array::make_view<double, 1>(fdm->getF());
        auto U_tdcy = atlas::array::make_view<double, 1>(pU_tdcy);
        atlas::idx_t i = 0;
        atlas::idx_t stride = fdm->getFunctionSpace().grid().nx(0);
        atlas::idx_t nx = fdm->getFunctionSpace().grid().nx(0);
        atlas::idx_t ny = fdm->getFunctionSpace().grid().ny();
        #pragma omp parallel for private(i)
        for(atlas::idx_t y=1;y<ny-1;++y)
        {
            double c1, c2, c3, c4;
            double kphi=0;
            for(atlas::idx_t x=1;x<nx-1;++x)
            {
                i = x+y*stride;

                c1 = K[i+1];
                c2 = phi[i+1];
                c3 = K[i-1];
                c4 = phi[i-1];

                kphi = (c1+c2-c3-c4)/(2*dx);

                U_tdcy[i] = (tourbillon[i]+f[i])*V[i] - kphi;
            }
        }
    }
    
    void AGridBarotropicDynamics::calcV_tdcy(
        const atlas::Field& pU, 
        const atlas::Field& pphi, 
        const atlas::Field& pzeta,
        const atlas::Field& pK,
        atlas::Field& pV_tdcy) const
    {
        auto U = atlas::array::make_view<double, 1>(pU);
        auto phi = atlas::array::make_view<double, 1>(pphi);
        auto tourbillon = atlas::array::make_view<double, 1>(pzeta);
        auto K = atlas::array::make_view<double, 1>(pK);
        auto f = atlas::array::make_view<double, 1>(fdm->getF());
        auto V_tdcy = atlas::array::make_view<double, 1>(pV_tdcy);
        atlas::idx_t i = 0;
        atlas::idx_t stride = fdm->getFunctionSpace().grid().nx(0);
        atlas::idx_t nx = fdm->getFunctionSpace().grid().nx(0);
        atlas::idx_t ny = fdm->getFunctionSpace().grid().ny();
        #pragma omp parallel for private(i)
        for(atlas::idx_t y=1;y<ny-1;++y)
        {
            double c1, c2, c3, c4;
            double kphi=0;
            for(atlas::idx_t x=1;x<nx-1;++x)
            {
                i = x+y*stride;

                c1 = K[i-stride];
                c2 = phi[i-stride];
                c3 = K[i+stride];
                c4 = phi[i+stride];

                kphi = (c1+c2-c3-c4)/(2*dy);

                V_tdcy[i] = -(tourbillon[i]+f[i])*U[i] - kphi;
            }
        }
    }

    void AGridBarotropicDynamics::calcphi_tdcy(
        const atlas::Field& pU,
        const atlas::Field& pV,
        const atlas::Field& pphi,
        atlas::Field& pphi_tdcy) const
    {
        auto U = atlas::array::make_view<double, 1>(pU);
        auto V = atlas::array::make_view<double, 1>(pV);
        auto phi = atlas::array::make_view<double, 1>(pphi);
        auto m = atlas::array::make_view<double, 1>(fdm->getM());
        auto phi_tdcy = atlas::array::make_view<double, 1>(pphi_tdcy);
        atlas::idx_t i;
        atlas::idx_t stride = fdm->getFunctionSpace().grid().nx(0);
        atlas::idx_t nx = fdm->getFunctionSpace().grid().nx(0);
        atlas::idx_t ny = fdm->getFunctionSpace().grid().ny();
        #pragma omp parallel for private(i)
        for(atlas::idx_t y=1;y<ny-1;++y)
        {
            double mm = 0;
            for(atlas::idx_t x=1;x<nx-1;++x)
            {
                i = x+y*stride;

                mm = m[i];

                phi_tdcy[i] = -(mm*mm)*(
                        (phi[i+1]*U[i+1] - phi[i-1]*U[i-1])/(dx*2)
                        +(phi[i-stride]*V[i-stride] - phi[i+stride]*V[i+stride])/(dy*2)
                    );
            }
        }
    }

    void AGridBarotropicDynamics::calcK(
        const atlas::Field& pU,
        const atlas::Field& pV,
        atlas::Field& pK) const
    {
        auto U = atlas::array::make_view<double, 1>(pU);
        auto V = atlas::array::make_view<double, 1>(pV);
        auto K = atlas::array::make_view<double, 1>(pK);
        auto m = atlas::array::make_view<double, 1>(fdm->getM());
        atlas::idx_t i = 0;
        atlas::idx_t stride = fdm->getFunctionSpace().grid().nx(0);
        atlas::idx_t nx = fdm->getFunctionSpace().grid().nx(0);
        atlas::idx_t ny = fdm->getFunctionSpace().grid().ny();
        #pragma omp parallel for private(i)
        for(atlas::idx_t y=1;y<ny-1;++y)
        {
            double u1 = 0;
            double v1 = 0;
            for(atlas::idx_t x=1;x<nx-1;++x)
            {
                i = x+y*stride;
                u1 = U[i];
                v1 = V[i];
                K[i] = m[i]*m[i]*0.5*(u1*u1+v1*v1);
            }
        } 
    }

    void AGridBarotropicDynamics::calcZeta(
        const atlas::Field& pU,
        const atlas::Field& pV,
        atlas::Field& pzeta) const
    {
        auto U = atlas::array::make_view<double, 1>(pU);
        auto V = atlas::array::make_view<double, 1>(pV);
        auto tourbillon = atlas::array::make_view<double, 1>(pzeta);
        auto m = atlas::array::make_view<double, 1>(fdm->getM());
        atlas::idx_t i = 0;
        atlas::idx_t stride = fdm->getFunctionSpace().grid().nx(0);
        atlas::idx_t nx = fdm->getFunctionSpace().grid().nx(0);
        atlas::idx_t ny = fdm->getFunctionSpace().grid().ny();
        #pragma omp parallel for private(i)
        for (atlas::idx_t y=1;y<ny-1;++y)
        {
            double m1 = 0;
            double u1 = 0, u2 = 0;
            double v1 = 0, v2 = 0;
            for(atlas::idx_t x=1;x<nx-1;++x)
            {
                i = x+y*stride;

                m1 = m[i];

                u1 = U[i-stride];
                u2 = U[i+stride];

                v1 = V[i+1];
                v2 = V[i-1];

                tourbillon[i] = m1*m1
                        *((v1-v2)/(2*dx)
                        - (u1-u2)/(2*dy)
                        );
            }
        }
    }    
}