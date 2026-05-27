/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\      /   O peration     |
    \\    /    A nd         |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::musclThincBvdAdvection

Description
    Port-ready MUSCL-THINC-BVD advection class.

    The class now carries the AMReX-style limiter and reconstruction helpers
    needed for the port. The face-assembly step will be filled in next; the
    current advect() implementation still falls back to the existing THINC
    transport so the solver branch remains stable while the port is completed.

SourceFiles
    musclThincBvdAdvection.C

\*---------------------------------------------------------------------------*/

#ifndef musclThincBvdAdvection_H
#define musclThincBvdAdvection_H

#include "fvMesh.H"
#include "surfaceFields.H"
#include "volFields.H"

namespace Foam
{

inline scalar mcLimiter(const scalar dqR, const scalar dqL) noexcept
{
    const scalar a = scalar(2.0)*dqR;
    const scalar b = scalar(2.0)*dqL;
    const scalar c = scalar(0.5)*(dqR + dqL);

    if (a > 0 && b > 0 && c > 0)
    {
        return Foam::min(Foam::min(a, b), c);
    }

    if (a < 0 && b < 0 && c < 0)
    {
        return Foam::max(Foam::max(a, b), c);
    }

    return 0;
}


inline void reconstruct1d
(
    const scalar qim1,
    const scalar qi,
    const scalar qip1,
    const scalar beta,
    const scalar eps,
    scalar& qimhM,
    scalar& qiphM,
    scalar& qimhT,
    scalar& qiphT
) noexcept
{
    const scalar dq = mcLimiter(qip1 - qi, qi - qim1);
    qiphM = qi + scalar(0.5)*dq;
    qimhM = qi - scalar(0.5)*dq;

    const scalar qmin = Foam::min(qim1, qip1);
    const scalar qmax = Foam::max(qim1, qip1) - qmin;
    const scalar tB = std::tanh(beta);

    if (qmax < eps || Foam::mag(tB) < SMALL)
    {
        qiphT = qi;
        qimhT = qi;
        return;
    }

    const scalar theta = (qip1 > qim1) ? scalar(1.0) : scalar(-1.0);
    const scalar Craw = (qi - qmin + eps)/(qmax + eps);
    const scalar C = Foam::max(Foam::min(Craw, scalar(1.0)), scalar(0.0));
    const scalar arg = Foam::max
    (
        Foam::min(beta*theta*(scalar(2.0)*C - scalar(1.0)), scalar(50.0)),
        scalar(-50.0)
    );
    const scalar B = std::exp(arg);
    const scalar A = (B/std::cosh(beta) - scalar(1.0))/tB;

    qiphT = qmin + scalar(0.5)*qmax
         * (scalar(1.0) + theta*(tB + A)/(scalar(1.0) + A*tB));
    qimhT = qmin + scalar(0.5)*qmax*(scalar(1.0) + theta*A);
}


inline bool bvdSelect1d
(
    const scalar qimhM,
    const scalar qiphM,
    const scalar qimhT,
    const scalar qiphT,
    const scalar qiphMLeft,
    const scalar qiphTLeft,
    const scalar qimhMRight,
    const scalar qimhTRight,
    const scalar qim1,
    const scalar qi,
    const scalar qip1,
    const scalar eps,
    const scalar delta
) noexcept
{
    const scalar tbvT = Foam::min
    (
        Foam::min
        (
            Foam::mag(qiphMLeft - qimhT) + Foam::mag(qiphT - qimhMRight),
            Foam::mag(qiphTLeft - qimhT) + Foam::mag(qiphT - qimhTRight)
        ),
        Foam::min
        (
            Foam::mag(qiphMLeft - qimhT) + Foam::mag(qiphT - qimhTRight),
            Foam::mag(qiphTLeft - qimhT) + Foam::mag(qiphT - qimhMRight)
        )
    );

    const scalar tbvM = Foam::min
    (
        Foam::min
        (
            Foam::mag(qiphMLeft - qimhM) + Foam::mag(qiphM - qimhMRight),
            Foam::mag(qiphTLeft - qimhM) + Foam::mag(qiphM - qimhTRight)
        ),
        Foam::min
        (
            Foam::mag(qiphMLeft - qimhM) + Foam::mag(qiphM - qimhTRight),
            Foam::mag(qiphTLeft - qimhM) + Foam::mag(qiphM - qimhMRight)
        )
    );

    const scalar qmin = Foam::min(qim1, qip1);
    const scalar qmx = Foam::max(qim1, qip1) - qmin;
    const scalar C = (qi - qmin + eps)/(qmx + eps);

    return (C > delta)
        && (C < scalar(1.0) - delta)
        && ((qip1 - qi)*(qi - qim1) > scalar(0.0))
        && (tbvT < tbvM);
}

class musclThincBvdAdvection
{
    const fvMesh& mesh_;

    volScalarField& alpha1_;

    const surfaceScalarField& phi_;

    const volVectorField& U_;

    scalar beta_;

    scalar eps_;

    scalar delta_;

    surfaceScalarField alphaPhi_;

    scalarField qimhM_O_, qiphM_O_, qimhT_O_, qiphT_O_;
    scalarField qimhM_N_, qiphM_N_, qimhT_N_, qiphT_N_;
    scalarField qimhM_L_, qiphM_L_, qimhT_L_, qiphT_L_;
    scalarField qimhM_R_, qiphM_R_, qimhT_R_, qiphT_R_;

    labelField stencilOm1_;
    labelField stencilOm2_;
    labelField stencilNp1_;
    labelField stencilNp2_;

    void buildStencilMap();

    void computeAlphaPhi();

public:

    musclThincBvdAdvection
    (
        volScalarField& alpha1,
        const surfaceScalarField& phi,
        const volVectorField& U,
        const dictionary& dict
    );

    ~musclThincBvdAdvection() = default;

    void advect();

    tmp<surfaceScalarField> getRhoPhi
    (
        const dimensionedScalar& rho1,
        const dimensionedScalar& rho2
    ) const;

    const surfaceScalarField& alphaPhi() const noexcept
    {
        return alphaPhi_;
    }
};

} // End namespace Foam

#endif

// ************************************************************************* //