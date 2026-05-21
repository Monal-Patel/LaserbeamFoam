/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
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
    Foam::thincAdvection

Description
    THINC/QQ (Tangent of Hyperbola INterface Capturing with Quadratic
    surface representation and Quick-switching) volume fraction advection
    for structured grids.

    Reference:
        Xie, B., & Xiao, F. (2017).
        Toward efficient and accurate interface capturing on arbitrary hybrid
        unstructured grids: The THINC method with quadratic surface
        representation and Gaussian quadrature.
        Journal of Computational Physics, 349, 415-440.

SourceFiles
    thincAdvection.C

\*---------------------------------------------------------------------------*/

#ifndef thincAdvection_H
#define thincAdvection_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class thincAdvection Declaration
\*---------------------------------------------------------------------------*/

class thincAdvection
{
    // Private Data

        const fvMesh& mesh_;

        volScalarField& alpha1_;

        const surfaceScalarField& phi_;

        const volVectorField& U_;

        scalar thincBeta_;

        label nGaussPoints_;

        scalar vofTol_;

        surfaceScalarField alphaPhi_;

        vectorField cellDims_;


    // Private Member Functions

        void computeCellDims();

        vector cellLocalCoord
        (
            const point& pt,
            label celli
        ) const;

        scalar evaluateTHINC
        (
            scalar beta,
            const vector& normal,
            scalar d,
            const vector& xLocal
        ) const;

        scalar computeIntercept
        (
            scalar alpha,
            const vector& normal
        ) const;

        scalar computeInterceptGauss
        (
            scalar alpha,
            const vector& normal
        ) const;

        void computeAlphaPhi();

        //- No copy construct
        thincAdvection(const thincAdvection&) = delete;

        //- No copy assignment
        void operator=(const thincAdvection&) = delete;


public:

    // Constructors

        thincAdvection
        (
            volScalarField& alpha1,
            const surfaceScalarField& phi,
            const volVectorField& U,
            const dictionary& dict
        );


    //- Destructor
    ~thincAdvection() = default;


    // Member Functions

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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
