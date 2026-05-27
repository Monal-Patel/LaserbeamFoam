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

\*---------------------------------------------------------------------------*/

#include "thincAdvection.H"
#include "fvc.H"
#include "fvcGrad.H"
#include "fvcDiv.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thincAdvection::thincAdvection
(
    volScalarField& alpha1,
    const surfaceScalarField& phi,
    const volVectorField& U,
    const dictionary& dict
)
:
    mesh_(alpha1.mesh()),
    alpha1_(alpha1),
    phi_(phi),
    U_(U),
    thincBeta_(dict.getOrDefault<scalar>("beta", 6.0)),
    nGaussPoints_(dict.getOrDefault<label>("nGaussPoints", 1)),
    vofTol_(dict.getOrDefault<scalar>("vofTol", 1e-3)),
    alphaPhi_
    (
        IOobject
        (
            "alphaPhiTHINC",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(phi.dimensions(), Zero)
    ),
    cellDims_(mesh_.nCells(), vector::zero)
{
    computeCellDims();

    Info<< "THINC/QQ advection initialised" << nl
        << "  beta = " << thincBeta_ << nl
        << "  nGaussPoints = " << nGaussPoints_ << nl
        << "  vofTol = " << vofTol_ << endl;
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::thincAdvection::computeCellDims()
{
    const cellList& cells = mesh_.cells();
    const vectorField& Sf = mesh_.faceAreas();
    const scalarField& V = mesh_.V();

    forAll(cells, celli)
    {
        const cell& c = cells[celli];

        scalar maxSfx = 0;
        scalar maxSfy = 0;
        scalar maxSfz = 0;

        forAll(c, fi)
        {
            const label facei = c[fi];
            const vector& sf = Sf[facei];

            maxSfx = Foam::max(maxSfx, mag(sf.x()));
            maxSfy = Foam::max(maxSfy, mag(sf.y()));
            maxSfz = Foam::max(maxSfz, mag(sf.z()));
        }

        cellDims_[celli] = vector
        (
            V[celli] / Foam::max(maxSfx, SMALL),
            V[celli] / Foam::max(maxSfy, SMALL),
            V[celli] / Foam::max(maxSfz, SMALL)
        );
    }
}


Foam::vector Foam::thincAdvection::cellLocalCoord
(
    const point& pt,
    label celli
) const
{
    const vector& cc = mesh_.C()[celli];
    const vector& dims = cellDims_[celli];

    return vector
    (
        (pt.x() - cc.x()) / dims.x() + 0.5,
        (pt.y() - cc.y()) / dims.y() + 0.5,
        (pt.z() - cc.z()) / dims.z() + 0.5
    );
}


Foam::scalar Foam::thincAdvection::evaluateTHINC
(
    scalar beta,
    const vector& normal,
    scalar d,
    const vector& xLocal
) const
{
    const scalar arg = beta*(normal & xLocal) + beta*d;

    if (arg > 20) return 1.0;
    if (arg < -20) return 0.0;

    return 0.5*(1.0 + Foam::tanh(arg));
}


Foam::scalar Foam::thincAdvection::computeIntercept
(
    scalar alpha,
    const vector& normal
) const
{
    // Midpoint rule: evaluate at cell centre (0.5, 0.5, 0.5)
    // H(0.5, 0.5, 0.5) = alpha
    // => d = atanh(2*alpha - 1)/beta - 0.5*(nx + ny + nz)

    const scalar clampedAlpha =
        Foam::max(Foam::min(alpha, 1.0 - SMALL), SMALL);

    return Foam::atanh(2.0*clampedAlpha - 1.0) / thincBeta_
         - 0.5*(normal.x() + normal.y() + normal.z());
}


Foam::scalar Foam::thincAdvection::computeInterceptGauss
(
    scalar alpha,
    const vector& normal
) const
{
    // 2-point Gaussian quadrature in 3D (2^3 = 8 points)
    static const scalar gp0 = 0.5 - Foam::sqrt(3.0)/6.0;
    static const scalar gp1 = 0.5 + Foam::sqrt(3.0)/6.0;
    static const scalar gw = 0.125; // (0.5)^3

    static const scalar gpArr[2] = {gp0, gp1};

    // Initial guess from midpoint rule
    scalar d = computeIntercept(alpha, normal);

    for (label iter = 0; iter < 8; ++iter)
    {
        scalar F = -alpha;
        scalar dFdd = 0;

        for (label i = 0; i < 2; ++i)
        for (label j = 0; j < 2; ++j)
        for (label k = 0; k < 2; ++k)
        {
            const scalar dotProd =
                normal.x()*gpArr[i]
              + normal.y()*gpArr[j]
              + normal.z()*gpArr[k];

            const scalar arg = thincBeta_*(dotProd + d);

            if (mag(arg) < 20)
            {
                const scalar th = Foam::tanh(arg);
                F += gw * 0.5*(1.0 + th);
                dFdd += gw * 0.5*thincBeta_*(1.0 - th*th);
            }
            else if (arg > 0)
            {
                F += gw;
            }
            // arg < -20: H = 0, contributes nothing
        }

        if (mag(dFdd) < SMALL) break;

        const scalar dd = -F / dFdd;
        d += dd;

        if (mag(dd) < 1e-12) break;
    }

    return d;
}


void Foam::thincAdvection::computeAlphaPhi()
{
    const scalarField& alphaIn = alpha1_.primitiveField();
    const vectorField& faceCentres = mesh_.Cf();
    const labelUList& own = mesh_.faceOwner();
    const labelUList& nei = mesh_.faceNeighbour();
    const scalarField& phiIn = phi_.primitiveField();

    scalarField& alphaPhiIn = alphaPhi_.primitiveFieldRef();

    alpha1_.correctBoundaryConditions();

    // Compute interface normals
    const volVectorField gradAlpha(fvc::grad(alpha1_));
    const vectorField& gradAlphaIn = gradAlpha.primitiveField();

    // Per-cell reconstruction data
    vectorField normals(mesh_.nCells(), vector::zero);
    scalarField intercepts(mesh_.nCells(), 0);
    boolList isInterface(mesh_.nCells(), false);

    forAll(alphaIn, celli)
    {
        if (alphaIn[celli] > vofTol_ && alphaIn[celli] < (1.0 - vofTol_))
        {
            const scalar magGrad = mag(gradAlphaIn[celli]);
            if (magGrad > SMALL)
            {
                isInterface[celli] = true;
                normals[celli] = gradAlphaIn[celli] / magGrad;

                if (nGaussPoints_ <= 1)
                {
                    intercepts[celli] =
                        computeIntercept(alphaIn[celli], normals[celli]);
                }
                else
                {
                    intercepts[celli] =
                        computeInterceptGauss(alphaIn[celli], normals[celli]);
                }
            }
        }
    }

    // BVD selection: evaluate THINC at faces and compare boundary variation
    const label nInternal = mesh_.nInternalFaces();
    scalarField thincFaceOwn(nInternal, 0);
    scalarField thincFaceNei(nInternal, 0);
    scalarField tbvTHINC(mesh_.nCells(), 0);
    scalarField tbvUpwind(mesh_.nCells(), 0);

    forAll(nei, facei)
    {
        const label o = own[facei];
        const label n = nei[facei];

        if (isInterface[o])
        {
            const vector xLocal = cellLocalCoord(faceCentres[facei], o);
            thincFaceOwn[facei] = evaluateTHINC
            (
                thincBeta_, normals[o], intercepts[o], xLocal
            );
            tbvTHINC[o] += mag(thincFaceOwn[facei] - alphaIn[n]);
            tbvUpwind[o] += mag(alphaIn[o] - alphaIn[n]);
        }

        if (isInterface[n])
        {
            const vector xLocal = cellLocalCoord(faceCentres[facei], n);
            thincFaceNei[facei] = evaluateTHINC
            (
                thincBeta_, normals[n], intercepts[n], xLocal
            );
            tbvTHINC[n] += mag(thincFaceNei[facei] - alphaIn[o]);
            tbvUpwind[n] += mag(alphaIn[n] - alphaIn[o]);
        }
    }

    boolList useTHINC(mesh_.nCells(), false);
    forAll(isInterface, celli)
    {
        if (isInterface[celli])
        {
            useTHINC[celli] = (tbvTHINC[celli] < tbvUpwind[celli]);
        }
    }

    // Internal face fluxes using BVD-selected reconstruction
    forAll(nei, facei)
    {
        scalar alphaFace;

        if (phiIn[facei] >= 0)
        {
            const label celli = own[facei];
            alphaFace = useTHINC[celli]
                ? thincFaceOwn[facei]
                : alphaIn[celli];
        }
        else
        {
            const label celli = nei[facei];
            alphaFace = useTHINC[celli]
                ? thincFaceNei[facei]
                : alphaIn[celli];
        }

        alphaPhiIn[facei] = phiIn[facei] * alphaFace;
    }

    // Boundary face fluxes (upwind)
    forAll(alphaPhi_.boundaryField(), patchi)
    {
        fvsPatchScalarField& alphaPhip =
            alphaPhi_.boundaryFieldRef()[patchi];
        const fvsPatchScalarField& phip = phi_.boundaryField()[patchi];
        const fvPatchScalarField& alphap = alpha1_.boundaryField()[patchi];
        const labelUList& fc = mesh_.boundary()[patchi].faceCells();

        if (mesh_.boundary()[patchi].coupled())
        {
            const scalarField alphaN(alphap.patchNeighbourField());

            forAll(alphaPhip, fi)
            {
                alphaPhip[fi] = (phip[fi] >= 0)
                    ? phip[fi] * alphaIn[fc[fi]]
                    : phip[fi] * alphaN[fi];
            }
        }
        else
        {
            forAll(alphaPhip, fi)
            {
                alphaPhip[fi] = phip[fi] * alphaIn[fc[fi]];
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::thincAdvection::advect()
{
    computeAlphaPhi();

    // Explicit Euler update with div(U) correction for numerical divergence
    alpha1_ = alpha1_.oldTime()
        - mesh_.time().deltaT()
        * (fvc::div(alphaPhi_) - alpha1_.oldTime() * fvc::div(phi_));

    // Bound alpha to [eps, 1-eps] to avoid division by zero elsewhere
    scalarField& alphaI = alpha1_.primitiveFieldRef();
    const scalar alphaMin = 1e-6;
    forAll(alphaI, i)
    {
        alphaI[i] = Foam::max(Foam::min(alphaI[i], 1.0 - alphaMin), alphaMin);
    }

    alpha1_.correctBoundaryConditions();
}


Foam::tmp<Foam::surfaceScalarField> Foam::thincAdvection::getRhoPhi
(
    const dimensionedScalar& rho1,
    const dimensionedScalar& rho2
) const
{
    return alphaPhi_*(rho1 - rho2) + phi_*rho2;
}


// ************************************************************************* //
