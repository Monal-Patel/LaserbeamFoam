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

\*---------------------------------------------------------------------------*/

#include "musclThincBvdAdvection.H"


void Foam::musclThincBvdAdvection::buildStencilMap()
{
    const label nIF = mesh_.nInternalFaces();
    stencilOm1_.setSize(nIF);
    stencilOm2_.setSize(nIF);
    stencilNp1_.setSize(nIF);
    stencilNp2_.setSize(nIF);

    const cellList& cells = mesh_.cells();
    const faceList& faces = mesh_.faces();
    const labelUList& own = mesh_.faceOwner();
    const labelUList& nei = mesh_.faceNeighbour();

    for (label fi = 0; fi < nIF; ++fi)
    {
        const label O = own[fi];
        const label N = nei[fi];

        label Om1 = O;
        label Om2 = O;
        const label fOppO = cells[O].opposingFaceLabel(fi, faces);
        if (fOppO >= 0 && fOppO < nIF)
        {
            Om1 = (own[fOppO] == O) ? nei[fOppO] : own[fOppO];
            const label fOppOm1 =
                cells[Om1].opposingFaceLabel(fOppO, faces);
            if (fOppOm1 >= 0 && fOppOm1 < nIF)
            {
                Om2 = (own[fOppOm1] == Om1)
                    ? nei[fOppOm1] : own[fOppOm1];
            }
            else
            {
                Om2 = Om1;
            }
        }

        label Np1 = N;
        label Np2 = N;
        const label fOppN = cells[N].opposingFaceLabel(fi, faces);
        if (fOppN >= 0 && fOppN < nIF)
        {
            Np1 = (own[fOppN] == N) ? nei[fOppN] : own[fOppN];
            const label fOppNp1 =
                cells[Np1].opposingFaceLabel(fOppN, faces);
            if (fOppNp1 >= 0 && fOppNp1 < nIF)
            {
                Np2 = (own[fOppNp1] == Np1)
                    ? nei[fOppNp1] : own[fOppNp1];
            }
            else
            {
                Np2 = Np1;
            }
        }

        stencilOm1_[fi] = Om1;
        stencilOm2_[fi] = Om2;
        stencilNp1_[fi] = Np1;
        stencilNp2_[fi] = Np2;
    }
}


void Foam::musclThincBvdAdvection::computeAlphaPhi()
{
    alpha1_.correctBoundaryConditions();

    const scalarField& alphaIn = alpha1_.primitiveField();
    const scalarField& phiIn   = phi_.primitiveField();
    const labelUList& own      = mesh_.faceOwner();
    const labelUList& nei      = mesh_.faceNeighbour();

    const label nInternalFaces = mesh_.nInternalFaces();

    auto clampAlpha = [](const scalar a) noexcept
    {
        return Foam::max(Foam::min(a, scalar(1.0)), scalar(0.0));
    };

    // Phase 1: per-face MUSCL+THINC reconstruction via direct stencil
    for (label facei = 0; facei < nInternalFaces; ++facei)
    {
        const label owni = own[facei];
        const label neii = nei[facei];

        const scalar alphaO   = clampAlpha(alphaIn[owni]);
        const scalar alphaN   = clampAlpha(alphaIn[neii]);
        const scalar alphaOm1 = clampAlpha(alphaIn[stencilOm1_[facei]]);
        const scalar alphaNp1 = clampAlpha(alphaIn[stencilNp1_[facei]]);
        const scalar alphaOm2 = clampAlpha(alphaIn[stencilOm2_[facei]]);
        const scalar alphaNp2 = clampAlpha(alphaIn[stencilNp2_[facei]]);

        reconstruct1d(alphaOm1, alphaO, alphaN, beta_, eps_,
                      qimhM_O_[facei], qiphM_O_[facei],
                      qimhT_O_[facei], qiphT_O_[facei]);

        reconstruct1d(alphaO, alphaN, alphaNp1, beta_, eps_,
                      qimhM_N_[facei], qiphM_N_[facei],
                      qimhT_N_[facei], qiphT_N_[facei]);

        {
            scalar tmp0, tmp1, tmp2, tmp3;
            reconstruct1d(alphaOm2, alphaOm1, alphaO, beta_, eps_,
                          tmp0, tmp1, tmp2, tmp3);
            qimhM_L_[facei] = tmp0;
            qiphM_L_[facei] = tmp1;
            qimhT_L_[facei] = tmp2;
            qiphT_L_[facei] = tmp3;
        }
        {
            scalar tmp0, tmp1, tmp2, tmp3;
            reconstruct1d(alphaN, alphaNp1, alphaNp2, beta_, eps_,
                          tmp0, tmp1, tmp2, tmp3);
            qimhM_R_[facei] = tmp0;
            qiphM_R_[facei] = tmp1;
            qimhT_R_[facei] = tmp2;
            qiphT_R_[facei] = tmp3;
        }
    }

    // Phase 2: BVD selection + upwind flux
    scalarField& alphaPhiIn = alphaPhi_.primitiveFieldRef();

    for (label facei = 0; facei < nInternalFaces; ++facei)
    {
        const label owni = own[facei];
        const label neii = nei[facei];

        const scalar alphaO   = clampAlpha(alphaIn[owni]);
        const scalar alphaN   = clampAlpha(alphaIn[neii]);
        const scalar alphaOm1 = clampAlpha(alphaIn[stencilOm1_[facei]]);
        const scalar alphaNp1 = clampAlpha(alphaIn[stencilNp1_[facei]]);

        const bool useTHINC_O = bvdSelect1d
        (
            qimhM_O_[facei], qiphM_O_[facei],
            qimhT_O_[facei], qiphT_O_[facei],
            qiphM_L_[facei], qiphT_L_[facei],
            qimhM_N_[facei], qimhT_N_[facei],
            alphaOm1, alphaO, alphaN,
            eps_, delta_
        );

        const bool useTHINC_N = bvdSelect1d
        (
            qimhM_N_[facei], qiphM_N_[facei],
            qimhT_N_[facei], qiphT_N_[facei],
            qiphM_O_[facei], qiphT_O_[facei],
            qimhM_R_[facei], qimhT_R_[facei],
            alphaO, alphaN, alphaNp1,
            eps_, delta_
        );

        const scalar alphaFacePos =
            useTHINC_O ? qiphT_O_[facei] : qiphM_O_[facei];
        const scalar alphaFaceNeg =
            useTHINC_N ? qimhT_N_[facei] : qimhM_N_[facei];

        alphaPhiIn[facei] = phiIn[facei]
            * ((phiIn[facei] >= 0) ? alphaFacePos : alphaFaceNeg);
    }

    // Boundary faces: upwind
    forAll(alphaPhi_.boundaryField(), patchi)
    {
        fvsPatchScalarField& alphaPhip =
            alphaPhi_.boundaryFieldRef()[patchi];
        const fvsPatchScalarField& phip = phi_.boundaryField()[patchi];
        const fvPatchScalarField& alphap =
            alpha1_.boundaryField()[patchi];
        const labelUList& fc = mesh_.boundary()[patchi].faceCells();

        if (mesh_.boundary()[patchi].coupled())
        {
            const scalarField alphaNbr(alphap.patchNeighbourField());

            forAll(alphaPhip, fi)
            {
                const label owni = fc[fi];
                alphaPhip[fi] = (phip[fi] >= 0)
                    ? phip[fi]*alphaIn[owni]
                    : phip[fi]*alphaNbr[fi];
            }
        }
        else
        {
            forAll(alphaPhip, fi)
            {
                alphaPhip[fi] = phip[fi]*alphaIn[fc[fi]];
            }
        }
    }
}


Foam::musclThincBvdAdvection::musclThincBvdAdvection
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
    beta_(dict.getOrDefault<scalar>("beta", 1.6)),
    eps_(dict.getOrDefault<scalar>("eps", 1.0e-20)),
    delta_(dict.getOrDefault<scalar>("delta", 1.0e-4)),
    alphaPhi_
    (
        IOobject
        (
            "alphaPhiMUSCL_THINC_BVD",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(phi.dimensions(), Zero)
    )
{
    const label nIF = mesh_.nInternalFaces();

    qimhM_O_.setSize(nIF);
    qiphM_O_.setSize(nIF);
    qimhT_O_.setSize(nIF);
    qiphT_O_.setSize(nIF);

    qimhM_N_.setSize(nIF);
    qiphM_N_.setSize(nIF);
    qimhT_N_.setSize(nIF);
    qiphT_N_.setSize(nIF);

    qimhM_L_.setSize(nIF);
    qiphM_L_.setSize(nIF);
    qimhT_L_.setSize(nIF);
    qiphT_L_.setSize(nIF);

    qimhM_R_.setSize(nIF);
    qiphM_R_.setSize(nIF);
    qimhT_R_.setSize(nIF);
    qiphT_R_.setSize(nIF);

    buildStencilMap();
}


void Foam::musclThincBvdAdvection::advect()
{
    computeAlphaPhi();

    const scalar dt = mesh_.time().deltaTValue();
    scalarField& alphaNew = alpha1_.primitiveFieldRef();
    const scalarField& alphaOld = alpha1_.oldTime().primitiveField();
    const auto& V = mesh_.V();
    const scalarField& alphaPhiIn = alphaPhi_.primitiveField();
    const scalarField& phiIn = phi_.primitiveField();
    const labelUList& own = mesh_.faceOwner();
    const labelUList& nei = mesh_.faceNeighbour();

    alphaNew = alphaOld;

    const label nIF = mesh_.nInternalFaces();
    for (label fi = 0; fi < nIF; ++fi)
    {
        const label o = own[fi];
        const label n = nei[fi];
        alphaNew[o] -= dt*(alphaPhiIn[fi] - alphaOld[o]*phiIn[fi]) / V[o];
        alphaNew[n] += dt*(alphaPhiIn[fi] - alphaOld[n]*phiIn[fi]) / V[n];
    }

    forAll(alphaPhi_.boundaryField(), patchi)
    {
        const auto& alphaPhip = alphaPhi_.boundaryField()[patchi];
        const auto& phip = phi_.boundaryField()[patchi];
        const labelUList& fc = mesh_.boundary()[patchi].faceCells();
        forAll(alphaPhip, fi)
        {
            const label c = fc[fi];
            alphaNew[c] -= dt*(alphaPhip[fi] - alphaOld[c]*phip[fi]) / V[c];
        }
    }

    alpha1_.primitiveFieldRef() =
        Foam::max
        (
            Foam::min(alpha1_.primitiveField(), scalar(1.0)),
            scalar(0.0)
        );

    alpha1_.correctBoundaryConditions();
}


Foam::tmp<Foam::surfaceScalarField> Foam::musclThincBvdAdvection::getRhoPhi
(
    const dimensionedScalar& rho1,
    const dimensionedScalar& rho2
) const
{
    return alphaPhi_*(rho1 - rho2) + phi_*rho2;
}


// ************************************************************************* //
