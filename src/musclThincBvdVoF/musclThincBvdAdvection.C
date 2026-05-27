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

#include "fvc.H"


void Foam::musclThincBvdAdvection::computeAlphaPhi()
{
    alpha1_.correctBoundaryConditions();

    const volVectorField gradAlpha(fvc::grad(alpha1_));
    const vectorField& gradAlphaIn = gradAlpha.primitiveField();

    const scalarField& alphaIn = alpha1_.primitiveField();
    const scalarField& phiIn   = phi_.primitiveField();
    const labelUList& own      = mesh_.faceOwner();
    const labelUList& nei      = mesh_.faceNeighbour();
    const vectorField& cellCentres = mesh_.C();

    const label nInternalFaces = mesh_.nInternalFaces();

    auto clampAlpha = [](const scalar a) noexcept
    {
        return Foam::max(Foam::min(a, scalar(1.0)), scalar(0.0));
    };

    // -----------------------------------------------------------------------
    // Phase 1 (AMReX bxg2 pass): per-cell MUSCL+THINC reconstruction.
    //
    // For each cell i we need qim1 and qip1 along the face-normal direction.
    // In an unstructured mesh there is no single "i-1 / i+1" cell.  We use
    // the gradient-extrapolation approach: qim1 and qip1 are approximated by
    // extrapolating the cell-centred alpha along its own gradient to the
    // two ends of the cell, mirroring each end point about the cell centre
    // to get a virtual upstream / downstream value.  This is equivalent to
    // the AMReX structured approach when the grid is uniform and Cartesian.
    //
    // Storage: 4 scalars per cell [qimhM, qiphM, qimhT, qiphT].
    // -----------------------------------------------------------------------

    // We build one reconstruction per internal face direction, computed from
    // the perspective of each cell in its own face-normal direction.
    // reco[cellI] holds the reconstruction along the direction toward the
    // face that will be looked up.  Because each face picks a specific
    // owner/neighbour pair, we pre-compute per-cell arrays indexed by face
    // for the four relevant values.
    //
    // Concretely, for each internal face f with owner O and neighbour N:
    //   direction d = unit(C_N - C_O)
    //   For cell O along d: qim1 = alpha(O) - grad(O)·d·hO,  qip1 = alpha(N)
    //   For cell N along d: qim1 = alpha(O),  qip1 = alpha(N) + grad(N)·d·hN
    // where hO = |C_N - C_O|, hN = |C_N - C_O| (same span).
    //
    // Second-layer values for the BVD neighbourhood are obtained by the same
    // gradient extrapolation: the "upstream of O" value is
    //   alpha_Om2 = alpha(O) - 2 * grad(O)·d·hO
    // and "downstream of N":
    //   alpha_Np2 = alpha(N) + 2 * grad(N)·d·hN.
    //
    // This maps exactly to the AMReX Algorithm 2 neighbourhood structure.
    // -----------------------------------------------------------------------

    // Per-face reconstruction arrays (owner-perspective and neighbour-perspective).
    // [face] → {qimhM, qiphM, qimhT, qiphT}
    scalarField qimhM_O(nInternalFaces);
    scalarField qiphM_O(nInternalFaces);
    scalarField qimhT_O(nInternalFaces);
    scalarField qiphT_O(nInternalFaces);

    scalarField qimhM_N(nInternalFaces);
    scalarField qiphM_N(nInternalFaces);
    scalarField qimhT_N(nInternalFaces);
    scalarField qiphT_N(nInternalFaces);

    // Second-layer reconstructions needed for BVD neighbourhood
    scalarField qimhM_L(nInternalFaces);
    scalarField qiphM_L(nInternalFaces);
    scalarField qimhT_L(nInternalFaces);
    scalarField qiphT_L(nInternalFaces);

    scalarField qimhM_R(nInternalFaces);
    scalarField qiphM_R(nInternalFaces);
    scalarField qimhT_R(nInternalFaces);
    scalarField qiphT_R(nInternalFaces);

    for (label facei = 0; facei < nInternalFaces; ++facei)
    {
        const label owni = own[facei];
        const label neii = nei[facei];

        const scalar alphaO = clampAlpha(alphaIn[owni]);
        const scalar alphaN = clampAlpha(alphaIn[neii]);

        // Face-normal direction (O → N), scaled to full cell span
        const vector d  = cellCentres[neii] - cellCentres[owni];
        const scalar hON = Foam::mag(d);
        const vector dn = (hON > SMALL) ? d/hON : vector::zero;

        // Gradient projections onto the O→N direction (one cell span)
        const scalar slopeO = (gradAlphaIn[owni] & dn) * hON;
        const scalar slopeN = (gradAlphaIn[neii] & dn) * hON;

        // Virtual cell values (gradient extrapolation, AMReX i-1 / i+1 / i+2)
        const scalar alphaOm1 = clampAlpha(alphaO - slopeO);  // left of O (i-1)
        const scalar alphaNp1 = clampAlpha(alphaN + slopeN);  // right of N (i+2 relative to O)

        // Double-step for second-layer cells (one more span each way)
        const scalar alphaOm2 = clampAlpha(alphaO - scalar(2.0)*slopeO);
        const scalar alphaNp2 = clampAlpha(alphaN + scalar(2.0)*slopeN);

        // --- Phase 1: reconstruct_1d for each relevant cell in the stencil ---

        // Owner cell O (stencil: alphaOm1, alphaO, alphaN)
        reconstruct1d(alphaOm1, alphaO, alphaN, beta_, eps_,
                      qimhM_O[facei], qiphM_O[facei],
                      qimhT_O[facei], qiphT_O[facei]);

        // Neighbour cell N (stencil: alphaO, alphaN, alphaNp1)
        reconstruct1d(alphaO, alphaN, alphaNp1, beta_, eps_,
                      qimhM_N[facei], qiphM_N[facei],
                      qimhT_N[facei], qiphT_N[facei]);

        {
            scalar tmp0, tmp1, tmp2, tmp3;
            reconstruct1d(alphaOm2, alphaOm1, alphaO, beta_, eps_,
                          tmp0, tmp1, tmp2, tmp3);
            qimhM_L[facei] = tmp0;
            qiphM_L[facei] = tmp1;
            qimhT_L[facei] = tmp2;
            qiphT_L[facei] = tmp3;
        }
        {
            scalar tmp0, tmp1, tmp2, tmp3;
            reconstruct1d(alphaN, alphaNp1, alphaNp2, beta_, eps_,
                          tmp0, tmp1, tmp2, tmp3);
            qimhM_R[facei] = tmp0;
            qiphM_R[facei] = tmp1;
            qimhT_R[facei] = tmp2;
            qiphT_R[facei] = tmp3;
        }
    }

    // -----------------------------------------------------------------------
    // Phase 2 (AMReX bxg pass): BVD selection + flux assembly for each face.
    //
    // AMReX Phase 3 (upwind flux):
    //   face (i-1/2): left = sel(i-1).qiph,  right = sel(i).qimh
    //
    // Mapping to our face data:
    //   At face f between owner O and neighbour N:
    //     "left"  upwind value = O's selected qiph  (owner  right face value)
    //     "right" upwind value = N's selected qimh  (neighbour left face value)
    //
    //   BVD for O (decides O's qiph):
    //     bvdSelect1d(O's {qimhM,qiphM,qimhT,qiphT},
    //                 L's qiph_{M,T},          // left  neighbor of O
    //                 N's qimh_{M,T},           // right neighbor of O
    //                 alphaOm1, alphaO, alphaN)
    //
    //   BVD for N (decides N's qimh):
    //     bvdSelect1d(N's {qimhM,qiphM,qimhT,qiphT},
    //                 O's qiph_{M,T},           // left  neighbor of N
    //                 R's qimh_{M,T},            // right neighbor of N
    //                 alphaO, alphaN, alphaNp1)
    // -----------------------------------------------------------------------

    scalarField& alphaPhiIn = alphaPhi_.primitiveFieldRef();

    for (label facei = 0; facei < nInternalFaces; ++facei)
    {
        const label owni = own[facei];
        const label neii = nei[facei];

        const scalar alphaO = clampAlpha(alphaIn[owni]);
        const scalar alphaN = clampAlpha(alphaIn[neii]);

        const vector d  = cellCentres[neii] - cellCentres[owni];
        const scalar hON = Foam::mag(d);
        const vector dn = (hON > SMALL) ? d/hON : vector::zero;

        const scalar slopeO = (gradAlphaIn[owni] & dn) * hON;
        const scalar slopeN = (gradAlphaIn[neii] & dn) * hON;

        const scalar alphaOm1 = clampAlpha(alphaO - slopeO);
        const scalar alphaNp1 = clampAlpha(alphaN + slopeN);

        const bool useTHINC_O = bvdSelect1d
        (
            qimhM_O[facei], qiphM_O[facei],
            qimhT_O[facei], qiphT_O[facei],
            qiphM_L[facei], qiphT_L[facei],
            qimhM_N[facei], qimhT_N[facei],
            alphaOm1, alphaO, alphaN,
            eps_, delta_
        );

        const bool useTHINC_N = bvdSelect1d
        (
            qimhM_N[facei], qiphM_N[facei],
            qimhT_N[facei], qiphT_N[facei],
            qiphM_O[facei], qiphT_O[facei],
            qimhM_R[facei], qimhT_R[facei],
            alphaO, alphaN, alphaNp1,
            eps_, delta_
        );

        const scalar alphaFacePos = useTHINC_O ? qiphT_O[facei] : qiphM_O[facei];
        const scalar alphaFaceNeg = useTHINC_N ? qimhT_N[facei] : qimhM_N[facei];

        const scalar alphaFace = (phiIn[facei] >= 0)
            ? alphaFacePos
            : alphaFaceNeg;

        alphaPhiIn[facei] = phiIn[facei] * alphaFace;
    }

    // -----------------------------------------------------------------------
    // Boundary faces: upwind with neighbour reconstruction where possible.
    // -----------------------------------------------------------------------
    forAll(alphaPhi_.boundaryField(), patchi)
    {
        fvsPatchScalarField& alphaPhip = alphaPhi_.boundaryFieldRef()[patchi];
        const fvsPatchScalarField& phip = phi_.boundaryField()[patchi];
        const fvPatchScalarField& alphap = alpha1_.boundaryField()[patchi];
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
{}


void Foam::musclThincBvdAdvection::advect()
{
    computeAlphaPhi();

    alpha1_ = alpha1_.oldTime()
        - mesh_.time().deltaT()
        * (fvc::div(alphaPhi_) - alpha1_.oldTime()*fvc::div(phi_));

    // Clip to physical bounds — the explicit unsplit scheme is not TVD in 3D
    // and can produce small overshoots that amplify without this guard.
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