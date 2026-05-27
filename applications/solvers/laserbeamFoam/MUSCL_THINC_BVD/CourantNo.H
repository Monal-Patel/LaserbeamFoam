/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description
    Unsplit Courant number for MUSCL-THINC-BVD explicit advection.

    Computes max per-direction CFL: Co = max_i(|u_i|*dt/dx_i)
    Differs from MULES split CFL which sums directions.

\*---------------------------------------------------------------------------*/

scalar CoNum = 0.0;
scalar meanCoNum = 0.0;

{
    const surfaceScalarField magPhiByS(mag(phi)/mesh.magSf());

    scalarField maxPhiByDx(mesh.nCells(), 0.0);

    const labelUList& own = mesh.faceOwner();
    const labelUList& nei = mesh.faceNeighbour();
    const label nInternalFaces = mesh.nInternalFaces();

    // For each cell, track max(|phi|/area) across all faces
    // Equivalent to max(|u|/dx) since area ~ dx*dy, |phi| ~ |u|*dy
    for (label facei = 0; facei < nInternalFaces; ++facei)
    {
        const scalar val = magPhiByS[facei];
        maxPhiByDx[own[facei]] = max(maxPhiByDx[own[facei]], val);
        maxPhiByDx[nei[facei]] = max(maxPhiByDx[nei[facei]], val);
    }

    // Boundary faces
    forAll(mesh.boundary(), patchi)
    {
        const labelUList& fc = mesh.boundary()[patchi].faceCells();
        const scalarField& magPhiByS_p = magPhiByS.boundaryField()[patchi];

        forAll(fc, fi)
        {
            maxPhiByDx[fc[fi]] = max(maxPhiByDx[fc[fi]], magPhiByS_p[fi]);
        }
    }

    CoNum = gMax(maxPhiByDx) * runTime.deltaTValue();
    meanCoNum = gSum(maxPhiByDx) / scalar(mesh.nCells()) * runTime.deltaTValue();
}

Info<< "Courant Number mean: " << meanCoNum
    << " max: " << CoNum << endl;

// ************************************************************************* //
