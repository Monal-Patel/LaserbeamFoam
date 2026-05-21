if (pimple.nCorrPIMPLE() > 1)
{
    if (pimple.firstIter())
    {
        alpha1.storePrevIter();
        phi.storePrevIter();
        U.storePrevIter();
    }
    else
    {
        alpha1 = alpha1.prevIter();

        U = 0.5*U.prevIter() + 0.5*U;
        phi = 0.5*phi.prevIter() + 0.5*phi;
    }
}

if (nAlphaSubCycles > 1)
{
    dimensionedScalar totalDeltaT = runTime.deltaT();
    surfaceScalarField rhoPhiSum
    (
        mesh.newIOobject("rhoPhiSum"),
        mesh,
        dimensionedScalar(rhoPhi.dimensions(), Zero)
    );

    for
    (
        subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
        !(++alphaSubCycle).end();
    )
    {
        #include "alphaEqn.H"
        rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
    }

    rhoPhi = rhoPhiSum;
}
else
{
    #include "alphaEqn.H"
}

if (!pimple.firstIter())
{
    U = 2.0*U - U.prevIter();
    phi = 2.0*phi - phi.prevIter();
}

rho == alpha1*rho1 + alpha2*rho2;
