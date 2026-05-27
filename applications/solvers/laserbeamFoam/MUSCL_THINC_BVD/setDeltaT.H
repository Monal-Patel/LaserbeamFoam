/*---------------------------------------------------------------------------*\
Description
    Timestep control for MUSCL-THINC-BVD unsplit advection.
    Uses only maxCo (velocity CFL), ignores maxAlphaCo (no subcycling).

\*---------------------------------------------------------------------------*/

if (adjustTimeStep)
{
    scalar maxDeltaTFact = maxCo/(CoNum + SMALL);

    const scalar deltaTFact =
        Foam::min(Foam::min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

    runTime.setDeltaT
    (
        Foam::min
        (
            deltaTFact*runTime.deltaTValue(),
            maxDeltaT
        )
    );

    Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
}

// ************************************************************************* //
