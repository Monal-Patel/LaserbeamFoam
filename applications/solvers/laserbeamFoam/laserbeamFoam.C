/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

Application
    laserbeamFoam

Group
    grpMultiphaseSolvers

Description
    Ray-Tracing heat source implementation with two phase incompressible VoF
    description of the metallic substrate and shielding gas phase,
    with optional mesh motion and mesh topology changes including adaptive
    re-meshing.
Authors

    Tom Flint, UoM.
    Philip Cardiff, UCD.
    Gowthaman Parivendhan, UCD.
    Joe Robson, UoM.
    Petar Cosic, UCD
    Simon Rodriguez, UCD

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "isoAdvection.H"
#include "thincAdvection.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "incompressibleInterPhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"
#include "dynamicRefineFvMesh.H"

#include "Polynomial.H"
#include "laserHeatSource.H"
#include "mthdModel.H"
#include "chrono.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for two incompressible, isothermal immiscible fluids"
        " using VOF phase-fraction based interface capturing.\n"
        "With optional mesh motion and mesh topology changes including"
        " adaptive re-meshing."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "MULES/createAlphaFluxes.H"
    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"

    if (interfaceTrackingScheme == "MULES")
    {
        if (!LTS)
        {
            #include "MULES/CourantNo.H"
            #include "setInitialDeltaT.H"
        }
    }
    else if (interfaceTrackingScheme == "isoAdvector")
    {
        #include "isoAdvector/porousCourantNo.H"
        #include "setInitialDeltaT.H"
    }
    else if (interfaceTrackingScheme == "THINC")
    {
        if (!LTS)
        {
            #include "MULES/CourantNo.H"
            #include "setInitialDeltaT.H"
        }
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

#ifdef LASERFOAM_PROFILING
    laserfoam::CumulativeTimer eqnTimer({"alpha", "laser", "U", "T", "p"});
#endif
    while (runTime.run())
    {
        #include "readControls.H"
        #include "readDyMControls.H"

        if (interfaceTrackingScheme == "MULES")
        {
            if (LTS)
            {
                #include "MULES/setRDeltaT.H"
            }
            else
            {
                #include "MULES/CourantNo.H"
                #include "MULES/alphaCourantNo.H"
                #include "MULES/setDeltaT.H"
            }
        }
        else if (interfaceTrackingScheme == "isoAdvector")
        {
            #include "isoAdvector/porousCourantNo.H"
            #include "isoAdvector/porousAlphaCourantNo.H"
            #include "isoAdvector/setDeltaT.H"
        }
        else if (interfaceTrackingScheme == "THINC")
        {
            if (LTS)
            {
                #include "MULES/setRDeltaT.H"
            }
            else
            {
                #include "MULES/CourantNo.H"
                #include "MULES/alphaCourantNo.H"
                #include "MULES/setDeltaT.H"
            }
        }

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {

#ifdef LASERFOAM_PROFILING
            auto alphaT0 = eqnTimer.tick();
#endif
            if (interfaceTrackingScheme == "MULES")
            {
                #include "MULES/firstIter.H"
                #include "MULES/alphaControls.H"
                #include "MULES/alphaEqnSubCycle.H"
            }
            else if (interfaceTrackingScheme == "isoAdvector")
            {
                #include "isoAdvector/firstIter.H"
                #include "isoAdvector/alphaControls.H"
                #include "isoAdvector/alphaEqnSubCycle.H"
            }
            else if (interfaceTrackingScheme == "THINC")
            {
                #include "THINC/firstIter.H"
                #include "THINC/alphaControls.H"
                #include "THINC/alphaEqnSubCycle.H"
            }
#ifdef LASERFOAM_PROFILING
            eqnTimer.accumulate("alpha", alphaT0);
#endif

            #include "updateProps.H"

#ifdef LASERFOAM_PROFILING
            auto laserT0 = eqnTimer.tick();
#endif
            // Update the laser deposition field
            laser.updateDeposition
            (
                alpha_filtered, n_filtered, electrical_resistivity
            );
#ifdef LASERFOAM_PROFILING
            eqnTimer.accumulate("laser", laserT0);
#endif

            mixture.correct();

            if (pimple.frozenFlow())
            {
                continue;
            }

#ifdef LASERFOAM_PROFILING
            auto uT0 = eqnTimer.tick();
#endif
            #include "UEqn.H"
            if (mthd.valid())
            {
                mthd->solve(phi, U);
            }
#ifdef LASERFOAM_PROFILING
            eqnTimer.accumulate("U", uT0);
#endif
#ifdef LASERFOAM_PROFILING
            auto tT0 = eqnTimer.tick();
#endif
            #include "TEqn.H"
#ifdef LASERFOAM_PROFILING
            eqnTimer.accumulate("T", tT0);
#endif

            // --- Pressure corrector loop
#ifdef LASERFOAM_PROFILING
            auto pT0 = eqnTimer.tick();
#endif
            while (pimple.correct())
            {
                #include "pEqn.H"
            }
#ifdef LASERFOAM_PROFILING
            eqnTimer.accumulate("p", pT0);
#endif

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        // Update the melt history
        const volScalarField& alphaMetal =
            mesh.lookupObject<volScalarField>("alpha.metal");
        condition = pos(alphaMetal - 0.5) * pos(epsilon1 - 0.5);
        meltHistory += condition;

        runTime.write();

        // Write ray paths to VTK files
        if (runTime.outputTime())
        {
            laser.writeRayPathsToVTK();
            laser.writeRayPathVTKSeriesFile();
        }

        runTime.printExecutionTime(Info);
    }

#ifdef LASERFOAM_PROFILING
    eqnTimer.report();
#endif
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
