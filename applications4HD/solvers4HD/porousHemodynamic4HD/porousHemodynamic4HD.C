/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    pimpleFoam

Description
    Large time-step transient solver for incompressible turbulent flow using
    the PIMPLE (merged PISO-SIMPLE) algorithm.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

    Consistent formulation without time-step and relaxation dependence by Jasak
    and Tukovic

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"
#include "pimpleControl.H"
#include "porousZones.H"
// Windkessel includes:
#include "fixedFluxPressureFvPatchScalarField.H"
#include "scalarIOList.H"
#include "WKFunctions.C"
#include "WKBCFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    pimpleControl pimple(mesh);

#   include "createWindkessel.H" //Windkessel header file (has to be placed in this order since p depends on scalar list storage)
#   include "createFields.H"
#   include "initContinuityErrs.H"
#   include "createTimeControls.H"

    // Read properties for haematocrit transp
#    include "hematocritProperties.H"

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
#       include "readTimeControls.H"
#       include "CourantNo.H"
#       include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- PIMPLE loop
        while (pimple.loop())
        {
#           include "UEqn.H"


            // --- PISO loop
            while (pimple.correct())
            {
#               include "pEqn.H"
            }

#           include "HEqn.H"

            if (haemoSwitch.value() == 0 || runTime < haemoSwitchTime)
            {
                Info<< "Not Solving for H, Migration Model is inactive 1" << nl << endl;
            }
            else
            {
                if (!pimple.finalIter())

                {
                    HEqn.relax();                   // use these two lines
                    HEqn.solve().initialResidual(); // for underrelaxed solver
                }
                else
                {
                    HEqn.solve(); // or this one for no underrelaxation (unstable for SIMPLE)
                    Info<< "Last PIMPLE iteration, no underrelaxation!" << nl << endl;
                }
            }



            if (TActive.value()==1)
            {
#                include "TEqn.H"
            }

            if (CActive.value()==1)
            {
#                include "CEqn.H"
            }


            turbulence->correct();
        }

        /* Updating the Windkessel struct data structure*/
        if (WK_FLAG >0)
        {
            Info << nl << "Updating Windkessel values ..." << endl;
            execute_at_end(mesh, phi, WKpressures, p, windkesselProperties, runTime.value());
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
