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
    TAHI4HD

Description
    Calculates and reports wall shear stress for all patches, for the
    specified times when using RAS turbulence models.

    Default behaviour assumes operating in incompressible mode.
    Use the -compressible option for compressible RAS cases.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "incompressible/RAS/RASModel/RASModel.H"

#include "basicPsiThermo.H"
#include "compressible/RAS/RASModel/RASModel.H"
//#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    timeSelector::addOptions();

    #include "addRegionOption.H"

    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"


    volScalarField TAHI4HD
    (
        IOobject
        (
            "TAHI4HD",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "TAHI4HD",
            dimLength/sqr(dimTime),
            0.0
        )
    );

    int nTime = 1;

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        mesh.readUpdate();

        volScalarField helicity4HD
        (
            IOobject
            (
                "helicity4HD",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar
            (
                "helicity4HD",
                dimLength/sqr(dimTime),
                0.0
            )
        );

        IOobject UHeader
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if (UHeader.headerOk())
        {
            Info<< "Reading field U\n" << endl;
            volVectorField U(UHeader, mesh);
            helicity4HD = U & (fvc::curl(U));
            TAHI4HD += mag(helicity4HD);
            TAHI4HD /= nTime;
            Info<< "Writing average helicity to field " << TAHI4HD.name() << nl << endl;
            TAHI4HD.write();
            TAHI4HD *= nTime;
            nTime++;

        }
        else
        {
            Info<< "    no U field" << endl;
        }

        /*Info<< "Writing wall shear stress to field " << WSS4HD.name()
            << nl << endl;
        WSS4HD.write();*/




    }

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
