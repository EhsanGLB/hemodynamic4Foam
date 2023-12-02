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
    magWSS4HD

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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void calcIncompressible
(
    const fvMesh& mesh,
    const Time& runTime,
    const volVectorField& U,
    volScalarField& magWSS4HD
)
{
    IOobject muHeader
    (
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (!muHeader.headerOk())
    {
        Info<< "    no mu field" << endl;
        return;
    }

    Info<< "Reading field mu\n" << endl;
    volScalarField mu(muHeader, mesh);

    #include "createPhi.H"

    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::RASModel> model
    (
        incompressible::RASModel::New(U, phi, laminarTransport)
    );

    const volSymmTensorField Reff(model->devReff());

    dimensionedScalar rho(laminarTransport.lookup("rho"));

    forAll(magWSS4HD.boundaryField(), patchI)
    {
        magWSS4HD.boundaryField()[patchI] = rho.value() * mag(
        (
           -mesh.Sf().boundaryField()[patchI]
           /mesh.magSf().boundaryField()[patchI]
        ) & Reff.boundaryField()[patchI] );
    }
}


void calcCompressible
(
    const fvMesh& mesh,
    const Time& runTime,
    const volVectorField& U,
    volScalarField& magWSS4HD
)
{
    IOobject rhoHeader
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (!rhoHeader.headerOk())
    {
        Info<< "    no rho field" << endl;
        return;
    }

    Info<< "Reading field rho\n" << endl;
    volScalarField rho(rhoHeader, mesh);

    #include "compressibleCreatePhi.H"

    autoPtr<basicPsiThermo> pThermo
    (
        basicPsiThermo::New(mesh)
    );
    basicPsiThermo& thermo = pThermo();

    autoPtr<compressible::RASModel> model
    (
        compressible::RASModel::New
        (
            rho,
            U,
            phi,
            thermo
        )
    );

    const volSymmTensorField Reff(model->devRhoReff());

    forAll(magWSS4HD.boundaryField(), patchI)
    {
        magWSS4HD.boundaryField()[patchI] = mag(
        (
           -mesh.Sf().boundaryField()[patchI]
           /mesh.magSf().boundaryField()[patchI]
        ) & Reff.boundaryField()[patchI] );
    }
}


int main(int argc, char *argv[])
{
    timeSelector::addOptions();

    #include "addRegionOption.H"

    argList::validOptions.insert("compressible","");

    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"

    bool compressible = args.optionFound("compressible");


    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        mesh.readUpdate();

        volScalarField magWSS4HD
        (
            IOobject
            (
                "magWSS4HD",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar
            (
                "magWSS4HD",
                dimMass/dimLength/sqr(dimTime),
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

            if (compressible)
            {
                calcCompressible(mesh, runTime, U, magWSS4HD);
            }
            else
            {
                calcIncompressible(mesh, runTime, U, magWSS4HD);
            }
        }
        else
        {
            Info<< "    no U field" << endl;
        }

        Info<< "Writing magnitude of wall shear stress to field " << magWSS4HD.name() << nl << endl;

        magWSS4HD.write();
    }

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
