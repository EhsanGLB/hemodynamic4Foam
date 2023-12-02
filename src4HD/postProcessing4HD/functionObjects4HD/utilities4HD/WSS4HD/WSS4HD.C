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

\*---------------------------------------------------------------------------*/

#include "WSS4HD.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "compressible/turbulenceModel/turbulenceModel.H"
#include "incompressible/turbulenceModel/turbulenceModel.H"
#include "wallPolyPatch.H"
#include "singlePhaseTransportModel.H"
#include "RASModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(WSS4HD, 0);
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::WSS4HD::writeFileHeader(const label i)
{
    // Add headers to output data
    writeHeader(file(), "Wall shear stress");
    writeCommented(file(), "Time");
    writeTabbed(file(), "patch");
    writeTabbed(file(), "min");
    writeTabbed(file(), "max");
    file() << endl;
}


void Foam::WSS4HD::calcShearStress
(
    const fvMesh& mesh,
    const volSymmTensorField& Reff,
    volVectorField& shearStress
)
{
    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchI = iter.key();
        const polyPatch& pp = mesh.boundaryMesh()[patchI];

        vectorField& ssp = shearStress.boundaryField()[patchI];
        const vectorField& Sfp = mesh.Sf().boundaryField()[patchI];
        const scalarField& magSfp = mesh.magSf().boundaryField()[patchI];
        const symmTensorField& Reffp = Reff.boundaryField()[patchI];

        ssp = (-Sfp/magSfp) & Reffp;

        vector minSsp = gMin(ssp);
        vector maxSsp = gMax(ssp);

        if (Pstream::master())
        {
            file() << mesh.time().value()
                << token::TAB << pp.name()
                << token::TAB << minSsp
                << token::TAB << maxSsp
                << endl;
        }

        if (log_) Info<< "    min/max(" << pp.name() << ") = "
            << minSsp << ", " << maxSsp << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::WSS4HD::WSS4HD
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    functionObjectFile(obr, name, typeName),
    name_(name),
    obr_(obr),
    active_(true),
    log_(true),
    patchSet_(),
    UName_("U"),
    phiName_("phi")
{
    // Check if the available mesh is an fvMesh, otherwise deactivate
    if (!isA<fvMesh>(obr_))
    {
        active_ = false;
        WarningIn
        (
            "WSS4HD::WSS4HD"
            "("
                "const word&, "
                "const objectRegistry&, "
                "const dictionary&, "
                "const bool"
            ")"
        )   << "No fvMesh available, deactivating " << name_ << nl
            << endl;
    }

    if (active_)
    {
        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        volVectorField* WSS4HDPtr
        (
            new volVectorField
            (
                IOobject
                (
                    type(),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedVector
                (
                    "0",
                    sqr(dimLength)/sqr(dimTime),
                    vector::zero
                )
            )
        );

        mesh.objectRegistry::store(WSS4HDPtr);
    }

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::WSS4HD::~WSS4HD()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::WSS4HD::read(const dictionary& dict)
{
    if (active_)
    {
        UName_ = dict.lookupOrDefault<word>("UName", "U");
        phiName_ = dict.lookupOrDefault<word>("phiName", "phi");

        log_ = dict.lookupOrDefault<Switch>("log", true);

        const fvMesh& mesh = refCast<const fvMesh>(obr_);
        const polyBoundaryMesh& pbm = mesh.boundaryMesh();

        patchSet_ =
            mesh.boundaryMesh().patchSet
            (
                wordReList(dict.lookupOrDefault("patches", wordReList()))
            );

        Info<< type() << " " << name_ << ":" << nl;

        if (patchSet_.empty())
        {
            forAll(pbm, patchI)
            {
                if (isA<wallPolyPatch>(pbm[patchI]))
                {
                    patchSet_.insert(patchI);
                }
            }

            Info<< "    processing all wall patches" << nl << endl;
        }
        else
        {
            Info<< "    processing wall patches: " << nl;
            labelHashSet filteredPatchSet;
            forAllConstIter(labelHashSet, patchSet_, iter)
            {
                label patchI = iter.key();
                if (isA<wallPolyPatch>(pbm[patchI]))
                {
                    filteredPatchSet.insert(patchI);
                    Info<< "        " << pbm[patchI].name() << endl;
                }
                else
                {
                    WarningIn("void WSS4HD::read(const dictionary&)")
                        << "Requested wall shear stress on non-wall boundary "
                        << "type patch: " << pbm[patchI].name() << endl;
                }
            }

            Info<< endl;

            patchSet_ = filteredPatchSet;
        }
    }
}


void Foam::WSS4HD::execute()
{
    typedef compressible::turbulenceModel cmpModel;
    typedef incompressible::turbulenceModel icoModel;

    if (active_)
    {
        const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);
        const surfaceScalarField& phi = obr_.lookupObject<surfaceScalarField>(phiName_);

        functionObjectFile::write();

        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        volVectorField& WSS4HD =
            const_cast<volVectorField&>
            (
                mesh.lookupObject<volVectorField>(type())
            );

        if (log_) Info<< type() << " " << name_ << " output:" << nl;


        tmp<volSymmTensorField> Reff;
        /*if (mesh.foundObject<cmpModel>("turbulenceProperties"))
        {
            const cmpModel& model =
                mesh.lookupObject<cmpModel>("turbulenceProperties");

            Reff = model.devRhoReff();
        }
        else if (mesh.foundObject<icoModel>("turbulenceProperties"))
        {
            const icoModel& model =
                mesh.lookupObject<icoModel>("turbulenceProperties");

            Reff = model.devReff();
        }
        else
        {
            FatalErrorIn("void Foam::WSS4HD::execute()")
                << "Unable to find turbulence model in the "
                << "database" << exit(FatalError);
        }*/
        singlePhaseTransportModel laminarTransport(U, phi);
        dimensionedScalar rho(laminarTransport.lookup("rho"));
        autoPtr<Foam::incompressible::RASModel> model(Foam::incompressible::RASModel::New(U, phi, laminarTransport));
        Reff = rho*model->devReff();

        calcShearStress(mesh, Reff(), WSS4HD);
    }
}


void Foam::WSS4HD::end()
{
    if (active_)
    {
        execute();
    }
}


void Foam::WSS4HD::timeSet()
{
    // Do nothing
}


void Foam::WSS4HD::write()
{
    if (active_)
    {
        functionObjectFile::write();

        const volVectorField& WSS4HD =
            obr_.lookupObject<volVectorField>(type());

        if (log_) Info<< type() << " " << name_ << " output:" << nl
            << "    writing field " << WSS4HD.name() << nl
            << endl;

        WSS4HD.write();
    }
}


// ************************************************************************* //
