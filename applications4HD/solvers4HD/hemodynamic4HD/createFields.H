Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field H\n" <<endl;
volScalarField H
(
    IOobject
    (
        "H",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field C\n" <<endl;
volScalarField C
(
    IOobject
    (
        "C",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);




#include "createPhi.H"



label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.schemesDict().setFluxRequired(p.name());

singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

Info<< "Reading field rAU if present\n" << endl;
volScalarField rAU
(
    IOobject
    (
        "rAU",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    runTime.deltaT(),
    zeroGradientFvPatchScalarField::typeName
);



dimensionedScalar rho(laminarTransport.lookup("rho"));
dimensionedScalar Prt(laminarTransport.lookup("Prt"));
dimensionedScalar k(laminarTransport.lookup("k"));
dimensionedScalar Cp(laminarTransport.lookup("Cp"));
dimensionedScalar alpha = k/(rho*Cp);
dimensionedScalar dp(laminarTransport.lookup("dp"));
dimensionedScalar kp(laminarTransport.lookup("kp"));
dimensionedScalar kB("kB", dimensionSet(1, 2, -2, -1, 0, 0, 0), 1.38064e-23);
dimensionedScalar TActive(laminarTransport.lookup("TActive"));
dimensionedScalar CActive(laminarTransport.lookup("CActive"));




