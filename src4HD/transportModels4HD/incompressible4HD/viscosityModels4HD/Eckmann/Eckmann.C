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

#include "Eckmann.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(Eckmann, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        Eckmann,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::Eckmann::calcNu() const
{
    const volScalarField& H = U_.mesh().lookupObject<volScalarField>("H");
    const volScalarField& T = U_.mesh().lookupObject<volScalarField>("T");
    dimensionedScalar dimMu("dimMu", dimMass/dimLength/dimTime, 1.0);
    dimensionedScalar Tc("Tc", dimTemperature, 273.15);

    return (1e-3*dimMu/rho_)*exp((lambda_/strainRate())+etha_*H) * ( alpha_ + phi_/(1+exp(beta_*(T-Tc-kesi_))) );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::Eckmann::Eckmann
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    EckmannCoeffs_(viscosityProperties.subDict(typeName + "Coeffs")),
    lambda_(EckmannCoeffs_.lookup("lambda")),
    etha_(EckmannCoeffs_.lookup("etha")),
    alpha_(EckmannCoeffs_.lookup("alpha")),
    phi_(EckmannCoeffs_.lookup("phi")),
    beta_(EckmannCoeffs_.lookup("beta")),
    kesi_(EckmannCoeffs_.lookup("kesi")),
    rho_(EckmannCoeffs_.lookup("rho")),
    nu_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcNu()
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::Eckmann::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    EckmannCoeffs_ = viscosityProperties.subDict(typeName + "Coeffs");

    EckmannCoeffs_.lookup("lambda") >> lambda_;
    EckmannCoeffs_.lookup("etha") >> etha_;
    EckmannCoeffs_.lookup("alpha") >> alpha_;
    EckmannCoeffs_.lookup("phi") >> phi_;
    EckmannCoeffs_.lookup("beta") >> beta_;
    EckmannCoeffs_.lookup("kesi") >> kesi_;
    EckmannCoeffs_.lookup("rho") >> rho_;

    return true;
}


// ************************************************************************* //
