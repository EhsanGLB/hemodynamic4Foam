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

#include "Ameenuddin.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(Ameenuddin, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        Ameenuddin,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::Ameenuddin::calcNu() const
{
    const volScalarField& H = U_.mesh().lookupObject<volScalarField>("H");

    volScalarField nu0_ = 0.0736*(a1_ + a2_*H + a3_*pow(H, 2))/rho_;
    volScalarField nuInf_ = 0.005*(b1_ + b2_*H + b3_*pow(H, 2))/rho_;
    volScalarField lambda_ = 14.81*(c1_*H + c2_*pow(H, 2));

    return nuInf_ + (nu0_ - nuInf_) * ((1+log(1+lambda_*strainRate()))/(1+lambda_*strainRate()));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::Ameenuddin::Ameenuddin
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    AmeenuddinCoeffs_(viscosityProperties.subDict(typeName + "Coeffs")),
    a1_(AmeenuddinCoeffs_.lookup("a1")),
    a2_(AmeenuddinCoeffs_.lookup("a2")),
    a3_(AmeenuddinCoeffs_.lookup("a3")),
    b1_(AmeenuddinCoeffs_.lookup("b1")),
    b2_(AmeenuddinCoeffs_.lookup("b2")),
    b3_(AmeenuddinCoeffs_.lookup("b3")),
    c1_(AmeenuddinCoeffs_.lookup("c1")),
    c2_(AmeenuddinCoeffs_.lookup("c2")),
    rho_(AmeenuddinCoeffs_.lookup("rho")),
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

bool Foam::viscosityModels::Ameenuddin::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    AmeenuddinCoeffs_ = viscosityProperties.subDict(typeName + "Coeffs");

    AmeenuddinCoeffs_.lookup("a1") >> a1_;
    AmeenuddinCoeffs_.lookup("a2") >> a2_;
    AmeenuddinCoeffs_.lookup("a3") >> a3_;
    AmeenuddinCoeffs_.lookup("b1") >> b1_;
    AmeenuddinCoeffs_.lookup("b2") >> b2_;
    AmeenuddinCoeffs_.lookup("b3") >> b3_;
    AmeenuddinCoeffs_.lookup("c1") >> c1_;
    AmeenuddinCoeffs_.lookup("c2") >> c2_;
    AmeenuddinCoeffs_.lookup("rho") >> rho_;

    return true;
}


// ************************************************************************* //
