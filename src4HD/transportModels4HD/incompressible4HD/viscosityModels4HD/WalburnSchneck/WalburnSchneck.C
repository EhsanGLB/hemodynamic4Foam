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

#include "WalburnSchneck.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(WalburnSchneck, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        WalburnSchneck,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::WalburnSchneck::calcNu() const
{
    const volScalarField& H = U_.mesh().lookupObject<volScalarField>("H");
    dimensionedScalar dimMu("dimMu", dimMass/dimLength/dimTime, 1.0);
    dimensionedScalar dimSR("dimMu", dimless/dimTime, 1.0);

    return (dimMu/rho_)*a1_*exp(a2_*H)*exp(a4_*(TPMA_/pow(H,2)))*pow(strainRate()/dimSR, -1*a3_*H);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::WalburnSchneck::WalburnSchneck
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    WalburnSchneckCoeffs_(viscosityProperties.subDict(typeName + "Coeffs")),
    a1_(WalburnSchneckCoeffs_.lookup("a1")),
    a2_(WalburnSchneckCoeffs_.lookup("a2")),
    a3_(WalburnSchneckCoeffs_.lookup("a3")),
    a4_(WalburnSchneckCoeffs_.lookup("a4")),
    TPMA_(WalburnSchneckCoeffs_.lookup("TPMA")),
    rho_(WalburnSchneckCoeffs_.lookup("rho")),
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

bool Foam::viscosityModels::WalburnSchneck::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    WalburnSchneckCoeffs_ = viscosityProperties.subDict(typeName + "Coeffs");

    WalburnSchneckCoeffs_.lookup("a1") >> a1_;
    WalburnSchneckCoeffs_.lookup("a2") >> a2_;
    WalburnSchneckCoeffs_.lookup("a3") >> a3_;
    WalburnSchneckCoeffs_.lookup("a4") >> a4_;
    WalburnSchneckCoeffs_.lookup("TPMA") >> TPMA_;
    WalburnSchneckCoeffs_.lookup("rho") >> rho_;

    return true;
}


// ************************************************************************* //
