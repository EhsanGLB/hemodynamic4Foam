/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v1812                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*\
Description
    Windkessel pressure outflow boundary condition, 3 element.

SourceFiles
    WKBC4HDFvPatchVectorField.C

Authors
    Andris Piebalgs, Imperial College London, 2017  All rights reserved
    Boram Gu, Emily Manchester, Imperial College London, 2020 All rights reserved
\*---------------------------------------------------------------------------*/

#include "WKBC4HD.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "mathematicalConstants.H"
#include "scalarIOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

WKBC4HD::WKBC4HD
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    index_(0)
{}


WKBC4HD::WKBC4HD
(
    const WKBC4HD& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    index_(ptf.index_)
{}


WKBC4HD::WKBC4HD
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    index_(readScalar(dict.lookup("index")))
{
    evaluate();
}


WKBC4HD::WKBC4HD
(
    const WKBC4HD& fcvpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(fcvpvf, iF),
    index_(fcvpvf.index_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void WKBC4HD::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    
    /*Creating a patch field of same size as the boundary field*/
    const fvPatch& p = this->patch();
    scalarField report(p.size());


    /* Accessing the variables stored in mesh */
    const fvMesh& mesh = patch().boundaryMesh().mesh();
    const scalarIOList& WKpressures = mesh.lookupObject<scalarIOList>("WKpressures");
    const scalar current_pressure = WKpressures[index_];

    //Info << "store in WKBC4HD " << store << nl << endl;
    //Info << "current_pressure" << current_pressure << nl << endl;

    /*Applying the pressure to each face on the outlet*/
    forAll(report, it)
    {
        report[it] = current_pressure/1060;
    }

    //Info << "BC P " << average(report)*1060 << nl << endl;

    /*Assigning the operator to the new patch field*/
    scalarField::operator=(report);

}


void WKBC4HD::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("index") << index_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, WKBC4HD);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
