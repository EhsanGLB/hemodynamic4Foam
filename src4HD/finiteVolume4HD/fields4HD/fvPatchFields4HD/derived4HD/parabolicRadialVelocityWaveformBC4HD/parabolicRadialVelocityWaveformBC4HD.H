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

Class
    parabolicRadialVelocityWaveformBC4HD

Description
    Boundary condition specifies a power law velocity inlet profile
    (fixed value), given mean velocity value,
    flow direction n and direction of the power law coordinate y.

SourceFiles
    parabolicRadialVelocityWaveformBC4HD.C

Author
    Ehsan Golab, SUT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef parabolicRadialVelocityWaveformBC4HD_H
#define parabolicRadialVelocityWaveformBC4HD_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "complex.H"
#include "complex4HD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class parabolicRadialVelocityFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class parabolicRadialVelocityWaveformBC4HD
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Amplitude
        scalar amplitude_;

        //- Frequency
        scalar frequency_;

        //- Flow direction
        vector lenDir_;

        //- Radius
        scalar radius_;

        //- fourierSeriesType
        word fourierSeriesType_;

        //- cosCoeffs component
        scalarList cosCoeffs_;

        //- sinCoeffs component
        scalarList sinCoeffs_;

        //- Real component
        scalarList real_;

        //- Imaginary component
        scalarList imaginary_;



public:

    //- Runtime type information
    TypeName("parabolicRadialVelocityWaveformBC4HD");


    // Constructors

        //- Construct from patch and internal field
        parabolicRadialVelocityWaveformBC4HD
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        parabolicRadialVelocityWaveformBC4HD
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given parabolicRadialVelocityWaveformBC4HD
        //  onto a new patch
        parabolicRadialVelocityWaveformBC4HD
        (
            const parabolicRadialVelocityWaveformBC4HD&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new parabolicRadialVelocityWaveformBC4HD(*this)
            );
        }

        //- Construct as copy setting internal field reference
        parabolicRadialVelocityWaveformBC4HD
        (
            const parabolicRadialVelocityWaveformBC4HD&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new parabolicRadialVelocityWaveformBC4HD(*this, iF)
            );
        }


    // Member functions

        //- Calculate the velocityReal(time)
        virtual scalar velocityRealFunction(scalar t);

        //- Calculate the velocityComplex(time)
        virtual complex velocityComplexFunction(scalar t);

        //- Calculate the velocity(time)
        virtual scalar velocityFunction(scalar t);

        //- Update coefficients
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
