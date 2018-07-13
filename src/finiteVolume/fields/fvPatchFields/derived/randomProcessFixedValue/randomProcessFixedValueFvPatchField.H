/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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
    Foam::randomProcessFixedValueFvPatchField

Description
    Foam::randomProcessFixedValueFvPatchField

SourceFiles
    randomProcessFixedValueFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef randomProcessFixedValueFvPatchField_H
#define randomProcessFixedValueFvPatchField_H

#include "Random.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class randomProcessFixedValueFvPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class randomProcessFixedValueFvPatchField
:
    public fixedValueFvPatchField<Type>
{
    // Private data

        //- Background flow
        Field<Type> U0_;

        //- Random gust direction
        Field<Type> gustDir_;

        //- Gust scale factor
        scalar scaleFac_;

        //- Power spectral density (w, S(w))
        List<scalar> circFreq_, psd_;

        //- Random number generator
        Random randGen_;

        //- Random phase container
        List<scalar> phi_;

        //- Current time index
        label curTimeIndex_;


    // Private member functions

        //- Return uniformly distributed random phase in [0, 2*pi]
        List<scalar> randomPhase();

        //- Return random process realisation
        scalar generate();


public:

    //- Runtime type information
    TypeName("randomProcessFixedValue");


    // Constructors

        //- Construct from patch and internal field
        randomProcessFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        randomProcessFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given randomProcessFixedValueFvPatchField
        //  onto a new patch
        randomProcessFixedValueFvPatchField
        (
            const randomProcessFixedValueFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        randomProcessFixedValueFvPatchField
        (
            const randomProcessFixedValueFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new randomProcessFixedValueFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        randomProcessFixedValueFvPatchField
        (
            const randomProcessFixedValueFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new randomProcessFixedValueFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the beackground flow
            const Field<Type>& U0() const
            {
                return U0_;
            }

            //- Return reference to the background flow to allow adjustment
            Field<Type>& U0()
            {
                return U0_;
            }

            //- Return the gust direction
            const Field<Type>& gustDir() const
            {
                return gustDir_;
            }

            Field<Type>& gustDir()
            {
                return gustDir_;
            }

            //- Return scale factor
            scalar scaleFac() const
            {
                return scaleFac_;
            }

            scalar& scaleFac()
            {
                return scaleFac_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<Type>&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "randomProcessFixedValueFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
