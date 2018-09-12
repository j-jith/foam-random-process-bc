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

\*---------------------------------------------------------------------------*/

#include "randomProcessPerfCorrFvPatchField.H"
#include "mathematicalConstants.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
List<scalar> randomProcessPerfCorrFvPatchField<Type>::randomPhase()
{
    List<scalar> phi(circFreq_.size(), 0.);

    forAll(phi, i)
    {
        phi[i] = randGen_.scalar01() * 2*mathematicalConstant::pi;
    }

    return phi;
}

template<class Type>
scalar randomProcessPerfCorrFvPatchField<Type>::generate()
{
    scalar x = 0.0;
    scalar dw = circFreq_[1]-circFreq_[0];

    forAll(circFreq_, i)
    {
        //x += sqrt(2*psd_[i]) * cos(circFreq_[i]*this->db().time().value() + randomPhase());
        x += sqrt(2*dw*psd_[i]) * cos(circFreq_[i]*this->db().time().value() + phi_[i]);
    }

    return x;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
randomProcessPerfCorrFvPatchField<Type>::randomProcessPerfCorrFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    U0_(p.size()),
    gustDir_(p.size()),
    scaleFac_(1.0),
    circFreq_(scalarList(2, 0.0)),
    psd_(scalarList(2, 0.0)),
    randGen_(label(0)),
    curTimeIndex_(-1)
{
    phi_ = randomPhase();
}


template<class Type>
randomProcessPerfCorrFvPatchField<Type>::randomProcessPerfCorrFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF),
    U0_("U0", dict, p.size()),
    gustDir_("gustDir", dict, p.size()),
    scaleFac_(readScalar(dict.lookup("scaleFac"))),
    circFreq_(dict.lookup("circFreq")),
    psd_(dict.lookup("psd")),
    randGen_(label(0)),
    curTimeIndex_(-1)
{
    phi_ = randomPhase();
    if (dict.found("value"))
    {
        fixedValueFvPatchField<Type>::operator==
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        fixedValueFvPatchField<Type>::operator==(U0_ + scaleFac_*gustDir_*generate());
    }
}


template<class Type>
randomProcessPerfCorrFvPatchField<Type>::randomProcessPerfCorrFvPatchField
(
    const randomProcessPerfCorrFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    U0_(ptf.U0_, mapper),
    gustDir_(ptf.gustDir_, mapper),
    scaleFac_(ptf.scaleFac_),
    circFreq_(ptf.circFreq_),
    psd_(ptf.psd_),
    randGen_(ptf.randGen_),
    curTimeIndex_(-1)
{
    phi_ = randomPhase();
}


template<class Type>
randomProcessPerfCorrFvPatchField<Type>::randomProcessPerfCorrFvPatchField
(
    const randomProcessPerfCorrFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    U0_(ptf.U0_),
    gustDir_(ptf.gustDir_),
    scaleFac_(ptf.scaleFac_),
    circFreq_(ptf.circFreq_),
    psd_(ptf.psd_),
    randGen_(ptf.randGen_),
    curTimeIndex_(-1)
{
    phi_ = randomPhase();
}


template<class Type>
randomProcessPerfCorrFvPatchField<Type>::randomProcessPerfCorrFvPatchField
(
    const randomProcessPerfCorrFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    U0_(ptf.U0_),
    gustDir_(ptf.gustDir_),
    scaleFac_(ptf.scaleFac_),
    circFreq_(ptf.circFreq_),
    psd_(ptf.psd_),
    randGen_(ptf.randGen_),
    curTimeIndex_(-1)
{
    phi_ = randomPhase();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void randomProcessPerfCorrFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<Type>::autoMap(m);
    U0_.autoMap(m);
    gustDir_.autoMap(m);
}


template<class Type>
void randomProcessPerfCorrFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<Type>::rmap(ptf, addr);

    const randomProcessPerfCorrFvPatchField<Type>& tiptf =
        refCast<const randomProcessPerfCorrFvPatchField<Type> >(ptf);

    U0_.rmap(tiptf.U0_, addr);
    gustDir_.rmap(tiptf.gustDir_, addr);
}


template<class Type>
void randomProcessPerfCorrFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        Field<Type>& patchField = *this;
        patchField = U0_ + scaleFac_*gustDir_*generate();

        //const vectorField& Cf = patch().Cf();

        //forAll(patchField, facei)
        //{
        //    patchField[facei] = generate() * gustDir_;
        //}

        //const vectorField Ui = generate() * ((Cf-Cf) + gustDir_);
        //vectorField::operator=(Ui);
        //patchField = Ui;

        curTimeIndex_ = this->db().time().timeIndex();
    }

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void randomProcessPerfCorrFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    U0_.writeEntry("U0", os);
    gustDir_.writeEntry("gustDir", os);
    os.writeKeyword("scaleFac")
        << scaleFac_ << token::END_STATEMENT << nl;
    os.writeKeyword("circFreq")
        << circFreq_ << token::END_STATEMENT << nl;
    os.writeKeyword("psd")
        << psd_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
