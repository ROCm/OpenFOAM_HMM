/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "waveMakerFlapPointPatchVectorField.H"
#include "mathematicalConstants.H"
#include "pointPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveMakerFlapPointPatchVectorField::
waveMakerFlapPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(p, iF),
    initialDepth_(0.0),
    wavePeriod_(0.0),
    waveHeigth_(0.0),
    waveLength_(0.0),
    wavePhase_(0.0),
    waveNumber_(0.0),
    rampTime_(0.0),
    g_(Zero),
    secondOrder_(false)
{}


Foam::waveMakerFlapPointPatchVectorField::
waveMakerFlapPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<vector>(p, iF, dict),
    initialDepth_(readScalar(dict.lookup("initialDepth"))),
    wavePeriod_(readScalar(dict.lookup("wavePeriod"))),
    waveHeigth_(readScalar(dict.lookup("waveHeigth"))),
    waveLength_(readScalar(dict.lookup("waveLength"))),
    wavePhase_(readScalar(dict.lookup("wavePhase"))),
    waveNumber_(readScalar(dict.lookup("waveNumber"))),
    rampTime_(readScalar(dict.lookup("rampTime"))),
    g_(dict.lookup("g")),
    secondOrder_(dict.lookupOrDefault<bool>("secondOrder",false))
{
    if (!dict.found("value"))
    {
        updateCoeffs();
    }
}


Foam::waveMakerFlapPointPatchVectorField::
waveMakerFlapPointPatchVectorField
(
    const waveMakerFlapPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<vector>(ptf, p, iF, mapper),
    initialDepth_(ptf.initialDepth_),
    wavePeriod_(ptf.wavePeriod_),
    waveHeigth_(ptf.waveHeigth_),
    waveLength_(ptf.waveLength_),
    wavePhase_(ptf.wavePhase_),
    waveNumber_(ptf.waveNumber_),
    rampTime_(ptf.rampTime_),
    g_(ptf.g_),
    secondOrder_(ptf.secondOrder_)
{}


Foam::waveMakerFlapPointPatchVectorField::
waveMakerFlapPointPatchVectorField
(
    const waveMakerFlapPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(ptf, iF),
    initialDepth_(ptf.initialDepth_),
    wavePeriod_(ptf.wavePeriod_),
    waveHeigth_(ptf.waveHeigth_),
    waveLength_(ptf.waveLength_),
    wavePhase_(ptf.wavePhase_),
    waveNumber_(ptf.waveNumber_),
    rampTime_(ptf.rampTime_),
    g_(ptf.g_),
    secondOrder_(ptf.secondOrder_)
{}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::scalar Foam::waveMakerFlapPointPatchVectorField::waveLength 
(
    const scalar h, 
    const scalar T
)
{
    const scalar L0 = mag(g_)*T*T/(constant::mathematical::twoPi);
    scalar L = L0;

    for(int i=1; i<=100; i++)
    {
            L = L0*tanh(constant::mathematical::twoPi*h/L);
    }
    return L;
}

Foam::scalar Foam::waveMakerFlapPointPatchVectorField::timeCoeff
(
    const scalar t
) const
{
    return max(0, min(t/rampTime_, 1));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::waveMakerFlapPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const polyMesh& mesh = this->internalField().mesh()();
    const Time& t = mesh.time();

    // Time ramp weight
    const scalar tCoeff = timeCoeff(t.value());

    vectorField localPoints_ = this->patch().localPoints();
    vectorField waveBoardMotion_ = 0.0*localPoints_;

    waveLength_ = waveLength (initialDepth_, wavePeriod_);    

    const scalar waveK = constant::mathematical::twoPi/waveLength_;
    const scalar sigma_ = (2.0*constant::mathematical::pi) / wavePeriod_;

    //First order
    if (secondOrder_== false)
    {
        forAll (localPoints_, labelP)
	{
	    scalar waveBoardStroke_ = (sinh(2.0*waveK*initialDepth_)
	        + 2.0*waveK*initialDepth_) 
		/ (4.0*sinh(waveK*initialDepth_)) 
		* (1.0 / ( sinh(waveK*initialDepth_) 
		+ (1.0-cosh(waveK*initialDepth_))
		/ (waveK*initialDepth_) ) ) * waveHeigth_;

	    waveBoardMotion_[labelP].component(0) = tCoeff
	        * (1.0+(localPoints_[labelP].component(2)-initialDepth_)
		/ (initialDepth_))*(waveBoardStroke_/2.0)*sin(sigma_*t.value()); 
	}

	Field<vector>::operator=
	(
	    waveBoardMotion_
	);
    }
    else if ( secondOrder_ == true)
    {
	scalar m1_;
	forAll (localPoints_, labelP)
	{
	    m1_ = ((4.0*sinh(waveK*initialDepth_))
	        / (sinh(2.0*waveK*initialDepth_)+2.0*waveK*initialDepth_))  
		* ( sinh(waveK*initialDepth_) + 1.0/(waveK*initialDepth_)
		* (1.0-cosh(waveK*initialDepth_)));

	    waveBoardMotion_[labelP].component(0) = waveHeigth_/(2.0*m1_)*sin(sigma_*t.value()) 
	        + pow(waveHeigth_,2)/(16.0*initialDepth_)*(3.0*cosh(waveK*initialDepth_)
		/ pow(sinh(waveK*initialDepth_),3) -2.0/m1_)*sin(2.0*sigma_*t.value());

	    waveBoardMotion_[labelP].component(0) = tCoeff
	        * (1.0+(localPoints_[labelP].component(2)-initialDepth_)
		/ (initialDepth_)) * waveBoardMotion_[labelP].component(0);
	}

 	Field<vector>::operator=
	(
	    waveBoardMotion_
	);
    }

    fixedValuePointPatchField<vector>::updateCoeffs();
}


void Foam::waveMakerFlapPointPatchVectorField::write(Ostream& os) const
{
    pointPatchField<vector>::write(os);
    os.writeEntry("initialDepth", initialDepth_);
    os.writeEntry("wavePeriod", wavePeriod_);
    os.writeEntry("waveHeigth", waveHeigth_);
    os.writeEntry("waveLength", waveLength_);
    os.writeEntry("wavePhase", wavePhase_);
    os.writeEntry("waveNumber", waveNumber_);
    os.writeEntry("rampTime", rampTime_);
    os.writeEntry("g", g_);
    os.writeEntry("secondOrder", secondOrder_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePointPatchTypeField
    (
        pointPatchVectorField,
        waveMakerFlapPointPatchVectorField
    );
}

// ************************************************************************* //
