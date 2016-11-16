/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2006-2010 OpenCFD Ltd.
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
/*---------------------------------------------------------------------------*\
   IH-Cantabria 2015 (http://www.ihcantabria.com/en/)
   IHFOAM 2015 (http://ihfoam.ihcantabria.com/) 

   Author(s):  Javier Lopez Lara (jav.lopez@unican.es)
               Gabriel Barajas   (barajasg@unican.es)
\*---------------------------------------------------------------------------*/

#include "IH_Waves_InletAlphaFvPatchScalarField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "Random.H"

#include "waveFun.H"

//C++
#include <stdio.h>
#include <unistd.h>
#include <errno.h>
#include "SquareMatrix.H"
#include "vector.H"
#include "Matrix.H"

#include "mpi.h"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
IH_Waves_InletAlphaFvPatchScalarField::
IH_Waves_InletAlphaFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    wavePeriod_(-1),
    waveHeight_(-1),
    waveLength_(-1),
    waterDepth_(-1),
    initialDepth_(-1),
    wavePhase_(3.0*PI()/2.0),
    lambdaStokesV_(-1),
    mCnoidal_(-1),
    nPaddles_(1),
    tSmooth_(-1),
    waveDictName_("IHWavesDict"),
    waveType_("aaa"),
    waveTheory_("aaa"),
    allCheck_(true),
    waveDir_(0),
    RealwaterDepth_(-1)
{}


Foam::
IH_Waves_InletAlphaFvPatchScalarField::
IH_Waves_InletAlphaFvPatchScalarField
(
    const IH_Waves_InletAlphaFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    wavePeriod_(ptf.wavePeriod_),
    waveHeight_(ptf.waveHeight_),
    waveLength_(ptf.waveLength_),
    waterDepth_(ptf.waterDepth_),
    initialDepth_(ptf.initialDepth_),
    wavePhase_(ptf.wavePhase_),
    lambdaStokesV_(ptf.lambdaStokesV_),
    mCnoidal_(ptf.mCnoidal_),
    nPaddles_(ptf.nPaddles_),
    tSmooth_(ptf.tSmooth_),
    waveDictName_(ptf.waveDictName_),
    waveType_(ptf.waveType_),
    waveTheory_(ptf.waveTheory_),
    allCheck_(ptf.allCheck_),
    waveDir_(ptf.waveDir_),
    RealwaterDepth_(ptf.RealwaterDepth_)
{}


Foam::
IH_Waves_InletAlphaFvPatchScalarField::
IH_Waves_InletAlphaFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p,iF,dict),
    wavePeriod_(dict.lookupOrDefault<scalar>("wavePeriod",-1)),
    waveHeight_(dict.lookupOrDefault<scalar>("waveHeight",-1)),
    waveLength_(dict.lookupOrDefault<scalar>("waveLength",-1)),
    waterDepth_(dict.lookupOrDefault<scalar>("waterDepth",-1)),
    initialDepth_(dict.lookupOrDefault<scalar>("initialDepth",-1)),
    wavePhase_(dict.lookupOrDefault<scalar>("wavePhase",3.0*PI()/2.0)),
    lambdaStokesV_(dict.lookupOrDefault<scalar>("lambdaStokesV",-1)),
    mCnoidal_(dict.lookupOrDefault<scalar>("mCnoidal",-1)),
    nPaddles_(dict.lookupOrDefault<label>("nPaddles",1)),
    tSmooth_(dict.lookupOrDefault<scalar>("tSmooth",-1)),
    waveDictName_(dict.lookupOrDefault<word>("waveDict","IHWavesDict")),
    waveType_(dict.lookupOrDefault<word>("waveType","aaa")),
    waveTheory_(dict.lookupOrDefault<word>("waveTheory","aaa")),
    allCheck_(dict.lookupOrDefault<bool>("allCheck",true)),
    waveDir_(dict.lookupOrDefault<scalar>("waveDir",0)),
    RealwaterDepth_(dict.lookupOrDefault<scalar>("RealwaterDepth",-1))
{}


Foam::
IH_Waves_InletAlphaFvPatchScalarField::
IH_Waves_InletAlphaFvPatchScalarField
(
    const IH_Waves_InletAlphaFvPatchScalarField& ptf
)
:
    fixedValueFvPatchField<scalar>(ptf),
    wavePeriod_(ptf.wavePeriod_),
    waveHeight_(ptf.waveHeight_),
    waveLength_(ptf.waveLength_),
    waterDepth_(ptf.waterDepth_),
    initialDepth_(ptf.initialDepth_),
    wavePhase_(ptf.wavePhase_),
    lambdaStokesV_(ptf.lambdaStokesV_),
    mCnoidal_(ptf.mCnoidal_),
    nPaddles_(ptf.nPaddles_),
    tSmooth_(ptf.tSmooth_),
    waveDictName_(ptf.waveDictName_),
    waveType_(ptf.waveType_),
    waveTheory_(ptf.waveTheory_),
    allCheck_(ptf.allCheck_),
    waveDir_(ptf.waveDir_),
    RealwaterDepth_(ptf.RealwaterDepth_)
{}


Foam::
IH_Waves_InletAlphaFvPatchScalarField::
IH_Waves_InletAlphaFvPatchScalarField
(
    const IH_Waves_InletAlphaFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ptf, iF),
    wavePeriod_(ptf.wavePeriod_),
    waveHeight_(ptf.waveHeight_),
    waveLength_(ptf.waveLength_),
    waterDepth_(ptf.waterDepth_),
    initialDepth_(ptf.initialDepth_),
    wavePhase_(ptf.wavePhase_),
    lambdaStokesV_(ptf.lambdaStokesV_),
    mCnoidal_(ptf.mCnoidal_),
    nPaddles_(ptf.nPaddles_),
    tSmooth_(ptf.tSmooth_),
    waveDictName_(ptf.waveDictName_),
    waveType_(ptf.waveType_),
    waveTheory_(ptf.waveTheory_),
    allCheck_(ptf.allCheck_),
    waveDir_(ptf.waveDir_),
    RealwaterDepth_(ptf.RealwaterDepth_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::IH_Waves_InletAlphaFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Aux. values
    scalar auxiliar = 0;
    scalar auxiliarTotal = 0;
    scalar X0 = 0;
    scalarField patchXsolit;

    const vector cMin = gMin(patch().patch().localPoints());
    const vector cMax = gMax(patch().patch().localPoints());
    const vector cSpan = cMax - cMin;
    const scalar zSpan = cSpan[2];

    scalar dMin = 0.0;
    scalar dSpan = 0.0;
    const scalarField patchD = patchDirection( cSpan, &dMin, &dSpan );

    // Variables & constants
    const volScalarField& alpha = 
        db().lookupObject<volScalarField>(alphaName());
    const volVectorField& U = db().lookupObject<volVectorField>("U");
    const fvMesh& mesh = alpha.mesh();
    const word& patchName = this->patch().name();
    const label patchID = mesh.boundaryMesh().findPatchID(patchName);
    const label nF = patch().faceCells().size();
    labelList cellGroup = Foam::labelList(nF, 1);
    const scalarField alphaCell = 
        alpha.boundaryField()[patchID].patchInternalField();
    const vectorField UCell = U.boundaryField()[patchID].patchInternalField();
    scalarField patchVOF = alphaCell; 
    const labelList celdas = patch().faceCells();
    const scalarField patchHeight = patch().Cf().component(2);

    // Calculate Z bounds of the faces
    scalarField zSup, zInf;
    faceBoundsZ( &zSup, &zInf );

    // Wave variables
    scalar waveOmega;
    scalarList waveOmegas;
    scalar waveK;
    scalarList waveKs;
    scalar waveAngle;
    scalar waveKx;
    scalarList waveKxs;
    scalar waveKy;

    // Define dictionary
    IOdictionary IHWavesDict
    (
        IOobject
        (
            waveDictName_,
            this->db().time().constant(),
            this->db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        )
    );

    scalar currTime = this->db().time().value();

    // Check for errors - Just the first time
    if (allCheck_)
    {
        waveType_ = (IHWavesDict.lookupOrDefault<word>("waveType","aaa")); 
        tSmooth_ = (IHWavesDict.lookupOrDefault<scalar>("tSmooth",-1.0));
	waterDepth_ = 
	    (IHWavesDict.lookupOrDefault<scalar>("waterDepth",-1.0));
	initialDepth_ = 
            (IHWavesDict.lookupOrDefault<scalar>("initialDepth",-1.0));
    	nPaddles_ = (IHWavesDict.lookupOrDefault<label>("nPaddles",1));
    } 

    // Grouping part
    scalarList dBreakPoints = Foam::scalarList(nPaddles_+1, dMin); 
    scalarList xGroup = Foam::scalarList(nPaddles_, 0.0);
    scalarList yGroup = Foam::scalarList(nPaddles_, 0.0);

    for (label i=0; i<nPaddles_; i++)
    {
        // Breakpoints, X & Y centre of the paddles
        dBreakPoints[i+1] = dMin + dSpan/(nPaddles_)*(i+1);
        xGroup[i] = cMin[0] + cSpan[0]/(2.0*nPaddles_) 
	    + cSpan[0]/(nPaddles_)*i;
        yGroup[i] = cMin[1] + cSpan[1]/(2.0*nPaddles_)
	    + cSpan[1]/(nPaddles_)*i;
    }

    forAll(patchD, patchCells) 
    {
        for (label i=0; i<nPaddles_; i++)
        {
            if ( (patchD[patchCells]>=dBreakPoints[i])
	        && (patchD[patchCells]<dBreakPoints[i+1]) )
            {
                cellGroup[patchCells] = i+1;
                continue;
            }
        }      
    }

    // Check for errors - Just the first time
    if (allCheck_)
    {
	if (RealwaterDepth_ == -1.0)
	{
	    if (waterDepth_ == -1.0)
	    {
	        RealwaterDepth_ = 
		    calcWL( alphaCell, cellGroup, zSpan )[0];
		if (initialDepth_ !=-1.0)
		{
			RealwaterDepth_ = RealwaterDepth_
			    + initialDepth_;
		}
	    }
	    else if ( waterDepth_ != -1.0 )
	    {
		RealwaterDepth_ = waterDepth_;
	    }
	}
    }

    if (allCheck_)
    {
	waveTheory_ = (IHWavesDict.lookupOrDefault<word>("waveTheory","aaa"));
    	waveHeight_ = (IHWavesDict.lookupOrDefault<scalar>("waveHeight",-1));
    	wavePeriod_ = (IHWavesDict.lookupOrDefault<scalar>("wavePeriod",-1));
    	waveDir_ = (IHWavesDict.lookupOrDefault<scalar>("waveDir",0));
    	wavePhase_ = 
	    (IHWavesDict.lookupOrDefault<scalar>("wavePhase",3.0*PI()/2.0));
	
	if ( waveType_ == "regular" )
	{
	   waveLength_ = StokesIFun::waveLength (RealwaterDepth_, wavePeriod_);
	}
    }

    if ( waveType_ == "regular" )
    {
        waveOmega = (2.0*PI())/wavePeriod_;
        waveK = 2.0*PI()/waveLength_;

        waveAngle = waveDir_*PI()/180.0;
        waveKx = waveK*cos(waveAngle);
        waveKy = waveK*sin(waveAngle);
    }
    else if ( waveType_ == "solitary" )
    {
        waveAngle = waveDir_*PI()/180.0;
        patchXsolit = patch().Cf().component(0)*cos(waveAngle) 
	    + patch().Cf().component(1)*sin(waveAngle);
        X0 = gMin(patchXsolit);
    }

    scalar timeMult = 1.0;

    if ( tSmooth_ > 0 && currTime < tSmooth_ )
    {
	timeMult = currTime/tSmooth_;
    }

    if (allCheck_)
    {
        if ( waveType_ == "regular" )
        {
            #include "checkInputErrorsRegular.H"
        }
        else if ( waveType_ == "solitary" )
        {
            #include "checkInputErrorsSolitary.H"
        }
        else
        {
            FatalError << "Wave type not supported, use:\n" 
	        << "regular, solitary" << exit(FatalError);
        }
	allCheck_ = false; 
    }

    scalarList calculatedLevel (nPaddles_,0.0);

    if ( waveType_ == "regular" )
    {
	if (waveDir_ == 0)
	{
           #include "calculatedLevelRegularNormal.H"
	}
	else
	{
	   #include "calculatedLevelRegular.H"
	}
    }
    else if ( waveType_ == "solitary" )
    {
        #include "calculatedLevelSolitary.H"
    }

    forAll(patchHeight, cellIndex)    
    {
	auxiliarTotal = 0;
	auxiliar = 0;
	
	if (zSup[cellIndex] < calculatedLevel[cellGroup[cellIndex]-1])
	{
	    patchVOF[cellIndex] = 1.0;
	}
	else if ( zInf[cellIndex] > calculatedLevel[cellGroup[cellIndex]-1])
	{
	    patchVOF[cellIndex] = 0.0;
	}
	else 
	{
	    auxiliar = calculatedLevel[cellGroup[cellIndex]-1]
	        - zInf[cellIndex];
	    auxiliarTotal = zSup[cellIndex]-zInf[cellIndex];
	    auxiliarTotal = auxiliar/auxiliarTotal;
	    patchVOF[cellIndex] = auxiliarTotal;
	}
    }

    operator == (patchVOF);

    fixedValueFvPatchField<scalar>::updateCoeffs();
}


void Foam::IH_Waves_InletAlphaFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);

    os.writeKeyword("waveDictName") << waveDictName_ << 
        token::END_STATEMENT << nl;

    os.writeKeyword("RealwaterDepth") << RealwaterDepth_ << 
        token::END_STATEMENT << nl;

    os.writeKeyword("initialDepth") << initialDepth_ << 
        token::END_STATEMENT << nl;

    writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchScalarField,
       IH_Waves_InletAlphaFvPatchScalarField
   );

}


// ************************************************************************* //
