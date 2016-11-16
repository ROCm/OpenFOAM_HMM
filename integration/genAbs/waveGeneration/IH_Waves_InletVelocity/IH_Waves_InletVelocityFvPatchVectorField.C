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

#include "IH_Waves_InletVelocityFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"

#include "waveFun.H"

#include <stdio.h>
#include <unistd.h>
#include <errno.h>
#include "SquareMatrix.H"
#include "vector.H"
#include "Matrix.H"

#include "mpi.h"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
IH_Waves_InletVelocityFvPatchVectorField::
IH_Waves_InletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    wavePeriod_(-1),
    waveHeight_(-1),
    waveLength_(-1),
    waterDepth_(-1),
    initialDepth_(-1),
    RealwaterDepth_(-1),
    wavePhase_(3.0*PI()/2.0),
    lambdaStokesV_(-1),
    mCnoidal_(-1),
    genAbs_(false),
    nPaddles_(1),
    tSmooth_(-1),
    leftORright_(-1),
    waveDictName_("IHWavesDict"),
    waveType_("aaa"),
    waveTheory_("aaa"),
    allCheck_(true),
    waveDir_(0)
{}


Foam::
IH_Waves_InletVelocityFvPatchVectorField::
IH_Waves_InletVelocityFvPatchVectorField
(
    const IH_Waves_InletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    wavePeriod_(ptf.wavePeriod_),
    waveHeight_(ptf.waveHeight_),
    waveLength_(ptf.waveLength_),
    waterDepth_(ptf.waterDepth_),
    initialDepth_(ptf.initialDepth_),
    RealwaterDepth_(ptf.RealwaterDepth_),
    wavePhase_(ptf.wavePhase_),
    lambdaStokesV_(ptf.lambdaStokesV_),
    mCnoidal_(ptf.mCnoidal_),
    genAbs_(ptf.genAbs_),
    nPaddles_(ptf.nPaddles_),
    tSmooth_(ptf.tSmooth_),
    leftORright_(ptf.leftORright_),
    waveDictName_(ptf.waveDictName_),
    waveType_(ptf.waveType_),
    waveTheory_(ptf.waveTheory_),
    allCheck_(ptf.allCheck_),
    waveDir_(ptf.waveDir_)
{}


Foam::
IH_Waves_InletVelocityFvPatchVectorField::
IH_Waves_InletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    wavePeriod_(dict.lookupOrDefault<scalar>("wavePeriod",-1)),
    waveHeight_(dict.lookupOrDefault<scalar>("waveHeight",-1)),
    waveLength_(dict.lookupOrDefault<scalar>("waveLength",-1)),
    waterDepth_(dict.lookupOrDefault<scalar>("waterDepth",-1)),
    initialDepth_(dict.lookupOrDefault<scalar>("initialDepth",-1)),
    RealwaterDepth_(dict.lookupOrDefault<scalar>("RealwaterDepth",-1)),
    wavePhase_(dict.lookupOrDefault<scalar>("wavePhase",3.0*PI()/2.0)),
    lambdaStokesV_(dict.lookupOrDefault<scalar>("lambdaStokesV",-1)),
    mCnoidal_(dict.lookupOrDefault<scalar>("mCnoidal",-1)),
    genAbs_(dict.lookupOrDefault<bool>("genAbs",false)),
    nPaddles_(dict.lookupOrDefault<label>("nPaddles",1)),
    tSmooth_(dict.lookupOrDefault<scalar>("tSmooth",-1)),
    leftORright_(dict.lookupOrDefault<scalar>("leftORright",-1)),
    waveDictName_(dict.lookupOrDefault<word>("waveDict","IHWavesDict")),
    waveType_(dict.lookupOrDefault<word>("waveType","aaa")),
    waveTheory_(dict.lookupOrDefault<word>("waveTheory","aaa")),
    allCheck_(dict.lookupOrDefault<bool>("allCheck",true)),
    waveDir_(dict.lookupOrDefault<scalar>("waveDir",0))
{}


Foam::
IH_Waves_InletVelocityFvPatchVectorField::
IH_Waves_InletVelocityFvPatchVectorField
(
    const IH_Waves_InletVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    wavePeriod_(ptf.wavePeriod_),
    waveHeight_(ptf.waveHeight_),
    waveLength_(ptf.waveLength_),
    waterDepth_(ptf.waterDepth_),
    initialDepth_(ptf.initialDepth_),
    RealwaterDepth_(ptf.RealwaterDepth_),
    wavePhase_(ptf.wavePhase_),
    lambdaStokesV_(ptf.lambdaStokesV_),
    mCnoidal_(ptf.mCnoidal_),
    genAbs_(ptf.genAbs_),
    nPaddles_(ptf.nPaddles_),
    tSmooth_(ptf.tSmooth_),
    leftORright_(ptf.leftORright_),
    waveDictName_(ptf.waveDictName_),
    waveType_(ptf.waveType_),
    waveTheory_(ptf.waveTheory_),
    allCheck_(ptf.allCheck_),
    waveDir_(ptf.waveDir_)
{}


Foam::
IH_Waves_InletVelocityFvPatchVectorField::
IH_Waves_InletVelocityFvPatchVectorField
(
    const IH_Waves_InletVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    wavePeriod_(ptf.wavePeriod_),
    waveHeight_(ptf.waveHeight_),
    waveLength_(ptf.waveLength_),
    waterDepth_(ptf.waterDepth_),
    initialDepth_(ptf.initialDepth_),
    RealwaterDepth_(ptf.RealwaterDepth_),
    wavePhase_(ptf.wavePhase_),
    lambdaStokesV_(ptf.lambdaStokesV_),
    mCnoidal_(ptf.mCnoidal_),
    genAbs_(ptf.genAbs_),
    nPaddles_(ptf.nPaddles_),
    tSmooth_(ptf.tSmooth_),
    leftORright_(ptf.leftORright_),
    waveDictName_(ptf.waveDictName_),
    waveType_(ptf.waveType_),
    waveTheory_(ptf.waveTheory_),
    allCheck_(ptf.allCheck_),
    waveDir_(ptf.waveDir_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::IH_Waves_InletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    } 

    // Auxiliar variables
    scalar auxiliar = 0; 
    scalar auxiliarTotal = 0;
    scalar X0 = 0;
    scalarField patchXsolit;

    // Variables stream function
    scalar celerity = 0;

    // 3D Variables
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
    const fvMesh& mesh = alpha.mesh();
    const word& patchName = this->patch().name();
    const label patchID = mesh.boundaryMesh().findPatchID(patchName);
    const label nF = patch().faceCells().size();
    labelList cellGroup = Foam::labelList(nF, 1);
    const scalarField alphaCell = 
        alpha.boundaryField()[patchID].patchInternalField();
    scalarField patchU = Foam::scalarField(nF, 0.0);
    scalarField patchUABS = Foam::scalarField(nF, 0.0);
    scalarField patchV = Foam::scalarField(nF, 0.0);
    scalarField patchVABS = Foam::scalarField(nF, 0.0);
    scalarField patchW = Foam::scalarField(nF, 0.0);
    const labelList celdas = patch().faceCells();
    const scalarField patchHeight = patch().Cf().component(2);
    const scalar g = 9.81;

    // Calculate Z bounds of the faces
    scalarField zSup, zInf;
    faceBoundsZ( &zSup, &zInf );

    // Waves variables
    scalar waveOmega;
    scalar waveK;
    scalar waveAngle;
    scalar waveKx;
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
	genAbs_ = (IHWavesDict.lookupOrDefault<bool>("genAbs",false));
	waterDepth_ = (IHWavesDict.lookupOrDefault<scalar>("waterDepth",-1));
	initialDepth_ = 
	    (IHWavesDict.lookupOrDefault<scalar>("initialDepth",-1));
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
                cellGroup[patchCells] = i+1; // Group of each face
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
    	genAbs_ = (IHWavesDict.lookupOrDefault<bool>("genAbs",false));
    	wavePhase_ = 
	    (IHWavesDict.lookupOrDefault<scalar>("wavePhase",3.0*PI()/2.0));
	
	if ( waveType_ == "regular" )
	{
		waveLength_ = StokesIFun::waveLength (
		    RealwaterDepth_, wavePeriod_ );
	}
    }

    if ( waveType_ == "regular" )
    {
        waveOmega = (2.0*PI())/wavePeriod_;
        waveK = 2.0*PI()/waveLength_;

        celerity = waveLength_/wavePeriod_;

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
	        << "regular, solitary."<< exit(FatalError);
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

    scalarList measuredLevels (nPaddles_,0.0);
    forAll(measuredLevels, iterMin)
    {
        measuredLevels[iterMin] = RealwaterDepth_;
    }

    scalarList measuredLevelsGENAB = calcWL( alphaCell, cellGroup, zSpan );

    if (initialDepth_ !=-1.0)
    {
    	forAll(measuredLevelsGENAB, iterMin)
    	{
            measuredLevelsGENAB[iterMin] = measuredLevelsGENAB[iterMin] 
	        + initialDepth_;
    	}
    }

    // Define heights as minimum of calculatedLevel and measuredLevels
    scalarList heights (nPaddles_,0.0);
    forAll(heights, iterMin)
    {
        heights[iterMin] = 
	    min(calculatedLevel[iterMin],measuredLevels[iterMin]);
    }
    
    forAll(patchHeight, cellIndex)    
    {
    	#include "velocityProfile.H"
    
    	patchU[cellIndex] = patchU[cellIndex]*alphaCell[cellIndex];
    	patchV[cellIndex] = patchV[cellIndex]*alphaCell[cellIndex];
    	patchW[cellIndex] = patchW[cellIndex]*alphaCell[cellIndex];
    }

    const vectorField n1 = Foam::vectorField(nF, vector(1.0, 0.0, 0.0));
    const vectorField n2 = Foam::vectorField(nF, vector(0.0, 1.0, 0.0));
    const vectorField n3 = Foam::vectorField(nF, vector(0.0, 0.0, 1.0));

    operator == (n1*patchU + n2*patchV + n3*patchW);

    fixedValueFvPatchField<vector>::updateCoeffs();
}


void Foam::IH_Waves_InletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);

    os.writeKeyword("waveDictName") << waveDictName_ 
        << token::END_STATEMENT << nl;

    os.writeKeyword("leftORright") << leftORright_ 
        << token::END_STATEMENT << nl;

    os.writeKeyword("RealwaterDepth") << RealwaterDepth_ 
        << token::END_STATEMENT << nl;

    os.writeKeyword("initialDepth") << initialDepth_ 
        << token::END_STATEMENT << nl;

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       IH_Waves_InletVelocityFvPatchVectorField
   );
}


// ************************************************************************* //
