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

#include "IH_3D_3DAbsorption_InletVelocityFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
IH_3D_3DAbsorption_InletVelocityFvPatchVectorField::
IH_3D_3DAbsorption_InletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    nPaddles_(1),
    leftORright_(-1),
    waterDepth_(-1),
    initialDepthABS_(-1),
    RealwaterDepth_(-1),
    allCheck_(true),
    waveDictName_("IHWavesDict")
{}


Foam::
IH_3D_3DAbsorption_InletVelocityFvPatchVectorField::
IH_3D_3DAbsorption_InletVelocityFvPatchVectorField
(
    const IH_3D_3DAbsorption_InletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    nPaddles_(ptf.nPaddles_),
    leftORright_(ptf.leftORright_),    
    waterDepth_(ptf.waterDepth_),
    initialDepthABS_(ptf.initialDepthABS_),
    RealwaterDepth_(ptf.RealwaterDepth_),
    allCheck_(ptf.allCheck_),
    waveDictName_(ptf.waveDictName_)
{}


Foam::
IH_3D_3DAbsorption_InletVelocityFvPatchVectorField::
IH_3D_3DAbsorption_InletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    nPaddles_(dict.lookupOrDefault<label>("nPaddles", 1)),
    leftORright_(dict.lookupOrDefault<scalar>("leftORright",-1)),
    waterDepth_(dict.lookupOrDefault<scalar>("waterDepth",-1 )),
    initialDepthABS_(dict.lookupOrDefault<scalar>("initialDepthABS",-1)),
    RealwaterDepth_(dict.lookupOrDefault<scalar>("RealwaterDepth",-1)),
    allCheck_(dict.lookupOrDefault<bool>("allCheck", true )),
    waveDictName_(dict.lookupOrDefault<word>("waveDict","IHWavesDict"))
{}


Foam::
IH_3D_3DAbsorption_InletVelocityFvPatchVectorField::
IH_3D_3DAbsorption_InletVelocityFvPatchVectorField
(
    const IH_3D_3DAbsorption_InletVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    nPaddles_(ptf.nPaddles_),
    leftORright_(ptf.leftORright_),
    waterDepth_(ptf.waterDepth_),
    initialDepthABS_(ptf.initialDepthABS_),
    RealwaterDepth_(ptf.RealwaterDepth_),
    allCheck_(ptf.allCheck_),
    waveDictName_(ptf.waveDictName_)
{}


Foam::
IH_3D_3DAbsorption_InletVelocityFvPatchVectorField::
IH_3D_3DAbsorption_InletVelocityFvPatchVectorField
(
    const IH_3D_3DAbsorption_InletVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    nPaddles_(ptf.nPaddles_),
    leftORright_(ptf.leftORright_),
    waterDepth_(ptf.waterDepth_),
    initialDepthABS_(ptf.initialDepthABS_),
    RealwaterDepth_(ptf.RealwaterDepth_),
    allCheck_(ptf.allCheck_),
    waveDictName_(ptf.waveDictName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::IH_3D_3DAbsorption_InletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Variables & constants
    const vector cMin = gMin(patch().patch().localPoints());
    const vector cMax = gMax(patch().patch().localPoints());
    const vector cSpan = cMax - cMin;
    const scalar zSpan = cSpan[2];
    scalar dMin = 0.0;
    scalar dSpan = 0.0;
    const scalarField patchD = patchDirection( cSpan, &dMin, &dSpan );
    const volScalarField& alpha = 
        db().lookupObject<volScalarField>(alphaName());
    const volVectorField& U = db().lookupObject<volVectorField>("U");
    const fvMesh& mesh = alpha.mesh();
    const word& patchName = this->patch().name();
    const label patchID = mesh.boundaryMesh().findPatchID(patchName);
    const label nF = patch().faceCells().size();
    const scalarField alphaCell = 
        alpha.boundaryField()[patchID].patchInternalField();
    const vectorField UCell = U.boundaryField()[patchID].patchInternalField();
    scalarField patchUc = Foam::scalarField(nF, 0.0);
    scalarField patchVc = Foam::scalarField(nF, 0.0);
    scalarField patchWc = Foam::scalarField(nF, 0.0);
    const scalarField patchHeight = patch().Cf().component(2);
    const scalar g = 9.81;

    // Calculate Z bounds of the faces
    scalarField zSup, zInf;
    faceBoundsZ( &zSup, &zInf );

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
        waterDepth_ = (IHWavesDict.lookupOrDefault<scalar>("waterDepth",-1.0));

	initialDepthABS_ = 
	    (IHWavesDict.lookupOrDefault<scalar>("initialDepthABS",-1.0));

        // Check if the value of nPaddles is correct for the number of columns
        if (nPaddles_ < 1)
        {
            FatalError << "Check nPaddles value." << exit(FatalError);
        }

        if ( nPaddles_ > 1 )
        {
            nPaddles_ = decreaseNPaddles( nPaddles_, patchD, dMin, dSpan );
            reduce(nPaddles_, minOp<label>());
        }
    }

    // Grouping part
    labelList faceGroup = Foam::labelList(nF, 0);
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
                faceGroup[patchCells] = i+1;
                continue;
            }
        }      
    }

    if (allCheck_)
    {
	if (RealwaterDepth_ == -1.0)
	{
	    if (waterDepth_ == -1.0)
	    {
	        RealwaterDepth_ = 
		calcWL( alphaCell, faceGroup, zSpan )[0];

		if (initialDepthABS_ !=-1.0)
		{
		    RealwaterDepth_ = RealwaterDepth_
		        + initialDepthABS_;
		}
	    }
	    else if ( waterDepth_ != -1.0 )
	    {
	        RealwaterDepth_ = waterDepth_;
	    }
	}

        allCheck_ = false;
    }

    // Calculate water measured levels
    scalarList measuredLevelsABS = calcWL( alphaCell, faceGroup, zSpan );

    if (initialDepthABS_ !=-1.0)
    {
    	forAll(measuredLevelsABS, iterMin)
    	{
            measuredLevelsABS[iterMin] = measuredLevelsABS[iterMin] 
	        + initialDepthABS_;
    	}
    }

    forAll(patchHeight, cellIndex)    
    {
    	if ( zInf[cellIndex] >= measuredLevelsABS[faceGroup[cellIndex]-1] )
	{
	    patchUc[cellIndex] = 0.0;
	    patchVc[cellIndex] = 0.0;
             patchWc[cellIndex] = 0.0;
	}
	else 
	{
	    patchUc[cellIndex] = 
	        ( -measuredLevelsABS[faceGroup[cellIndex]-1]
		+ RealwaterDepth_ )
		* sqrt(g/measuredLevelsABS[faceGroup[cellIndex]-1]);

	    patchUc[cellIndex] = leftORright_ * patchUc[cellIndex];
            patchVc[cellIndex] = 0.0;
            patchWc[cellIndex] = 0.0;
	}
    }

    const vectorField n1 = Foam::vectorField(nF, vector(1.0, 0.0, 0.0));
    const vectorField n2 = Foam::vectorField(nF, vector(0.0, 1.0, 0.0));
    const vectorField n3 = Foam::vectorField(nF, vector(0.0, 0.0, 1.0));

    operator == (n1*patchUc + n2*patchVc + n3*patchWc);

    fixedValueFvPatchField<vector>::updateCoeffs();
}


void Foam::IH_3D_3DAbsorption_InletVelocityFvPatchVectorField::
write(Ostream& os) const
{
    fvPatchField<vector>::write(os);

    os.writeKeyword("waveDictName") << waveDictName_ << 
        token::END_STATEMENT << nl;

    os.writeKeyword("leftORright") << leftORright_ << 
        token::END_STATEMENT << nl;

    os.writeKeyword("nPaddles") << nPaddles_ << token::END_STATEMENT << nl;

    os.writeKeyword("RealwaterDepth") << RealwaterDepth_ << 
        token::END_STATEMENT << nl;

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       IH_3D_3DAbsorption_InletVelocityFvPatchVectorField
   );
}


// ************************************************************************* //
