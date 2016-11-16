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

Class
    Foam::IH_Waves_InletVelocityFvPatchVectorField

Description
    Describes a volumetric/mass flow normal vector boundary condition by its
    magnitude as an integral over its area.

    The basis of the patch (volumetric or mass) is determined by the
    dimensions of the flux, phi.
    The current density is used to correct the velocity when applying the
    mass basis.

    Example of the boundary condition specification:
    @verbatim
    inlet
    {
        type            IH_Waves_InletVelocityV2;
        waveDict        IHWavesDict;
        value           uniform (0 0 0);
	leftORright	1.0;

    }
    @endverbatim

    \heading Function object usage
    \table
        Property     | Description             | Required    | Default value
        type         | type:    IH_Waves_InletVelocityV2 | yes
        waveDict     | Dictionary where variables for generation/absorption are defined | yes | IHWavesDict 
	leftORright  | Define location of Boundary condition: Left(1) or Right (-1) | yes | -1 
    \endtable

Note
    - The value is positive inwards
    - May not work correctly for transonic inlets
    - Strange behaviour with potentialFoam since the U equation is not solved

SourceFiles
    IH_Waves_InletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*\
   IH-Cantabria 2015 (http://www.ihcantabria.com/en/)
   IHFOAM 2015 (http://ihfoam.ihcantabria.com/) 

   Author(s):  Javier Lopez Lara (jav.lopez@unican.es)
               Gabriel Barajas   (barajasg@unican.es)
\*---------------------------------------------------------------------------*/

#ifndef IH_Waves_InletVelocityFvPatchVectorField_H
#define IH_Waves_InletVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
               Class IH_Waves_InletVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class IH_Waves_InletVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Wave period (seconds)
        scalar wavePeriod_;

        //- Wave height (meters)
        scalar waveHeight_;

        //- Wave length (meters)
        scalar waveLength_;

        //- Water depth (meters)
        scalar waterDepth_;

        //- Wave phase (radians)
        scalar wavePhase_;

        //- Lambda StokesV parameter
        scalar lambdaStokesV_;

        //- m Cnoidal parameter
        scalar mCnoidal_;

        //- Generation + Absorption at the same time in the boundary (on=1, off=0)
        bool genAbs_;

        //- Number of different paddles (for absorption)
        label nPaddles_;

        //- Smoothing factor for the wave(s) (seconds)
        scalar tSmooth_;

        //Aborption Left(1) or Right(-1)
	scalar leftORright_;

        //- Dictionary name
        word waveDictName_;

        //- Regular or Irregular
        word waveType_;

        //- Name of the theory
        word waveTheory_;

        //- BC check and read values just the first time
        bool allCheck_;

        //- Direction of propagation (degrees, from X axis)
        scalar waveDir_;

        //- Initial water depth (referece)
	scalar RealwaterDepth_;

        //- Water depth (meters)
        scalar initialDepth_;

public:

   //- Runtime type information
   TypeName("IH_Waves_InletVelocityV2");


   // Constructors

        //- Construct from patch and internal field
        IH_Waves_InletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        IH_Waves_InletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  IH_Waves_InletVelocityFvPatchVectorField
        //  onto a new patch
        IH_Waves_InletVelocityFvPatchVectorField
        (
            const IH_Waves_InletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        IH_Waves_InletVelocityFvPatchVectorField
        (
            const IH_Waves_InletVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new IH_Waves_InletVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        IH_Waves_InletVelocityFvPatchVectorField
        (
            const IH_Waves_InletVelocityFvPatchVectorField&,
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
                new IH_Waves_InletVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;

    // Other common member functions
        #include "memberFun.H"

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
