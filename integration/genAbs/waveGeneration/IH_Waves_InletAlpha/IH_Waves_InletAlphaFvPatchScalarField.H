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
    Foam::IH_Waves_InletAlphaFvPatchScalarField

Description
    Describes a volumetric/mass flow normal Scalar boundary condition by its
    magnitude as an integral over its area.

    The basis of the patch (volumetric or mass) is determined by the
    dimensions of the flux, phi.
    The current density is used to correct the Alpha when applying the
    mass basis.

    Example of the boundary condition specification:
    @verbatim
    inlet
    {
        type            IH_Waves_InletAlphaV2;
        waveDict        IHWavesDict;
    }
    @endverbatim

    \heading Function object usage
    \table
        Property     | Description             | Required    | Default value
        type         | type:    IH_Waves_InletAlphaV2 | yes
        waveDict     | Dictionary where variables for generation/absorption are defined | yes | IHWavesDict  
    \endtable

Note
    - 

SourceFiles
    IH_Waves_InletAlphaFvPatchScalarField.C

\*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*\
   IH-Cantabria 2015 (http://www.ihcantabria.com/en/)
   IHFOAM 2015 (http://ihfoam.ihcantabria.com/) 

   Author(s):  Javier Lopez Lara (jav.lopez@unican.es)
               Gabriel Barajas   (barajasg@unican.es)
\*---------------------------------------------------------------------------*/

#ifndef IH_Waves_InletAlphaFvPatchScalarField_H
#define IH_Waves_InletAlphaFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
               Class IH_Waves_InletAlphaFvPatch Declaration
\*---------------------------------------------------------------------------*/

class IH_Waves_InletAlphaFvPatchScalarField
:
    public fixedValueFvPatchScalarField
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

        //- Number of different paddles (for absorption)
        label nPaddles_;

        //- Smoothing factor for the wave(s) (seconds)
        scalar tSmooth_;

        //- Dictionary name
        word waveDictName_;

        //- Regular or Irregular wave(s)
        word waveType_;

        //- Name of the theory
        word waveTheory_;

        //- BC check and read values just the first time
        bool allCheck_;

        //- Direction of propagation (degrees, from X axis)
        scalar waveDir_;

        //- Numerical Water depth (meters)
	scalar RealwaterDepth_;

        //- Initial Depth on BC for wave generation (meters)
        scalar initialDepth_;

public:

   //- Runtime type information
   TypeName("IH_Waves_InletAlphaV2");


   // Constructors

        //- Construct from patch and internal field
        IH_Waves_InletAlphaFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        IH_Waves_InletAlphaFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  IH_Waves_InletAlphaFvPatchScalarField
        //  onto a new patch
        IH_Waves_InletAlphaFvPatchScalarField
        (
            const IH_Waves_InletAlphaFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        IH_Waves_InletAlphaFvPatchScalarField
        (
            const IH_Waves_InletAlphaFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new IH_Waves_InletAlphaFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        IH_Waves_InletAlphaFvPatchScalarField
        (
            const IH_Waves_InletAlphaFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new IH_Waves_InletAlphaFvPatchScalarField(*this, iF)
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
