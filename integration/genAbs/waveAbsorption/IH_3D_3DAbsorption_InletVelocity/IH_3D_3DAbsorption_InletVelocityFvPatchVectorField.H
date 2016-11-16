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
    Foam::IH_3D_3DAbsorption_InletVelocityFvPatchVectorField

Description
    Wave absorption boundary condition based on shallow water theory and on a
    3D approach. Works both in 2D and 3D and for waves out of the shallow water
    regime.

    Example of the boundary condition specification:
    @verbatim
    outlet
    {
        type            IH_3D_3DAbsorption_InletVelocityV2;
        nPaddles        1;  
        value           uniform (0 0 0);
	leftORright	-1.0;
    }
    @endverbatim

    \heading Function object usage
    \table
        Property     | Description             | Required    | Default value
        type         | type:    IH_3D_3DAbsorption_InletVelocityV2 | yes
        nPaddles     | number of wavepaddles for absorption | yes | 1  
	leftORright  | Define location of Boundary condition: Left(1) or Right (-1) | yes | -1
    \endtable

Note
    -

SourceFiles
    IH_3D_3DAbsorption_InletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*\
   IH-Cantabria 2015 (http://www.ihcantabria.com/en/)
   IHFOAM 2015 (http://ihfoam.ihcantabria.com/) 

   Author(s):  Javier Lopez Lara (jav.lopez@unican.es)
               Gabriel Barajas   (barajasg@unican.es)
\*---------------------------------------------------------------------------*/

#ifndef IH_3D_3DAbsorption_InletVelocityFvPatchVectorField_H
#define IH_3D_3DAbsorption_InletVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
               Class IH_3D_3DAbsorption_InletVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class IH_3D_3DAbsorption_InletVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Number of paddles
        label nPaddles_;
 
        //Aborption: Left(1) or Right(-1)
	scalar leftORright_;

        //- Theoretical Water depth (meters)
        scalar waterDepth_;        

        //- Numerical Water depth (meters)
        scalar RealwaterDepth_;

        //- BC check and read values just the first time
        bool allCheck_;

        //- Initial Depth on BC for absorption (meters)
        scalar initialDepthABS_;

        //- Dictionary name
        word waveDictName_;

public:

   //- Runtime type information
   TypeName("IH_3D_3DAbsorption_InletVelocityV2");


   // Constructors

        //- Construct from patch and internal field
        IH_3D_3DAbsorption_InletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        IH_3D_3DAbsorption_InletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  IH_3D_3DAbsorption_InletVelocityFvPatchVectorField
        //  onto a new patch
        IH_3D_3DAbsorption_InletVelocityFvPatchVectorField
        (
            const IH_3D_3DAbsorption_InletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        IH_3D_3DAbsorption_InletVelocityFvPatchVectorField
        (
            const IH_3D_3DAbsorption_InletVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new IH_3D_3DAbsorption_InletVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        IH_3D_3DAbsorption_InletVelocityFvPatchVectorField
        (
            const IH_3D_3DAbsorption_InletVelocityFvPatchVectorField&,
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
                new IH_3D_3DAbsorption_InletVelocityFvPatchVectorField
                    (*this, iF)
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
