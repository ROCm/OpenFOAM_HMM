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

Class
    Foam::waveMakerPistonPointPatchVectorField

Description
    Foam::waveMakerPistonPointPatchVectorField

SourceFiles
    waveMakerPistonPointPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef waveMakerPistonPointPatchVectorField_H
#define waveMakerPistonPointPatchVectorField_H

#include "fixedValuePointPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class waveMakerPistonPointPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class waveMakerPistonPointPatchVectorField
:
    public fixedValuePointPatchField<vector>
{
    // Private data

	//inital water depth
        scalar initialDepth_;

	//wave period
        scalar wavePeriod_;

	//wave heigth
        scalar waveHeigth_;

	//wavelength
        scalar waveLength_;

	//wave phase
        scalar wavePhase_;

	//wave number
        scalar waveNumber_;

	//ramp time
        scalar rampTime_;

	// gravity vector
        vector g_;

	//on/off second order generation
        scalar secondOrder_;

    // Protected Member Functions

        // Dispersion equation
            virtual scalar waveLength(const scalar h, const scalar T);

        //- Return the time scaling coefficient
            virtual scalar timeCoeff(const scalar t) const;

public:

    //- Runtime type information
    TypeName("waveMakerPiston");


    // Constructors

        //- Construct from patch and internal field
        waveMakerPistonPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        waveMakerPistonPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<vector> onto a new patch
        waveMakerPistonPointPatchVectorField
        (
            const waveMakerPistonPointPatchVectorField&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<vector>> clone() const
        {
            return autoPtr<pointPatchField<vector>>
            (
                new waveMakerPistonPointPatchVectorField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        waveMakerPistonPointPatchVectorField
        (
            const waveMakerPistonPointPatchVectorField&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<vector>> clone
        (
            const DimensionedField<vector, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchField<vector>>
            (
                new waveMakerPistonPointPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
