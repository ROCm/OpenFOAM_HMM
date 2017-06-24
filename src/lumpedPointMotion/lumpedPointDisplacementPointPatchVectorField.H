/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::lumpedPointDisplacementPointPatchVectorField

Description

    Interpolates pre-specified motion with motion specified as
    pointVectorFields.

    This is the point-patch responsible for managing the force
    integration on a 'lumped-point' basis, waiting for the external
    application, reading back the response from the external program
    and updating the locations of the associated patch points
    accordingly.

    The internal patch type name is 'lumpedPointDisplacement'.

    \heading Patch usage

    Example:
    \verbatim
    walls
    {
        type                lumpedPointDisplacement;
        value               uniform (0 0 0);
        fieldName           wantedDisplacement;
        interpolationScheme linear;
    }
    \endverbatim

    This will scan the case for \a wantedDisplacement pointVectorFields
    and interpolate those in time (using \c linear interpolation) to
    obtain the current displacement.
    The advantage of specifying displacement in this way is that it
    automatically works through decomposePar.

SourceFiles
    lumpedPointDisplacementPointPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef lumpedPointDisplacementPointPatchVectorField_H
#define lumpedPointDisplacementPointPatchVectorField_H

#include "fixedValuePointPatchField.H"
#include "lumpedPointMovement.H"
#include "lumpedPointState.H"
#include "lumpedPointIOMovement.H"
#include "labelList.H"
#include "tmp.H"
#include "pointField.H"
#include "pointFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class interpolationWeights;

/*---------------------------------------------------------------------------*\
        Class lumpedPointDisplacementPointPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class lumpedPointDisplacementPointPatchVectorField
:
    public fixedValuePointPatchField<vector>
{
    // Private data

    //- Convenience typedefs
    typedef lumpedPointDisplacementPointPatchVectorField patchType;
    typedef DimensionedField<vector, pointMesh> fieldType;

    // Private Member Functions

protected:

    //- The starting locations (obtained from the motionSolver).
    const pointField& points0() const;

    //- The auto-vivifying singleton for movement.
    const lumpedPointMovement& movement() const;


public:

    //- Runtime type information
    TypeName("lumpedPointDisplacement");


    // Constructors

        //- Construct from patch and internal field
        lumpedPointDisplacementPointPatchVectorField
        (
            const pointPatch& p,
            const DimensionedField<vector, pointMesh>& iF
        );

        //- Construct from patch, internal field and dictionary
        lumpedPointDisplacementPointPatchVectorField
        (
            const pointPatch& p,
            const DimensionedField<vector, pointMesh>& iF,
            const dictionary& dict
        );

        //- Construct by mapping given patchField<vector> onto a new patch
        lumpedPointDisplacementPointPatchVectorField
        (
            const lumpedPointDisplacementPointPatchVectorField& pf,
            const pointPatch& p,
            const DimensionedField<vector, pointMesh>& iF,
            const pointPatchFieldMapper& mapper
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<vector>> clone() const
        {
            return autoPtr<pointPatchField<vector>>
            (
                new lumpedPointDisplacementPointPatchVectorField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        lumpedPointDisplacementPointPatchVectorField
        (
            const lumpedPointDisplacementPointPatchVectorField& pf,
            const DimensionedField<vector, pointMesh>& iF
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<vector>> clone
        (
            const DimensionedField<vector, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchField<vector>>
            (
                new lumpedPointDisplacementPointPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }

    //- Destructor
    virtual ~lumpedPointDisplacementPointPatchVectorField();


    // Member functions

        //- The ids for all patches of this type
        static labelList patchIds(const pointVectorField& pvf);

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
