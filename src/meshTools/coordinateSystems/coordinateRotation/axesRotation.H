/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::axesRotation

Description
    A coordinate rotation specified using global axes

    The rotation is defined by a combination of vectors (e1/e2), (e2/e3)
    or (e3/e1). Any nonorthogonality will be absorbed into the second
    vector. In terms of cylindrical coordinates, the 'axis' would
    correspond to the \a z-axis and the 'direction' to the \a r-axis.

    \verbatim
    axesRotation
    {
        type        axesRotation;
        e1          (1 0 0);
        e2          (0 1 0);
    }
    \endverbatim

SourceFiles
    axesRotation.C

\*---------------------------------------------------------------------------*/

#ifndef axesRotation_H
#define axesRotation_H

#include "vector.H"
#include "coordinateRotation.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class axesRotation Declaration
\*---------------------------------------------------------------------------*/

class axesRotation
:
    public coordinateRotation
{
public:

    //- The order/combination of local axes for the axes-rotation definition
    //  Note that these follow the right-hand rule.
    enum axisOrder
    {
        E1_E2,  //!< The axis is X-dominant, the direction is Y-dominant
        E2_E3,  //!< The axis is Y-dominant, the direction is Z-dominant
        E3_E1   //!< The axis is Z-dominant, the direction is X-dominant
    };


private:

    // Private data

        //- Local-to-Global transformation tensor
        tensor R_;

        //- Global-to-Local transformation tensor
        tensor Rtr_;


public:

    //- Runtime type information
    TypeName("axesRotation");

    // Constructors

        //- Construct null
        axesRotation();

        //- Construct as copy
        axesRotation(const axesRotation& r);

        //- Construct from local to global rotation matrix
        explicit axesRotation(const tensor& R);

        //- Construct from two axes (axis and direction)
        axesRotation
        (
            const vector& axis,
            const vector& dir,
            const axisOrder& order = E3_E1
        );

        //- Construct from a single axis using a best-guess for the second axis
        //  For the best-guess, the largest component value and sign of the
        //  axis determines the direction orientation.
        explicit axesRotation(const vector& axis);

        //- Construct from dictionary
        explicit axesRotation(const dictionary& dict);

        //- Construct from dictionary and a registry (typically a mesh)
        axesRotation
        (
            const dictionary& dict,
            const objectRegistry& unused
        );

        //- Return clone
        autoPtr<coordinateRotation> clone() const
        {
            return autoPtr<coordinateRotation>(new axesRotation(*this));
        }


    //- Destructor
    virtual ~axesRotation()
    {}


    // Member Functions

        //- Reset rotation to an identity rotation
        virtual void clear()
        {
            R_ = sphericalTensor::I;
            Rtr_ = sphericalTensor::I;
        }

        //- Set the transformation tensors from two axes (axis and direction)
        void setTransform
        (
            const vector& axis1,
            const vector& axis2,
            const axisOrder& order = E3_E1
        );

        //- Update the rotation for a list of cells
        virtual void updateCells(const polyMesh&, const labelList&)
        {}

        //- Return local-to-global transformation tensor
        virtual const tensor& R() const
        {
            return R_;
        }

        //- Return global-to-local transformation tensor
        virtual const tensor& Rtr() const
        {
            return Rtr_;
        }

        //- Return local Cartesian x-axis in global coordinates
        virtual const vector e1() const
        {
            return Rtr_.x();
        }

        //- Return local Cartesian y-axis in global coordinates
        virtual const vector e2() const
        {
            return Rtr_.y();
        }

        //- Return local Cartesian z-axis in global coordinates
        virtual const vector e3() const
        {
            return Rtr_.z();
        }

        //- Return transformation tensor field
        virtual const tensorField& Tr() const;

        //- Transform vectorField using transformation tensor field
        virtual tmp<vectorField> transform(const vectorField& st) const;

        //- Transform vector using transformation tensor
        virtual vector transform(const vector& st) const;

        //- Inverse transform vectorField using transformation tensor field
        virtual tmp<vectorField> invTransform(const vectorField& st) const;

        //- Inverse transform vector using transformation tensor
        virtual vector invTransform(const vector& st) const;

        //- Transform tensor field using transformation tensorField
        virtual tmp<tensorField> transformTensor(const tensorField& st) const;

        //- Transform tensor using transformation tensorField
        virtual tensor transformTensor(const tensor& st) const;

        //- Transform tensor sub-field using transformation tensorField
        virtual tmp<tensorField> transformTensor
        (
            const tensorField& st,
            const labelList& cellMap
        ) const;

        //- Transform vectorField using transformation tensorField and return
        //  symmetric tensorField
        virtual tmp<symmTensorField> transformVector
        (
            const vectorField& st
        ) const;

        //- Transform vector using transformation tensor and return
        //  symmetric tensor
        virtual symmTensor transformVector(const vector& st) const;


    // Member Operators

        //- Assign from dictionary
        void operator=(const dictionary& dict);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
