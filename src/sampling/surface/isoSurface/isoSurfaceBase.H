/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::isoSurfaceBase

Description
    Low-level components common to various iso-surface algorithms.

SourceFiles
    isoSurfaceBase.C

\*---------------------------------------------------------------------------*/

#ifndef isoSurfaceBase_H
#define isoSurfaceBase_H

#include "scalar.H"
#include "Enum.H"
#include "boundBox.H"
#include "MeshedSurface.H"
#include "MeshedSurfacesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class isoSurfaceBase Declaration
\*---------------------------------------------------------------------------*/

class isoSurfaceBase
:
    public meshedSurface
{
protected:

    // Protected typedefs for convenience
    typedef meshedSurface MeshStorage;

    // Protected Data

        //- Iso value
        const scalar iso_;

        //- Optional bounds
        const boundBox bounds_;

        //- For every face, the original cell in mesh
        labelList meshCells_;


public:

    //- The algorithm types
    enum algorithmType : uint8_t
    {
        ALGO_POINT,
        ALGO_CELL,
        ALGO_TOPO
    };


    //- The filtering (regularization) to apply
    enum class filterType : uint8_t
    {
        NONE = 0,   //!< No filtering
        CELL,       //!< Remove points from pyramid edges
        DIAGCELL,   //!< Remove points from pyramid edges and face-diagonals
        PARTIAL = CELL,     //!< Same as CELL
        FULL = DIAGCELL,    //!< Same as DIAGCELL
    };


    // Public Data

        //- Names for the iso-surface algorithms
        static const Enum<algorithmType> algorithmNames;

        //- Names for the filtering types
        static const Enum<filterType> filterNames;



    // Static Member Functions

        //- Get 'regularise' as bool or enumeration
        static filterType getFilterType
        (
            const dictionary& dict,
            const isoSurfaceBase::filterType deflt
        );


    // Constructors

        //- Construct with iso value
        explicit isoSurfaceBase
        (
            const scalar iso,
            const boundBox& bounds = boundBox::invertedBox
        );


    // Member Functions

        //- The iso-value associated with the surface
        inline scalar isoValue() const
        {
            return iso_;
        }

        //- For each face, the original cell in mesh
        const labelList& meshCells() const
        {
            return meshCells_;
        }

        //- For each face, the original cell in mesh
        labelList& meshCells()
        {
            return meshCells_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
