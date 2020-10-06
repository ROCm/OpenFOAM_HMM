/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

Namespace
    Foam::blockEdges

Description
    A namespace for various blockEdge types.

Class
    Foam::blockEdge

Description
    Define a curved edge that is parameterized for 0<lambda<1
    between the start/end points.

SourceFiles
    blockEdge.C

\*---------------------------------------------------------------------------*/

#ifndef blockEdge_H
#define blockEdge_H

#include "searchableSurfaces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class blockEdge;
Ostream& operator<<(Ostream& os, const blockEdge& e);

/*---------------------------------------------------------------------------*\
                         Class blockEdge Declaration
\*---------------------------------------------------------------------------*/

class blockEdge
{
protected:

    // Protected Data

        //- The referenced point field
        const pointField& points_;

        //- Index of the start point
        const label start_;

        //- Index of the end point
        const label end_;


    // Protected Member Functions

        //- Return a complete point field by appending the start/end points
        //- to the given list
        //  \deprecated(2020-10) use polyLine::concat
        static pointField appendEndPoints
        (
            const pointField& points,   //!< Referenced point field
            const label start,      //!< Start point in referenced point field
            const label end,        //!< End point in referenced point field
            const pointField& intermediate  //!< Intermediate points (knots)
        );


public:

    //- Runtime type information
    TypeName("blockEdge");

    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            blockEdge,
            Istream,
            (
                const dictionary& dict,
                const label index,
                const searchableSurfaces& geometry,
                const pointField& points,
                Istream& is
            ),
            (dict, index, geometry, points, is)
        );


    // Constructors

        //- Construct from components
        blockEdge
        (
            const pointField& points,   //!< Referenced point field
            const label start,      //!< Start point in referenced point field
            const label end         //!< End point in referenced point field
        );

        //- Construct from Istream and point field.
        blockEdge
        (
            const dictionary& dict,
            const label index,
            const pointField& points,   //!< Referenced point field
            Istream& is
        );

        //- Clone function
        virtual autoPtr<blockEdge> clone() const;

        //- New function which constructs and returns pointer to a blockEdge
        static autoPtr<blockEdge> New
        (
            const dictionary& dict,
            const label index,
            const searchableSurfaces& geometry,
            const pointField& points,
            Istream& is
        );

        //- Class used for the read-construction of
        //  PtrLists of blockEdge
        class iNew
        {
            const dictionary& dict_;
            const searchableSurfaces& geometry_;
            const pointField& points_;
            mutable label index_;

        public:

            iNew
            (
                const dictionary& dict,
                const searchableSurfaces& geometry,
                const pointField& points
            )
            :
                dict_(dict),
                geometry_(geometry),
                points_(points),
                index_(0)
            {}

            autoPtr<blockEdge> operator()(Istream& is) const
            {
                return blockEdge::New(dict_, index_++, geometry_, points_, is);
            }
        };


    //- Destructor
    virtual ~blockEdge() = default;


    // Member Functions

        //- Index of start point
        inline label start() const;

        //- Index of end point
        inline label end() const;

        //- Compare the given start/end points with this block edge
        //  Return:
        //  -  0: different
        //  - +1: identical
        //  - -1: same edge, but different orientation
        inline int compare(const blockEdge& e) const;

        //- Compare the given start/end points with this block edge
        //  Return:
        //  -  0: different
        //  - +1: identical
        //  - -1: same edge, but different orientation
        inline int compare(const edge& e) const;

        //- Compare the given start/end points with this block edge
        //  Return:
        //  -  0: different
        //  - +1: identical
        //  - -1: same edge, but different orientation
        inline int compare(const label start, const label end) const;

        //- The point position corresponding to the curve parameter
        //  0 <= lambda <= 1
        virtual point position(const scalar lambda) const = 0;

        //- The point positions corresponding to the curve parameters
        //  0 <= lambda <= 1
        virtual tmp<pointField> position(const scalarList& lambdas) const;

        //- The length of the curve
        virtual scalar length() const = 0;

        //- Write edge with variable back-substitution
        void write(Ostream& os, const dictionary& dict) const;


    // Ostream Operator

        friend Ostream& operator<<(Ostream& os, const blockEdge& e);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "blockEdgeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
