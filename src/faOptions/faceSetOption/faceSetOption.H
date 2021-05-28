/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
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
    Foam::fa::faceSetOption

Description
    Intermediate abstract class for handling
    face-set options for the derived faOptions.

Usage
    Minimal example by using \c constant/faOptions:
    \verbatim
    <userDefinedName1>
    {
        // Mandatory/Optional (inherited) entries
        ...

        // Mandatory entries (unmodifiable)
        selectionMode     all;

        // Optional entries (runtime modifiable)
        timeStart         1.0;

        // Conditional mandatory entries (runtime modifiable)

            // when timeStart entry is present
            duration          1.4;

            // when selectionMode=volFaceZone
            faceZone          <faceZoneName>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property      | Description                        | Type  | Reqd | Dflt
      selectionMode | Mode of face selection - see below | word  | yes  | -
      timeStart     | Start time of faOption             | scalar | no  | -1
      duration      | Duration of faOption execution <!--
                          --> starting from timeStart    | scalar | cndtnl | 0
      faceZone      | Name of operand faceZone           | word  | cndtnl | -
    \endtable

    Options for the \c selectionMode entry:
    \verbatim
      all       | Use all faces in the computational domain
      faceZone  | Use a given faceZone
    \endverbatim

    The inherited entries are elaborated in:
     - \link faOption.H \endlink

Note
  - Source/sink options are to be added to the right-hand side of equations.

SourceFiles
    faceSetOption.C

\*---------------------------------------------------------------------------*/

#ifndef faceSetOption_H
#define faceSetOption_H

#include "faOption.H"
#include "faceSet.H"
#include "faMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fa
{

/*---------------------------------------------------------------------------*\
                         Class faceSetOption Declaration
\*---------------------------------------------------------------------------*/

class faceSetOption
:
    public fa::option
{
public:

    // Public Enumeration

        //- Enumeration for selection mode types
        enum selectionModeType
        {
            smAll,
            smVolFaceZone
        };

        //- List of selection mode type names
        static const Enum<selectionModeType> selectionModeTypeNames_;


protected:

    // Protected Data

        //- Time start
        scalar timeStart_;

        //- Duration
        scalar duration_;

        //- Face selection mode
        selectionModeType selectionMode_;

        //- Name of zone for "faceZone" selectionMode
        word faceSetName_;

        //- Set of faces to apply source to
        labelList faces_;

        //- Sum of face area
        scalar A_;


    // Protected Functions

        //- Set the face selection
        void setSelection(const dictionary& dict);

        //- Set the face set based on the user input selection mode
        void setFaceSet();

        //- Recalculate the area
        void setArea();


public:

    //- Runtime type information
    TypeName("faceSetOption");


    // Constructors

        //- Construct from components
        faceSetOption
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvPatch& patch
        );


    //- Destructor
    virtual ~faceSetOption() = default;


    // Member Functions

        // Access

            //- Return const access to the time start
            inline scalar timeStart() const;

            //- Return const access to the duration
            inline scalar duration() const;

            //- Return true if within time limits
            inline bool inTimeLimits(const scalar time) const;

            //- Return const access to the face selection mode
            inline const selectionModeType& selectionMode() const;

            //- Return const access to the name of face set for "faceZone"
            //- selectionMode
            inline const word& faceSetName() const;

            //- Return const access to the total face area
            inline scalar A() const;

            //- Return const access to the face set
            inline const labelList& faces() const;


        // Edit

            //- Return access to the time start
            inline scalar& timeStart();

            //- Return access to the duration
            inline scalar& duration();


        // Checks

            //- Is the source active?
            virtual bool isActive();


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "faceSetOptionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
