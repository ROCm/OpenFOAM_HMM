/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd.
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
    Foam::fv::cellSetOption

Description
    Intermediate abstract class for handling
    cell-set options for the derived fvOptions.

Usage
    Minimal example by using \c constant/fvOptions:
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

            // when selectionMode=cellSet
            cellSet           <name>;

            // when selectionMode=cellZone
            cellZone          <name>;

            // when selectionMode=points
            points            (<point1> <point2> ... <pointN>);

            // when selectionMode=geometric
            selection
            {
                box1
                {
                    action  use;
                    source  box;
                    min     (-0.1 -0.01 -0.1);
                    max     (0.1 0.30 0.1);
                }
                ball
                {
                    action  add;
                    source  sphere;
                    origin  (-0.1 -0.01 -0.1);
                    radius  0.25;
                }
                ...
            }

        // Mandatory/Optional (derived) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property      | Description                        | Type  | Reqd | Dflt
      selectionMode | Mode of cell selection - see below | word  | yes  | -
      timeStart     | Start time of fvOption             | scalar | no  | -1
      duration      | Duration of fvOption execution <!--
                          --> starting from timeStart    | scalar | cndtnl | 0
      cellSet       | Name of operand cellSet            | word  | cndtnl | -
      cellZone      | Name of operand cellZone           | word  | cndtnl | -
      points        | Set of points in global coordinate <!--
                            --> system              | vectorList | cndtnl | -
      selection     | Dictionary of geometric selections | dict | cndtnl | -
    \endtable

    Options for the \c selectionMode entry:
    \verbatim
      all       | Use all cells in the computational domain
      cellZone  | Use specified cellZone
      cellSet   | Use specified cellSet
      points    | Use cells containing a given set of points
      geometric | Select cells based on topoSetCellSource actions
    \endverbatim

    The inherited entries are elaborated in:
     - \link fvOption.H \endlink

    The geometric selection uses topoSetCellSource to select cells.
    Any searchableSurface selections must describe a closed volume.
    Ie, its hasVolumeType() method must be true.

Note
  - Source/sink options are to be added to the right-hand side of equations.

See also
    Foam::cellBitSet::select

SourceFiles
    cellSetOption.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_fv_cellSetOption_H
#define Foam_fv_cellSetOption_H

#include "fvOption.H"
#include "fvMesh.H"
#include "dictionary.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                         Class cellSetOption Declaration
\*---------------------------------------------------------------------------*/

class cellSetOption
:
    public fv::option
{
public:

    // Public Data

        //- Enumeration for selection mode types
        enum selectionModeType
        {
            smAll,          //!< "all" cells
            smCellSet,      //!< "cellSet"
            smCellZone,     //!< "cellZone"
            smPoints,       //!< "points"
            smGeometric     //!< "geometric"
        };

        //- List of selection mode type names
        static const Enum<selectionModeType> selectionModeTypeNames_;


protected:

    // Protected Data

        //- Start time of fvOption
        scalar timeStart_;

        //- Duration of fvOption execution starting from timeStart
        scalar duration_;

        //- Cell selection mode
        selectionModeType selectionMode_;

        //- Name of set/zone for "cellSet" and "cellZone" selectionMode
        wordRe zoneName_;

        //- List of points for "points" selectionMode
        List<point> points_;

        //- Dictionary entries for "geometric" (topoSetCellSource) selection
        dictionary geometricSelection_;

        //- Set of cells to apply source to
        labelList cells_;

        //- Sum of cell volumes
        scalar V_;


    // Protected Functions

        //- Set cell selection name or points selection from dictionary input
        void setSelection(const dictionary& dict);

        //- Set the cell selection based on user input selection mode
        void setCellSelection();

        //- Recalculate the volume
        void setVol();


public:

    //- Runtime type information
    TypeName("cellSetOption");


    // Constructors

        //- Construct from components
        cellSetOption
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~cellSetOption() = default;


    // Member Functions

        // Access

            //- Return const access to the time start
            inline scalar timeStart() const noexcept;

            //- Return const access to the duration
            inline scalar duration() const noexcept;

            //- True if within time limits
            inline bool inTimeLimits(const scalar timeValue) const;

            //- Return the cell selection mode
            inline selectionModeType selectionMode() const noexcept;

            //- True if sub-selection should be used
            inline bool useSubMesh() const noexcept;

            //- Return const access to the name of cell set/zone
            //- for "cellSet" / "cellZone" selection modes
            const wordRe& zoneName() const noexcept { return zoneName_; }

            //- Return const access to the total cell volume
            inline scalar V() const noexcept;

            //- Return const access to the cell selection
            inline const labelList& cells() const noexcept;


        // Edit

            //- Adjust the time start, return the old value
            inline scalar timeStart(scalar val) noexcept;

            //- Adjust the duration, return the old value
            inline scalar duration(scalar val) noexcept;


        // Checks

            //- Is the source active?
            virtual bool isActive();


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);


    // Housekeeping

        //- The name of the cell set/zone [as a word]
        //- for "cellSet" / "cellZone" selection modes)
        const word& cellSetName() const noexcept { return zoneName_; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cellSetOptionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
