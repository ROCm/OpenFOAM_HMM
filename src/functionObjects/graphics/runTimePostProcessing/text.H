/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::functionObjects::runTimePostPro::text

Description
    Example of text object specification:
    \verbatim
    text1
    {
        string      "text to display";
        position    (0.1 0.05);
        size        18;
        bold        yes;
        visible     yes;

        // Optionally override default colour
        // colour    (0 1 1);

        // Optional entry
        timeStamp   yes; //Append solution time to string
    }
    \endverbatim

SourceFiles
    text.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_runTimePostPro_text_H
#define functionObjects_runTimePostPro_text_H

#include "geometryBase.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class vtkRenderer;

namespace Foam
{
namespace functionObjects
{
namespace runTimePostPro
{

/*---------------------------------------------------------------------------*\
                          Class text Declaration
\*---------------------------------------------------------------------------*/

class text
:
    public geometryBase
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        text(const text&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const text&) = delete;


protected:

    // Protected data

        //- Text
        string string_;

        //- Position
        Tuple2<scalar, scalar> position_;

        //- Size
        scalar size_;

        //- Colour
        autoPtr<Function1<vector>> colour_;

        //- Bold flag
        bool bold_;

        //- Time stamp flag
        bool timeStamp_;


public:

    // Constructors

        //- Construct from dictionary
        text
        (
            const runTimePostProcessing& parent,
            const dictionary& dict,
            const HashPtrTable<Function1<vector>, word>& colours
        );


    //- Destructor
    virtual ~text();


    // Member Functions

        //- Add surface(s) to scene
        virtual void addGeometryToScene
        (
            const scalar position,
            vtkRenderer* renderer
        );

        //- Update actors
        virtual void updateActors(const scalar position);

        //- Clear files used to create the object(s)
        virtual bool clear();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace runTimePostPro
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
