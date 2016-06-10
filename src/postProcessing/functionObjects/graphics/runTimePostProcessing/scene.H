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
    Foam::scene

Description

SourceFiles
    scene.C

\*---------------------------------------------------------------------------*/

#ifndef scene_H
#define scene_H

// OpenFOAM includes
#include "dictionary.H"
#include "DataEntry.H"
#include "vector.H"
#include "point.H"
#include "boundBox.H"
#include "NamedEnum.H"
#include "HashPtrTable.H"
#include "vector.H"

// VTK includes
#include "vtkSmartPointer.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class vtkRenderer;
class vtkRenderWindow;

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class scene Declaration
\*---------------------------------------------------------------------------*/

class scene
{
public:

    enum modeType{mtStatic, mtFlightPath};

    NamedEnum<modeType, 2> modeTypeNames_;


private:

    // Private data

        //- Reference to the object registry
        const objectRegistry& obr_;

        //- Object name
        const word name_;


    // Private Member Functions

        //- Read camera properties
        void readCamera(const dictionary& dict);

        //- Read solour properties
        void readColours(const dictionary& dict);

        //- Disallow default bitwise copy construct
        scene(const scene&);

        //- Disallow default bitwise assignment
        void operator=(const scene&);


protected:

    // Protected data

        //- Colours
        HashPtrTable<DataEntry<vector>, word> colours_;


        // Camera settings

            //- Mode
            modeType mode_;

            //- Position
            autoPtr<DataEntry<point> > cameraPosition_;

            //- Focal point
            autoPtr<DataEntry<point> > cameraFocalPoint_;

            //- Up direction
            autoPtr<DataEntry<vector> > cameraUp_;

            //- View angle
            autoPtr<DataEntry<scalar> > cameraViewAngle_;


        // Scene management

            //- Clipping box
            boundBox clipBox_;

            //- Parallel projection flag
            bool parallelProjection_;

            //- Number of frames
            label nFrameTotal_;

            //- Start position [0-1]
            scalar startPosition_;

            //- Position [0-1]
            scalar position_;

            //- Change in position per frame
            scalar dPosition_;

            //- Index of current frame
            label currentFrameI_;

            //- Name prefix of output
            word outputName_;


    // Protected Member Functions

        void setCamera(vtkRenderer* renderer) const;

        string frameIndexStr() const;


public:

    // Constructors

        //- Construct from components
        scene(const objectRegistry& obr, const word& name);


    //- Destructor
    virtual ~scene();


    // Member Functions

        // Access

            //- Return the colours
            const HashPtrTable<DataEntry<vector>, word>& colours() const;

            //- Return the current frame index
            label frameIndex() const;

            //- Return the current position (in range 0-1)
            scalar position() const;

            void read(const dictionary& dict);

            void initialise(vtkRenderer* renderer, const word& outputName);

            //- Main control loop
            bool loop(vtkRenderer* renderer);

            //- Save image to file
            void saveImage(vtkRenderWindow* renderWindow) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
