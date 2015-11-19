/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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
    Foam::runTimePostProcessing

Description
    Function object to generate images during run-time.

    The functionality makes use of the VTK libraries (see http://www.vtk.org)
    which provide a broad set of functionality for scene composition and
    manipulation.

    Images are generated using a combination of function object output, and
    additional data e.gg triangulated surfaces and text.  Current capabilities
    include support for:
    - Camera
    - Objects
    - Points
    - Lines
    - Surfaces
    - Scalar bars
    - Annotations
    - Selection of colour maps
    Scene configuration is performed using standard OpenFOAM dictionaries, using
    the main headings of: output=, camera, colours, points, lines,
    surfaces and text.

SourceFiles
    runTimePostProcessing.C

\*---------------------------------------------------------------------------*/

#ifndef runTimePostProcessing_H
#define runTimePostProcessing_H

#include "functionObjectState.H"
#include "objectRegistry.H"
#include "mapPolyMesh.H"
#include "PtrList.H"
#include "scene.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class vtkRenderer;
class vtkRenderWindow;

namespace Foam
{

class dictionary;
class pointData;
class pathline;
class surface;
class text;

/*---------------------------------------------------------------------------*\
                   Class runTimePostProcessing Declaration
\*---------------------------------------------------------------------------*/

class runTimePostProcessing
:
    public functionObjectState
{
private:

    // Private data

        // Output
        struct outputType
        {
            word name_;
            label width_;
            label height_;
        };

        //- Output instance
        outputType output_;

        //- Scene manager
        scene scene_;

        //- List of points
        PtrList<pointData> points_;

        //- List of lines
        PtrList<pathline> lines_;

        //- List of surfaces
        PtrList<surface> surfaces_;

        //- List of text
        PtrList<text> text_;


        // Private Member Functions

            //- Helper function to read scene objects
            template<class Type>
            void readObjects
            (
                const dictionary& dict,
                PtrList<Type>& objects
            ) const;


protected:

    // Protected data

        //- Reference to the database
        const objectRegistry& obr_;

        //- on/off switch
        bool active_;


public:

    //- Runtime type information
    TypeName("runTimePostProcessing");


    // Constructors

        //- Construct from dictionary
        runTimePostProcessing
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );

    //- Desructor
    virtual ~runTimePostProcessing();


    // Member Functions

        virtual const objectRegistry& obr() const
        {
            return obr_;
        }

        //- Read the field min/max data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Write
        virtual void write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const polyMesh&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "runTimePostProcessingTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
