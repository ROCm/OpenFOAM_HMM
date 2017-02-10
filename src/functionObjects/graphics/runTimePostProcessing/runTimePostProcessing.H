/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd.
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
    Foam::functionObjects::runTimePostPro::runTimePostProcessing

Group
    grpGraphicsFunctionObjects

Description
    Generate images during run-time.

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

#ifndef functionObjects_runTimePostProcessing_H
#define functionObjects_runTimePostProcessing_H

#include "fvMeshFunctionObject.H"
#include "objectRegistry.H"
#include "mapPolyMesh.H"
#include "PtrList.H"
#include "scene.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class vtkRenderer;
class vtkRenderWindow;

namespace Foam
{
namespace functionObjects
{
namespace runTimePostPro
{
    class pointData;
    class pathline;
    class surface;
    class text;
}

/*---------------------------------------------------------------------------*\
                   Class runTimePostProcessing Declaration
\*---------------------------------------------------------------------------*/

class runTimePostProcessing
:
    public fvMeshFunctionObject
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
        runTimePostPro::scene scene_;

        //- List of points
        PtrList<runTimePostPro::pointData> points_;

        //- List of lines
        PtrList<runTimePostPro::pathline> lines_;

        //- List of surfaces
        PtrList<runTimePostPro::surface> surfaces_;

        //- List of text
        PtrList<runTimePostPro::text> text_;


        // Private Member Functions

            //- Helper function to read scene objects
            template<class Type>
            void readObjects
            (
                const dictionary& dict,
                PtrList<Type>& objects
            ) const;


public:

    //- Runtime type information
    TypeName("runTimePostProcessing");


    // Constructors

        //- Construct from dictionary
        runTimePostProcessing
        (
            const word& name,
            const Time& runTime,
            const dictionary&dict
        );


    //- Desructor
    virtual ~runTimePostProcessing();


    // Member Functions

        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Read the field min/max data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "runTimePostProcessingTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
