/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::vtk::lagrangianWriter

Description
    Write fields (internal).

SourceFiles
    lagrangianWriter.C
    lagrangianWriterTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef foamVtkLagrangianWriter_H
#define foamVtkLagrangianWriter_H

#include "OFstream.H"
#include "Cloud.H"
#include "volFields.H"
#include "pointFields.H"
#include "foamVtkOutputOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
class volPointInterpolation;

namespace vtk
{

/*---------------------------------------------------------------------------*\
                      Class lagrangianWriter Declaration
\*---------------------------------------------------------------------------*/

class lagrangianWriter
{
    // Private Member Data

        //- Reference to the OpenFOAM mesh (or subset)
        const fvMesh& mesh_;

        //- Commonly used query
        const bool legacy_;

        //- Write lagrangian as cell data (verts) or point data?
        const bool useVerts_;

        autoPtr<vtk::formatter> format_;

        const word cloudName_;

        std::ofstream os_;

        label nParcels_;


    // Private Member Functions

        //- Begin piece
        void beginPiece();

        //- Write positions
        void writePoints();

        //- Write vertex (cells)
        void writeVertsLegacy();

        //- Write vertex (cells)
        void writeVerts();


        //- Disallow default bitwise copy construct
        lagrangianWriter(const lagrangianWriter&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const lagrangianWriter&) = delete;


public:

    // Constructors

        //- Construct from components
        lagrangianWriter
        (
            const fvMesh& mesh,
            const word& cloudName,
            const fileName& baseName,
            const vtk::outputOptions outOpts,
            const bool dummyCloud = false
        );


    //- Destructor
    ~lagrangianWriter();


    // Member Functions

        inline std::ofstream& os()
        {
            return os_;
        }

        inline vtk::formatter& format()
        {
            return format_();
        }

        inline label nParcels() const
        {
            return nParcels_;
        }

        void beginParcelData(label nFields);
        void endParcelData();

        //- Write file footer
        void writeFooter();

        //- Write IOField
        template<class Type>
        void writeIOField(const wordList& objectNames);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "foamVtkLagrangianWriterTemplates.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
