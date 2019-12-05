/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
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
    Foam::faGlobalMeshData

Description
    Various mesh related information for a parallel run

Author
    Zeljko Tukovic, FMENA
    Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#ifndef faGlobalMeshData_H
#define faGlobalMeshData_H

#include "faProcessorTopology.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class faMesh;

/*---------------------------------------------------------------------------*\
                        Class faGlobalMeshData Declaration
\*---------------------------------------------------------------------------*/

class faGlobalMeshData
:
    public faProcessorTopology
{
    // Private data

        //- Reference to mesh
        const faMesh& mesh_;

        // Globally shared point addressing

            //- Total number of global points
            label nGlobalPoints_;

            //- Indices of local points that are globally shared
            labelList sharedPointLabels_;

            //- Indices of globally shared points in the master list
            //  This list contains all the shared points in the mesh
            labelList sharedPointAddr_;

    // Private Member Functions

        //- No copy construct
        faGlobalMeshData(const faGlobalMeshData&) = delete;

        //- No copy assignment
        void operator=(const faGlobalMeshData&) = delete;


public:

    //- Runtime type information
    ClassName("faGlobalMeshData");

    //- Construct from mesh
    faGlobalMeshData(const faMesh& mesh);

    //- Destructor
    ~faGlobalMeshData();


    // Member Functions

        // Access

            //- Return mesh reference
            const faMesh& mesh() const;

        // Globally shared point addressing

            //- Return number of globally shared points
            label nGlobalPoints() const
            {
                return nGlobalPoints_;
            }

            //- Return indices of local points that are globally shared
            const labelList& sharedPointLabels() const
            {
                return sharedPointLabels_;
            }

            //- Return addressing into the complete globally shared points
            //  list
            const labelList& sharedPointAddr() const
            {
                return sharedPointAddr_;
            }

            //- Change global mesh data given a topological change.
            void updateMesh();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
