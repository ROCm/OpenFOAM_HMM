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
    Foam::vtkMesh

Description
    Encapsulation of VTK mesh data. Holds mesh or meshsubset and
    polyhedral-cell decomposition on it.

SourceFiles
    vtkMesh.C

\*---------------------------------------------------------------------------*/

#ifndef vtkMesh_H
#define vtkMesh_H

#include "meshSubsetHelper.H"
#include "vtkTopo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class vtkMesh Declaration
\*---------------------------------------------------------------------------*/

class vtkMesh
:
    public meshSubsetHelper
{
    // Private data

        //- Current decomposition of topology
        mutable autoPtr<vtkTopo> topoPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        vtkMesh(const vtkMesh&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const vtkMesh&) = delete;


public:

    // Constructors

        //- Construct from components
        vtkMesh(fvMesh& baseMesh, const word& setName = word::null)
        :
            meshSubsetHelper(baseMesh, meshSubsetHelper::SET, setName)
        {}


    // Member Functions

        // Access

            //- Topology
            const vtkTopo& topo() const
            {
                if (topoPtr_.empty())
                {
                    topoPtr_.reset(new vtkTopo(mesh()));
                }
                return topoPtr_();
            }


            //- Number of field cells
            label nFieldCells() const
            {
                return topo().cellTypes().size();
            }


            //- Number of field points
            label nFieldPoints() const
            {
                return mesh().nPoints() + topo().addPointCellLabels().size();
            }


        // Edit

            //- Read mesh, forcing topo update if necessary
            polyMesh::readUpdateState readUpdate()
            {
                polyMesh::readUpdateState meshState
                    = meshSubsetHelper::readUpdate();

                if (meshState != polyMesh::UNCHANGED)
                {
                    topoPtr_.clear();
                }

                return meshState;
            }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
