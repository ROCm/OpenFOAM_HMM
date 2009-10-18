/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::vtkPV3blockMesh

Description
    Provides a reader interface for OpenFOAM blockMesh to VTK interaction

SourceFiles
    vtkPV3blockMesh.C
    vtkPV3blockMeshConvert.C
    vtkPV3blockMeshUpdate.C
    vtkPV3blockMeshUtils.C

    // Needed by VTK:
    vtkDataArrayTemplateImplicit.txx

\*---------------------------------------------------------------------------*/

#ifndef vtkPV3blockMesh_H
#define vtkPV3blockMesh_H

// do not include legacy strstream headers
#ifndef  VTK_EXCLUDE_STRSTREAM_HEADERS
# define VTK_EXCLUDE_STRSTREAM_HEADERS
#endif

#include "className.H"
#include "fileName.H"
#include "stringList.H"
#include "wordList.H"

#include "primitivePatch.H"

// * * * * * * * * * * * * * Forward Declarations  * * * * * * * * * * * * * //

class vtkDataArraySelection;
class vtkDataSet;
class vtkPoints;
class vtkPV3blockMeshReader;
class vtkRenderer;
class vtkTextActor;
class vtkMultiBlockDataSet;
class vtkPolyData;
class vtkUnstructuredGrid;
class vtkIndent;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Foam class forward declarations
class argList;
class Time;
class blockMesh;

template<class Type> class List;

/*---------------------------------------------------------------------------*\
                     Class vtkPV3blockMesh Declaration
\*---------------------------------------------------------------------------*/

class vtkPV3blockMesh
{
    // Private classes

        //- Bookkeeping for GUI checklists and the multi-block organization
        class partInfo
        {
            const char *name_;
            int block_;
            int start_;
            int size_;

        public:

            partInfo(const char *name, const int blockNo=0)
            :
                name_(name),
                block_(blockNo),
                start_(-1),
                size_(0)
            {}

            //- Return the block holding these datasets
            int block() const
            {
                return block_;
            }

            //- Assign block number, return previous value
            int block(int blockNo)
            {
                int prev = block_;
                block_ = blockNo;
                return prev;
            }

            const char* name() const
            {
                return name_;
            }

            int start() const
            {
                return start_;
            }

            int end() const
            {
                return start_ + size_;
            }

            int size() const
            {
                return size_;
            }

            bool empty() const
            {
                return !size_;
            }

            void reset()
            {
                start_ = -1;
                size_ = 0;
            }

            //- Assign new start and reset the size
            void operator=(const int i)
            {
                start_ = i;
                size_ = 0;
            }

            //- Increment the size
            void operator+=(const int n)
            {
                size_ += n;
            }
        };


    // Private Data

        //- Access to the controlling vtkPV3blockMeshReader
        vtkPV3blockMeshReader* reader_;

        //- Foam time control
        autoPtr<Time> dbPtr_;

        //- Foam mesh
        blockMesh* meshPtr_;

        //- Selected geometrical parts
        boolList partStatus_;

        //- Selected curved edges
        boolList edgeStatus_;

        //- First instance and size of bleckMesh blocks
        //  used to index into partStatus_
        partInfo partInfoBlocks_;

        //- First instance and size of CurvedEdges
        //  only partially used
        partInfo partInfoEdges_;

        //- First instance and size of block corners
        //  only partially used
        partInfo partInfoCorners_;

        //- List of point numbers for rendering to window
        List<vtkTextActor*> pointNumberTextActorsPtrs_;

    // Private Member Functions

        // Convenience method use to convert the readers from VTK 5
        // multiblock API to the current composite data infrastructure
        static void AddToBlock
        (
            vtkMultiBlockDataSet* output,
            vtkDataSet* dataset,
            const partInfo&,
            const label datasetNo,
            const std::string& datasetName
        );

        // Convenience method use to convert the readers from VTK 5
        // multiblock API to the current composite data infrastructure
        static vtkDataSet* GetDataSetFromBlock
        (
            vtkMultiBlockDataSet* output,
            const partInfo&,
            const label datasetNo
        );

        // Convenience method use to convert the readers from VTK 5
        // multiblock API to the current composite data infrastructure
        static label GetNumberOfDataSets
        (
            vtkMultiBlockDataSet* output,
            const partInfo&
        );

        //- Update boolList from GUI selection
        static void updateBoolListStatus
        (
            boolList&,
            vtkDataArraySelection*
        );

        //- Reset data counters
        void resetCounters();

        // Update information helper functions

            //- Internal block info
            void updateInfoBlocks();

            //- block curved edges info
            void updateInfoEdges();

        // Update helper functions

            //- Foam mesh
            void updateFoamMesh();

        // Mesh conversion functions

            //- mesh blocks
            void convertMeshBlocks(vtkMultiBlockDataSet*, int& blockNo);

            //- mesh curved edges
            void convertMeshEdges(vtkMultiBlockDataSet*, int& blockNo);

            //- mesh corners
            void convertMeshCorners(vtkMultiBlockDataSet*, int& blockNo);


       // GUI selection helper functions

            //- Retrieve the current selections
            static wordHashSet getSelected(vtkDataArraySelection*);

            //- Retrieve a sub-list of the current selections
            static wordHashSet getSelected
            (
                vtkDataArraySelection*,
                const partInfo&
            );

            //- Retrieve the current selections
            static stringList getSelectedArrayEntries(vtkDataArraySelection*);

            //- Retrieve a sub-list of the current selections
            static stringList getSelectedArrayEntries
            (
                vtkDataArraySelection*,
                const partInfo&
            );

            //- Set selection(s)
            static void setSelectedArrayEntries
            (
                vtkDataArraySelection*,
                const stringList&
            );


        //- Disallow default bitwise copy construct
        vtkPV3blockMesh(const vtkPV3blockMesh&);

        //- Disallow default bitwise assignment
        void operator=(const vtkPV3blockMesh&);


public:

    //- Static data members

        ClassName("vtkPV3blockMesh");


    // Constructors

        //- Construct from components
        vtkPV3blockMesh
        (
            const char* const FileName,
            vtkPV3blockMeshReader* reader
        );


    //- Destructor

        ~vtkPV3blockMesh();


    // Member Functions

        //- Update
        void updateInfo();

        void Update(vtkMultiBlockDataSet* output);

        //- Clean any storage
        void CleanUp();

        //- Add/remove point numbers to/from the view
        void renderPointNumbers(vtkRenderer*, const bool show);

     // Access

        //- Debug information
        void PrintSelf(ostream&, vtkIndent) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
