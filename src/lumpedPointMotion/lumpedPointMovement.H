/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenCFD Ltd.
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
    Foam::lumpedPointMovement

Description
    The movement 'driver' that describes initial point locations, the
    segmentation for pressure integration, the current state of the
    points/rotations, and forwarding to the communication
    coordinator. The 'lumpedPointIOMovement' class is simply a
    registered version of the same.

SourceFiles
    lumpedPointMovement.C

\*---------------------------------------------------------------------------*/

#ifndef lumpedPointMovement_H
#define lumpedPointMovement_H

#include "dictionary.H"
#include "scalarList.H"
#include "scalarField.H"
#include "pointField.H"
#include "vectorField.H"
#include "tensorField.H"
#include "vector.H"
#include "interpolationWeights.H"
#include "IOobject.H"
#include "tmp.H"
#include "faceZoneMeshFwd.H"
#include "externalFileCoupler.H"
#include "lumpedPointState.H"
#include "boundBox.H"
#include "Enum.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// Forward declarations
class polyMesh;

/*---------------------------------------------------------------------------*\
                     Class lumpedPointMovement Declaration
\*---------------------------------------------------------------------------*/

class lumpedPointMovement
{
public:

    //- Output format types
    enum class outputFormatType
    {
        PLAIN,
        DICTIONARY
    };

    // Static data

    //- Names for the output format types
    static const Enum<outputFormatType> formatNames;


private:

    // Private data

        //- Reference axis for the locations
        vector axis_;

        //- The locations of the lumped points within the reference axis
        //  The interpolator needs scalarField, not scalarList.
        scalarField locations_;

        //- The division when calculating pressure forces (0-1)
        scalar division_;

        //- The relaxation factor when moving points (default: 1)
        scalar relax_;

        //- The interpolation type (linear|spline)
        word interpolationScheme_;

        //- Optional owner information
        label ownerId_;

        //- The bounding box (if set)
        boundBox boundBox_;

        //- The offset of patch points to compared to the locations
        point centre_;

        //- Calculate centre based on the bounding box
        bool autoCentre_;

        //- Dictionary of controls for force calculation
        dictionary forcesDict_;

        //- Communication control
        externalFileCoupler coupler_;

        //- File io
        word inputName_;
        word outputName_;
        lumpedPointState::inputFormatType inputFormat_;
        outputFormatType outputFormat_;


    // Demand-driven private data

        //- The initial state (positions, rotations)
        lumpedPointState state0_;

        //- The current state (positions, rotations)
        //  Eg, as response from external application
        lumpedPointState state_;

        //- Threshold locations for pressure forces
        mutable scalarField* thresholdPtr_;

        //- User-specified interpolator
        mutable autoPtr<interpolationWeights> interpolatorPtr_;

        //- Pressure zones (only used from the master patch)
        mutable List<labelList> faceZones_;


    // Private Member Functions

        //- Calculate threshold locations
        void calcThresholds() const;

        //- Classify the position to be located in one of the threshold zones
        label threshold(scalar pos) const;


        //- Disallow default bitwise copy constructor
        lumpedPointMovement(const lumpedPointMovement&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const lumpedPointMovement&) = delete;

public:

    // Static data members

        //- The standard dictionary name
        static const word dictionaryName;

        //- Helper: return IOobject with optionally absolute path provided
        static IOobject selectIO(const IOobject&, const fileName&);


    // Constructors

        //- Construct null
        lumpedPointMovement();

        //- Construct from dictionary, optionally with some owner information
        lumpedPointMovement(const dictionary& dict, label ownerId=-1);


    //- Destructor
    virtual ~lumpedPointMovement();


    // Member Functions

        //- Update settings from dictionary
        void readDict(const dictionary& dict);


        //- If no number of lumped points (locations) were specified
        inline bool empty() const;

        //- The number of lumped points (number of locations)
        inline label size() const;

        //- The normalized reference axis
        inline const vector& axis() const;

        //- Read access to the locations
        inline const scalarField& locations() const;

        //- The division (0-1) when calculating pressure forces
        inline scalar division() const;

        //- An owner Id, if needed for bookkeeping purposes
        inline label ownerId() const;

        //- Change the owner id, if needed for bookkeeping purposes
        inline void ownerId(label id);


        //- Threshold locations for pressure forces
        inline const scalarField& thresholds() const;

        //- Classify the position to be located in one of the threshold zones
        inline label threshold(const point& position) const;


        //- Communication control
        inline const externalFileCoupler& coupler() const;

        //- Communication control
        inline externalFileCoupler& coupler();

        //- The initial state (positions/rotations)
        inline const lumpedPointState& state0() const;

        //- The current state (positions/rotations)
        inline const lumpedPointState& state() const;

        //- The relaxation factor when changing states
        inline scalar relax() const;

        //- The relaxation factor when changing states
        inline scalar& relax();

        //- The input (state) file name
        inline const word& inputName() const;

        //- The output (forces) file name
        inline const word& outputName() const;

        //- The input (state) file format
        inline lumpedPointState::inputFormatType inputFormat() const;

        //- The output (forces) file format
        inline lumpedPointMovement::outputFormatType outputFormat() const;


        //- Define the bounding-box required to enclose the specified patches.
        //  Calculates the centre as required.
        //
        //  \param mesh     The volume mesh reference
        //  \param patchIds The patch ids to be included in the bounding box.
        //  \param points0  The initial mesh points, prior to movement
        void setBoundBox
        (
            const polyMesh& mesh,
            const labelUList& patchIds,
            const pointField& points0
        );


        //- Define the pressure-zones mapping for faces in the specified
        //  patches.
        //  The face centres are compared to the threshold positions,
        //  which are determined by locations along the defined axis.
        //
        //  \param mesh     The volume mesh reference
        //  \param patchIds The patch ids to be included in the mapping
        //  \param points0  The initial mesh points, prior to movement
        void setMapping
        (
            const polyMesh& mesh,
            const labelUList& patchIds,
            const pointField& points0
        );


        //- True if the pressure-zones mapping has already been performed
        inline bool hasMapping() const;

        //- Return the pressure-zones mapping with the associated
        //  patch face ids.
        inline const List<labelList>& zones() const;


        //- The forces and moments acting on each pressure-zone.
        //  The zones must be previously defined via setMapping.
        bool forcesAndMoments
        (
            const polyMesh& pmesh,
            List<vector>& forces,
            List<vector>& moments
        ) const;


        //- Displace points according to the current state
        tmp<pointField> displacePoints
        (
            const pointField& points0,
            const labelList& pointLabels
        ) const;


        //- Displace points according to specified state
        tmp<pointField> displacePoints
        (
            const lumpedPointState& state,
            const pointField& points0,
            const labelList& pointLabels
        ) const;


        //- Interpolation weights
        const interpolationWeights& interpolator() const;

        //- Write axis, locations, division as a dictionary
        void writeDict(Ostream& os) const;


        //- Write points, forces
        bool writeData
        (
            const UList<vector>& forces
        ) const;

        //- Write points, forces, moments
        bool writeData
        (
            const UList<vector>& forces,
            const UList<vector>& moments
        ) const;


        //- Read state from file, applying relaxation as requested
        bool readState();

        //- Write state as VTK PolyData format.
        void writeStateVTP(const fileName& file) const;

        //- Write forces on points as VTK PolyData format.
        void writeForcesAndMomentsVTP
        (
            const fileName& file,
            const UList<vector>& forces,
            const UList<vector>& moments
        ) const;


        //- Write pressure-zones geometry, write as VTK PolyData format.
        void writeZonesVTP
        (
            const fileName& file,
            const polyMesh& mesh,
            const pointField& points0
        ) const;


        //- Write displaced geometry according to the current state,
        //  write as VTK PolyData format.
        void writeVTP
        (
            const fileName& file,
            const polyMesh& mesh,
            const labelUList& patchIds,
            const pointField& points0
        ) const;

        //- Write displaced geometry according to the specified state,
        //  write as VTK PolyData format.
        void writeVTP
        (
            const fileName& file,
            const lumpedPointState& state,
            const polyMesh& mesh,
            const labelUList& patchLst,
            const pointField& points0
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "lumpedPointMovementI.H"

#endif

// ************************************************************************* //
