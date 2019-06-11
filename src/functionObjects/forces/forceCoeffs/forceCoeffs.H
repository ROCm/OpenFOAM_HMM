/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::functionObjects::forceCoeffs

Group
    grpForcesFunctionObjects

Description
    Extends the forces functionObject by providing coefficients for:
    - drag, side and lift forces
    - roll, pitch and yaw moments
    The data can optionally be output into bins, defined in a given direction.

    The binned data provides the total and consitituent components per bin:
    - total coefficient
    - pressure coefficient contribution
    - viscous coefficient contribution
    - porous coefficient contribution

    Data is written into multiple files in the
    postProcessing/\<functionObjectName\> directory:
    - coefficient.dat   : integrated coefficients over all geometries
    - CdBin.dat         : drag coefficient bins
    - CsBin.dat         : side coefficient bins
    - ClBin.dat         : lift coefficient bins
    - CmRollBin.dat     : roll moment coefficient bins
    - CmPitchBin.dat    : pitch moment coefficient bins
    - CmYawBin.dat      : yaw moment coefficient bins

Usage
    Example of function object specification:
    \verbatim
    forceCoeffs1
    {
        type        forceCoeffs;
        libs        ("libforces.so");
        ...
        log         yes;
        writeFields yes;
        patches     (walls);

        coordinateSystem
        {
            type        xxx;
            dragDir     (1 0 0);
            liftDir     (0 0 1);
        }

        magUInf     100;
        lRef        3.5;
        Aref        2.2;

        binData
        {
            nBin        20;
            direction   (1 0 0);
            cumulative  yes;
        }
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description                          | Required | Default
        type         | Type name: forceCoeffs               | yes |
        log          | Write force data to standard output  | no  | no
        writeFields  | Write force,moment coefficient fields | no | no
        patches      | Patches included in the forces calculation | yes |
        dragDir      | Drag direction                       | no  | (1 0 0)
        sideDir      | Side force direction                 | no  | (0 1 0)
        liftDir      | Lift direction                       | no  | (0 0 1)
        rollAxis     | Roll axis                            | no  | (1 0 0)
        pitchAxis    | Pitch axis                           | no  | (0 1 0)
        yawAxis      | Yaw axis                             | no  | (0 0 1)
        magUInf      | Free stream velocity magnitude       | yes |
        lRef         | Reference length scale for moment calculations | yes |
        Aref         | Reference area                       | yes |
        porosity     | Include porosity contributions       | no  | false
    \endtable

    Bin data is optional, but if the dictionary is present, the entries must
    be defined according to following:
    \table
        nBin         | Number of data bins                    | yes |
        direction    | Direction along which bins are defined | yes |
        cumulative   | Bin data accumulated with incresing distance | yes |
    \endtable

See also
    Foam::functionObject
    Foam::functionObjects::timeControl
    Foam::functionObjects::forces

SourceFiles
    forceCoeffs.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_forceCoeffs_H
#define functionObjects_forceCoeffs_H

#include "forces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class forceCoeffs Declaration
\*---------------------------------------------------------------------------*/

class forceCoeffs
:
    public forces
{
    // Private data

        // Free-stream conditions

            //- Free-stream velocity magnitude
            scalar magUInf_;


        // Reference scales

            //- Reference length [m]
            scalar lRef_;

            //- Reference area [m^2]
            scalar Aref_;


        // File streams

            //- Integrated coefficients
            autoPtr<OFstream> coeffFilePtr_;

            //- Drag coefficient
            autoPtr<OFstream> CdBinFilePtr_;

            //- Side coefficient
            autoPtr<OFstream> CsBinFilePtr_;

            //- Lift coefficient
            autoPtr<OFstream> ClBinFilePtr_;

            //- Roll moment coefficient
            autoPtr<OFstream> CmRollBinFilePtr_;

            //- Pitch moment coefficient
            autoPtr<OFstream> CmPitchBinFilePtr_;

            //- Yaw moment coefficient
            autoPtr<OFstream> CmYawBinFilePtr_;


    // Private Member Functions

        //- No copy construct
        forceCoeffs(const forceCoeffs&) = delete;

        //- No copy assignment
        void operator=(const forceCoeffs&) = delete;


protected:

    // Protected Member Functions

        //- Create the output files
        void createFiles();

        //- Write header for integrated data
        void writeIntegratedHeader(const word& header, Ostream& os) const;

        //- Write header for binned data
        void writeBinHeader(const word& header, Ostream& os) const;

        //- Write integrated data
        void writeIntegratedData
        (
            const word& title,
            const List<Field<scalar>>& coeff
        ) const;

        //- Write binned data
        void writeBinData(const List<Field<scalar>> coeffs, Ostream& os) const;


public:

    //- Runtime type information
    TypeName("forceCoeffs");


    // Constructors

        //- Construct from Time and dictionary
        forceCoeffs
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );


    //- Destructor
    virtual ~forceCoeffs() = default;


    // Member Functions

        //- Read the forces data
        virtual bool read(const dictionary&);

        //- Execute
        virtual bool execute();

        //- Write the forces
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

