/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2016 OpenCFD Ltd.
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
    Foam::CloudToVTK

Group
    grpLagrangianIntermediateFunctionObjects

Description
    Write cloud data in VTK format

SourceFiles
    CloudToVTK.C

\*---------------------------------------------------------------------------*/

#ifndef CloudToVTK_H
#define CloudToVTK_H

#include "CloudFunctionObject.H"
#include "volFields.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class CloudToVTK Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class CloudToVTK
:
    public CloudFunctionObject<CloudType>
{

protected:

    // Protected data

        //- Ascii/binary output flag
        bool binary_;


    // Protected Member Functions

        //- Write post-processing info
        virtual void write();

        //- Helper function to write VTK data
        void writeData
        (
            std::ostream& vtkOs,
            const bool binary,
            const List<floatScalar>& data
        ) const;

        //- Helper function to write VTK field data
        template<class Type>
        void writeFieldData
        (
            std::ostream& vtkOs,
            const bool binary,
            const List<floatScalar>& data,
            const word& title,
            const label nParcels
        ) const;


public:

    //- Runtime type information
    TypeName("cloudToVTK");


    // Constructors

        //- Construct from dictionary
        CloudToVTK
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        CloudToVTK(const CloudToVTK<CloudType>& c);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObject<CloudType>> clone() const
        {
            return autoPtr<CloudFunctionObject<CloudType>>
            (
                new CloudToVTK<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~CloudToVTK();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "CloudToVTK.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
