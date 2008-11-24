/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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
    Foam::minMaxFields

Description
    Calculates scalar minimim and maximum field values.

    For variables with rank > 0, computes the magnitude of the min/max
    values.

    Data written to the file \<timeDir\>/minMaxFields.dat


SourceFiles
    minMaxFields.C
    IOminMaxFields.H

\*---------------------------------------------------------------------------*/

#ifndef minMaxFields_H
#define minMaxFields_H

#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "HashSet.H"
#include "OFstream.H"
#include "Switch.H"
#include "pointFieldFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                       Class minMaxFields Declaration
\*---------------------------------------------------------------------------*/

class minMaxFields
{
protected:

    // Private data

        //- Name of this set of forces,
        //  Also used as the name of the probes directory.
        word name_;

        const objectRegistry& obr_;

        //- on/off switch
        bool active_;

        //- Switch to send output to Info as well as to file
        Switch log_;

        //- Patches to integrate forces over
        wordList fieldSet_;


        //- Forces/moment file ptr
        autoPtr<OFstream> minMaxFieldsFilePtr_;


    // Private Member Functions

        //- If the forces file has not been created create it
        void makeFile();

        //- Disallow default bitwise copy construct
        minMaxFields(const minMaxFields&);

        //- Disallow default bitwise assignment
        void operator=(const minMaxFields&);

        //- Output file header information
        virtual void writeFileHeader();


public:

    //- Runtime type information
    TypeName("minMaxFields");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        minMaxFields
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    // Destructor

        virtual ~minMaxFields();


    // Member Functions

        //- Return name of the set of forces
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the forces data
        virtual void read(const dictionary&);

        //- Execute
        virtual void execute();

        //- Calculate the field min/max
        template<class Type>
        void calcMinMaxFields(const word& fieldName);

        //- Write the minMaxFields
        virtual void write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const pointField&)
        {}
};


// Template specialisation for scalar fields
template<>
void minMaxFields::calcMinMaxFields<scalar>(const word& fieldName);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "minMaxFieldsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
