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
    Foam::patchProbes

Description
    Set of locations to sample at patches

    Call write() to sample and write files.
    - find nearest location on nearest face
    - update *this with location (so header contains 'snapped' locations
    - use *this as the sampling location

    Example of function object specification:
    \verbatim
    patchProbes
    {
        type            patchProbes;
        libs            ( "libsampling.so" );

        // Name of the directory for probe data
        name            patchProbes;

        // Patches to sample (wildcards allowed)
        patches         (".*inl.*");

        // Write at same frequency as fields
        writeControl    writeTime;
        writeInterval   1;

        // Fields to be probed
        fields          (p U);

        // Locations to probe. These get snapped onto the nearest point
        // on the selected patches
        probeLocations
        (
            ( -100 0 0.01 )      // at inlet
        );
    }
    \endverbatim


SourceFiles
    patchProbes.C

\*---------------------------------------------------------------------------*/

#ifndef patchProbes_H
#define patchProbes_H

#include "probes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class fvMesh;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                          Class patchProbes Declaration
\*---------------------------------------------------------------------------*/

class patchProbes
:
    public probes
{
protected:

    // Protected data

        //- Patches to sample
        wordReList patchNames_;


    // Protected Member Functions

        //- Sample and write a particular volume field
        template<class Type>
        void sampleAndWrite
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        //- Sample and write a particular surface field
        template<class Type>
        void sampleAndWrite
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>&
        );

        //- Sample and write all the fields of the given type
        template<class Type>
        void sampleAndWrite(const fieldGroup<Type>&);

         //- Sample and write all the surface fields of the given type
        template<class Type>
        void sampleAndWriteSurfaceFields(const fieldGroup<Type>&);

        //- Sample a volume field at all locations
        template<class Type>
        tmp<Field<Type>> sample
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        //- Sample a surface field at all locations
        template<class Type>
        tmp<Field<Type>> sample
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>&
        ) const;

        //- Sample a single field on all sample locations
        template<class Type>
        tmp<Field<Type>> sample(const word& fieldName) const;

        //- Find elements containing patchProbes
        virtual void findElements(const fvMesh&);


private:

        //- Disallow default bitwise copy construct
        patchProbes(const patchProbes&);

        //- Disallow default bitwise assignment
        void operator=(const patchProbes&);


public:

    //- Runtime type information
    TypeName("patchProbes");


    // Constructors

        //- Construct from Time and dictionary
        patchProbes
        (
            const word& name,
            const Time& time,
            const dictionary& dict,
            const bool loadFromFiles = false,
            const bool readFields = true
        );

    //- Destructor
    virtual ~patchProbes();


    //- Public members

        //- Sample and write
        virtual bool write();

        //- Read
        virtual bool read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "patchProbesTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
