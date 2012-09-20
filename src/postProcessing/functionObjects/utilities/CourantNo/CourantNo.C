/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "CourantNo.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "dictionary.H"
#include "zeroGradientFvPatchFields.H"
#include "fvcSurfaceIntegrate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::CourantNo, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::CourantNo::rho
(
    const surfaceScalarField& phi
) const
{
    if (phi.dimensions() == dimMass/dimTime)
    {
        return (obr_.lookupObject<volScalarField>(rhoName_));
    }
    else
    {
        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "rho",
                    mesh.time().timeName(),
                    mesh
                ),
                mesh,
                dimensionedScalar("rho", dimless, 1.0)
            )
        );
    }
}


void Foam::CourantNo::makeFile()
{
    // Create the CourantNo file if not already created
    if (CourantNoFilePtr_.empty())
    {
        if (debug)
        {
            Info<< "Creating CourantNo file." << endl;
        }

        // File update
        if (Pstream::master())
        {
            fileName CourantNoDir;
            if (Pstream::parRun())
            {
                // Put in undecomposed case (Note: gives problems for
                // distributed data running)
                CourantNoDir =
                    obr_.time().path()/".."/name_/obr_.time().timeName();
            }
            else
            {
                CourantNoDir =
                    obr_.time().path()/name_/obr_.time().timeName();
            }

            // Create directory if does not exist.
            mkDir(CourantNoDir);

            // Open new file at start up
            CourantNoFilePtr_.reset
            (
                new OFstream(CourantNoDir/(type() + ".dat"))
            );

            // Add headers to output data
            writeFileHeader();
        }
    }
}


void Foam::CourantNo::writeFileHeader()
{
    if (CourantNoFilePtr_.valid())
    {
        CourantNoFilePtr_()
            << "# Time" << token::TAB << "min" << token::TAB << "position(min)";

        if (Pstream::parRun())
        {
            CourantNoFilePtr_() << token::TAB << "proc";
        }

        CourantNoFilePtr_()
            << token::TAB << "max" << token::TAB << "position(max)";

        if (Pstream::parRun())
        {
            CourantNoFilePtr_() << token::TAB << "proc";
        }

        CourantNoFilePtr_() << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CourantNo::CourantNo
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    active_(true),
    log_(false),
    phiName_("phi"),
    rhoName_("rho"),
    CourantNoFilePtr_(NULL)
{
    // Check if the available mesh is an fvMesh, otherwise deactivate
    if (!isA<fvMesh>(obr_))
    {
        active_ = false;
        WarningIn
        (
            "CourantNo::CourantNo"
            "("
                "const word&, "
                "const objectRegistry&, "
                "const dictionary&, "
                "const bool"
            ")"
        )   << "No fvMesh available, deactivating." << nl
            << endl;
    }

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::CourantNo::~CourantNo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::CourantNo::read(const dictionary& dict)
{
    if (active_)
    {
        phiName_ = dict.lookupOrDefault<word>("phiName", "phi");
        rhoName_ = dict.lookupOrDefault<word>("rhoName", "rho");
        log_ = dict.lookupOrDefault<Switch>("log", false);
    }
}


void Foam::CourantNo::execute()
{
    // Do nothing - only valid on write
}


void Foam::CourantNo::end()
{
    // Do nothing - only valid on write
}


void Foam::CourantNo::write()
{
    if (active_)
    {
        makeFile();

        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        const surfaceScalarField& phi =
            mesh.lookupObject<surfaceScalarField>(phiName_);

        volScalarField CourantNo
        (
            IOobject
            (
                "CourantNo",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ
            ),
            mesh,
            dimensionedScalar("0", dimless, 0.0),
            zeroGradientFvPatchScalarField::typeName
        );

        scalarField& iField = CourantNo.internalField();

        const scalarField sumPhi
        (
            fvc::surfaceSum(mag(phi))().internalField()
           /rho(phi)().internalField()
        );

        iField = 0.5*sumPhi/mesh.V().field()*mesh.time().deltaTValue();

        CourantNo.correctBoundaryConditions();

        const label procI = Pstream::myProcNo();

        labelList minIs(Pstream::nProcs());
        scalarList minVs(Pstream::nProcs());
        List<vector> minCs(Pstream::nProcs());
        minIs[procI] = findMin(iField);
        minVs[procI] = iField[minIs[procI]];
        minCs[procI] = mesh.C()[minIs[procI]];

        Pstream::gatherList(minIs);
        Pstream::gatherList(minVs);
        Pstream::gatherList(minCs);

        labelList maxIs(Pstream::nProcs());
        scalarList maxVs(Pstream::nProcs());
        List<vector> maxCs(Pstream::nProcs());
        maxIs[procI] = findMax(iField);
        maxVs[procI] = iField[maxIs[procI]];
        maxCs[procI] = mesh.C()[maxIs[procI]];

        Pstream::gatherList(maxIs);
        Pstream::gatherList(maxVs);
        Pstream::gatherList(maxCs);

        if (Pstream::master())
        {
            label minI = findMin(minVs);
            scalar minValue = minVs[minI];
            const vector& minC = minCs[minI];

            label maxI = findMax(maxVs);
            scalar maxValue = maxVs[maxI];
            const vector& maxC = maxCs[maxI];

            CourantNoFilePtr_()
                << obr_.time().value() << token::TAB
                << minValue << token::TAB << minC;

            if (Pstream::parRun())
            {
                CourantNoFilePtr_() << token::TAB << minI;
            }

            CourantNoFilePtr_()
                << token::TAB << maxValue << token::TAB << maxC;

            if (Pstream::parRun())
            {
                CourantNoFilePtr_() << token::TAB << maxI;
            }

            CourantNoFilePtr_() << endl;

            if (log_)
            {
                Info<< type() << " output:" << nl
                    << "    min = " << minValue << " at position " << minC;

                if (Pstream::parRun())
                {
                    Info<< " on processor " << minI;
                }

                Info<< nl << "    max = " << maxValue << " at position "
                    << maxC;

                if (Pstream::parRun())
                {
                    Info<< " on processor " << maxI;
                }

                Info<< nl << endl;
            }
        }


        CourantNo.write();
    }
}


// ************************************************************************* //
