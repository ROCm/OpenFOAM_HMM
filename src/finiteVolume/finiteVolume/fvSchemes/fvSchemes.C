/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "fvSchemes.H"
#include "Time.H"
#include "steadyStateDdtScheme.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::fvSchemes::debug(Foam::debug::debugSwitch("fvSchemes", 0));


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::fvSchemes::clear()
{
    ddtSchemes_.clear();
    ddtSchemeDefault_.clear();

    d2dt2Schemes_.clear();
    d2dt2SchemeDefault_.clear();

    interpolationSchemes_.clear();
    interpolationSchemeDefault_.clear();

    divSchemes_.clear();
    divSchemeDefault_.clear();

    gradSchemes_.clear();
    gradSchemeDefault_.clear();

    snGradSchemes_.clear();
    snGradSchemeDefault_.clear();

    laplacianSchemes_.clear();
    laplacianSchemeDefault_.clear();

    // Do not clear fluxRequired settings
}


void Foam::fvSchemes::read(const dictionary& dict)
{
    if (dict.found("ddtSchemes"))
    {
        ddtSchemes_ = dict.subDict("ddtSchemes");
    }
    else
    {
        ddtSchemes_.set("default", "none");
    }

    if
    (
        ddtSchemes_.found("default")
     && word(ddtSchemes_.lookup("default")) != "none"
    )
    {
        ddtSchemeDefault_ = ddtSchemes_.lookup("default");
        steady_ =
        (
            word(ddtSchemeDefault_)
         == fv::steadyStateDdtScheme<scalar>::typeName
        );
    }


    if (dict.found("d2dt2Schemes"))
    {
        d2dt2Schemes_ = dict.subDict("d2dt2Schemes");
    }
    else
    {
        d2dt2Schemes_.set("default", "none");
    }

    if
    (
        d2dt2Schemes_.found("default")
     && word(d2dt2Schemes_.lookup("default")) != "none"
    )
    {
        d2dt2SchemeDefault_ = d2dt2Schemes_.lookup("default");
    }


    if (dict.found("interpolationSchemes"))
    {
        interpolationSchemes_ = dict.subDict("interpolationSchemes");
    }
    else if (!interpolationSchemes_.found("default"))
    {
        interpolationSchemes_.add("default", "linear");
    }

    if
    (
        interpolationSchemes_.found("default")
     && word(interpolationSchemes_.lookup("default")) != "none"
    )
    {
        interpolationSchemeDefault_ =
            interpolationSchemes_.lookup("default");
    }


    divSchemes_ = dict.subDict("divSchemes");

    if
    (
        divSchemes_.found("default")
     && word(divSchemes_.lookup("default")) != "none"
    )
    {
        divSchemeDefault_ = divSchemes_.lookup("default");
    }


    gradSchemes_ = dict.subDict("gradSchemes");

    if
    (
        gradSchemes_.found("default")
     && word(gradSchemes_.lookup("default")) != "none"
    )
    {
        gradSchemeDefault_ = gradSchemes_.lookup("default");
    }


    if (dict.found("snGradSchemes"))
    {
        snGradSchemes_ = dict.subDict("snGradSchemes");
    }
    else if (!snGradSchemes_.found("default"))
    {
        snGradSchemes_.add("default", "corrected");
    }

    if
    (
        snGradSchemes_.found("default")
     && word(snGradSchemes_.lookup("default")) != "none"
    )
    {
        snGradSchemeDefault_ = snGradSchemes_.lookup("default");
    }


    laplacianSchemes_ = dict.subDict("laplacianSchemes");

    if
    (
        laplacianSchemes_.found("default")
     && word(laplacianSchemes_.lookup("default")) != "none"
    )
    {
        laplacianSchemeDefault_ = laplacianSchemes_.lookup("default");
    }


    if (dict.found("fluxRequired"))
    {
        fluxRequired_.merge(dict.subDict("fluxRequired"));

        if
        (
            fluxRequired_.found("default")
         && fluxRequired_.get<word>("default") != "none"
        )
        {
            fluxRequiredDefault_ = fluxRequired_.get<bool>("default");
        }
    }
}


void Foam::fvSchemes::writeDicts(Ostream& os) const
{
    // Write dictionaries
    ddtSchemes_.writeEntry(os);
    d2dt2Schemes_.writeEntry(os);
    interpolationSchemes_.writeEntry(os);
    divSchemes_.writeEntry(os);
    gradSchemes_.writeEntry(os);
    snGradSchemes_.writeEntry(os);
    laplacianSchemes_.writeEntry(os);
    fluxRequired_.writeEntry(os);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvSchemes::fvSchemes
(
    const objectRegistry& obr,
    const word& dictName,
    const dictionary* fallback
)
:
    IOdictionary
    (
        IOobject
        (
            dictName,
            obr.time().system(),
            obr,
            (
                obr.readOpt() == IOobject::MUST_READ
             || obr.readOpt() == IOobject::READ_IF_PRESENT
              ? IOobject::MUST_READ_IF_MODIFIED
              : obr.readOpt()
            ),
            IOobject::NO_WRITE
        ),
        fallback
    ),

    // Named, but empty dictionaries and default schemes
    ddtSchemes_
    (
        objectPath() + ".ddtSchemes"
    ),
    ddtSchemeDefault_
    (
        zero{},
        ddtSchemes_.name() + ".default"
    ),
    d2dt2Schemes_
    (
        objectPath() + ".d2dt2Schemes"
    ),
    d2dt2SchemeDefault_
    (
        zero{},
        d2dt2Schemes_.name() + ".default"
    ),
    interpolationSchemes_
    (
        objectPath() + ".interpolationSchemes"
    ),
    interpolationSchemeDefault_
    (
        zero{},
        interpolationSchemes_.name() + ".default"
    ),
    divSchemes_
    (
        objectPath() + ".divSchemes"
    ),
    divSchemeDefault_
    (
        zero{},
        divSchemes_.name() + ".default"
    ),
    gradSchemes_
    (
        objectPath() + ".gradSchemes"
    ),
    gradSchemeDefault_
    (
        zero{},
        gradSchemes_.name() + ".default"
    ),
    snGradSchemes_
    (
        objectPath() + ".snGradSchemes"
    ),
    snGradSchemeDefault_
    (
        zero{},
        snGradSchemes_.name() + ".default"
    ),
    laplacianSchemes_
    (
        objectPath() + ".laplacianSchemes"
    ),
    laplacianSchemeDefault_
    (
        zero{},
        laplacianSchemes_.name() + ".default"
    ),
    fluxRequired_
    (
        objectPath() + ".fluxRequired"
    ),
    fluxRequiredDefault_(false),
    steady_(false)
{
    if
    (
        readOpt() == IOobject::MUST_READ
     || readOpt() == IOobject::MUST_READ_IF_MODIFIED
     || (readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        read(schemesDict());
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fvSchemes::read()
{
    if (regIOobject::read())
    {
        // Clear current settings except fluxRequired
        clear();

        read(schemesDict());

        return true;
    }

    return false;
}


const Foam::dictionary& Foam::fvSchemes::schemesDict() const
{
    if (found("select"))
    {
        return subDict(word(lookup("select")));
    }

    return *this;
}


Foam::ITstream& Foam::fvSchemes::ddtScheme(const word& name) const
{
    DebugInfo<< "Lookup ddtScheme for " << name << endl;

    if (ddtSchemes_.found(name) || ddtSchemeDefault_.empty())
    {
        return ddtSchemes_.lookup(name);
    }
    else
    {
        // Default scheme
        ITstream& is = const_cast<ITstream&>(ddtSchemeDefault_);
        is.rewind();
        return is;
    }
}


Foam::ITstream& Foam::fvSchemes::d2dt2Scheme(const word& name) const
{
    DebugInfo<< "Lookup d2dt2Scheme for " << name << endl;

    if (d2dt2Schemes_.found(name) || d2dt2SchemeDefault_.empty())
    {
        return d2dt2Schemes_.lookup(name);
    }
    else
    {
        // Default scheme
        ITstream& is = const_cast<ITstream&>(d2dt2SchemeDefault_);
        is.rewind();
        return is;
    }
}


Foam::ITstream& Foam::fvSchemes::interpolationScheme(const word& name) const
{
    DebugInfo<< "Lookup interpolationScheme for " << name << endl;

    if
    (
        interpolationSchemes_.found(name)
     || interpolationSchemeDefault_.empty()
    )
    {
        return interpolationSchemes_.lookup(name);
    }
    else
    {
        // Default scheme
        ITstream& is = const_cast<ITstream&>(interpolationSchemeDefault_);
        is.rewind();
        return is;
    }
}


Foam::ITstream& Foam::fvSchemes::divScheme(const word& name) const
{
    DebugInfo<< "Lookup divScheme for " << name << endl;

    if (divSchemes_.found(name) || divSchemeDefault_.empty())
    {
        return divSchemes_.lookup(name);
    }
    else
    {
        // Default scheme
        ITstream& is = const_cast<ITstream&>(divSchemeDefault_);
        is.rewind();
        return is;
    }
}


Foam::ITstream& Foam::fvSchemes::gradScheme(const word& name) const
{
    DebugInfo<< "Lookup gradScheme for " << name << endl;

    if (gradSchemes_.found(name) || gradSchemeDefault_.empty())
    {
        return gradSchemes_.lookup(name);
    }
    else
    {
        // Default scheme
        ITstream& is = const_cast<ITstream&>(gradSchemeDefault_);
        is.rewind();
        return is;
    }
}


Foam::ITstream& Foam::fvSchemes::snGradScheme(const word& name) const
{
    DebugInfo<< "Lookup snGradScheme for " << name << endl;

    if (snGradSchemes_.found(name) || snGradSchemeDefault_.empty())
    {
        return snGradSchemes_.lookup(name);
    }
    else
    {
        // Default scheme
        ITstream& is = const_cast<ITstream&>(snGradSchemeDefault_);
        is.rewind();
        return is;
    }
}


Foam::ITstream& Foam::fvSchemes::laplacianScheme(const word& name) const
{
    DebugInfo<< "Lookup laplacianScheme for " << name << endl;

    if (laplacianSchemes_.found(name) || laplacianSchemeDefault_.empty())
    {
        return laplacianSchemes_.lookup(name);
    }
    else
    {
        // Default scheme
        ITstream& is = const_cast<ITstream&>(laplacianSchemeDefault_);
        is.rewind();
        return is;
    }
}


void Foam::fvSchemes::setFluxRequired(const word& name) const
{
    DebugInfo<< "Setting fluxRequired for " << name << endl;

    fluxRequired_.add(name, true, true);
}


bool Foam::fvSchemes::fluxRequired(const word& name) const
{
    DebugInfo<< "Lookup fluxRequired for " << name << endl;

    if (fluxRequired_.found(name))
    {
        return true;
    }

    return fluxRequiredDefault_;
}


// ************************************************************************* //
