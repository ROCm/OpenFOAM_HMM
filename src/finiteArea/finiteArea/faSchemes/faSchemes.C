/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
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

#include "faSchemes.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::faSchemes::debug(Foam::debug::debugSwitch("faSchemes", 0));


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::faSchemes::clear()
{
    ddtSchemes_.clear();
    ddtSchemeDefault_.clear();

    d2dt2Schemes_.clear();
    d2dt2SchemeDefault_.clear();

    interpolationSchemes_.clear();
    interpolationSchemeDefault_.clear();

    divSchemes_.clear(); // optional
    divSchemeDefault_.clear();

    gradSchemes_.clear(); // optional
    gradSchemeDefault_.clear();

    lnGradSchemes_.clear();
    lnGradSchemeDefault_.clear();

    laplacianSchemes_.clear(); // optional
    laplacianSchemeDefault_.clear();

    fluxRequired_.clear();
    fluxRequiredDefault_ = false;
}


void Foam::faSchemes::read(const dictionary& dict)
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


    if (dict.found("lnGradSchemes"))
    {
        lnGradSchemes_ = dict.subDict("lnGradSchemes");
    }
    else if (!lnGradSchemes_.found("default"))
    {
        lnGradSchemes_.add("default", "corrected");
    }

    if
    (
        lnGradSchemes_.found("default")
     && word(lnGradSchemes_.lookup("default")) != "none"
    )
    {
        lnGradSchemeDefault_ = lnGradSchemes_.lookup("default");
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


void Foam::faSchemes::writeDicts(Ostream& os) const
{
    // Write dictionaries
    ddtSchemes_.writeEntry(os);
    d2dt2Schemes_.writeEntry(os);
    interpolationSchemes_.writeEntry(os);
    divSchemes_.writeEntry(os);
    gradSchemes_.writeEntry(os);
    lnGradSchemes_.writeEntry(os);
    laplacianSchemes_.writeEntry(os);
    fluxRequired_.writeEntry(os);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faSchemes::faSchemes
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
    lnGradSchemes_
    (
        objectPath() + ".lnGradSchemes"
    ),
    lnGradSchemeDefault_
    (
        zero{},
        lnGradSchemes_.name() + ".default"
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
    fluxRequiredDefault_(false)
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

bool Foam::faSchemes::read()
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


const Foam::dictionary& Foam::faSchemes::schemesDict() const
{
    if (found("select"))
    {
        return subDict(word(lookup("select")));
    }

    return *this;
}


Foam::ITstream& Foam::faSchemes::ddtScheme(const word& name) const
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


Foam::ITstream& Foam::faSchemes::d2dt2Scheme(const word& name) const
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


Foam::ITstream& Foam::faSchemes::interpolationScheme(const word& name) const
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


Foam::ITstream& Foam::faSchemes::divScheme(const word& name) const
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


Foam::ITstream& Foam::faSchemes::gradScheme(const word& name) const
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


Foam::ITstream& Foam::faSchemes::lnGradScheme(const word& name) const
{
    DebugInfo<< "Lookup lnGradScheme for " << name << endl;

    if (lnGradSchemes_.found(name) || lnGradSchemeDefault_.empty())
    {
        return lnGradSchemes_.lookup(name);
    }
    else
    {
        // Default scheme
        ITstream& is = const_cast<ITstream&>(lnGradSchemeDefault_);
        is.rewind();
        return is;
    }
}


Foam::ITstream& Foam::faSchemes::laplacianScheme(const word& name) const
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


void Foam::faSchemes::setFluxRequired(const word& name) const
{
    DebugInfo<< "Setting fluxRequired for " << name << endl;

    fluxRequired_.add(name, true, true);
}


bool Foam::faSchemes::fluxRequired(const word& name) const
{
    DebugInfo<< "Lookup fluxRequired for " << name << endl;

    if (fluxRequired_.found(name))
    {
        return true;
    }

    return fluxRequiredDefault_;
}


bool Foam::faSchemes::writeData(Ostream& os) const
{
    this->writeDicts(os);
    return true;
}


// ************************************************************************* //
