/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2016-2017 Wikki Ltd
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

int Foam::faSchemes::debug(Foam::debug::debugSwitch("faSchemes", false));


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::faSchemes::clear()
{
    ddtSchemes_.clear();
    defaultDdtScheme_.clear();
    d2dt2Schemes_.clear();
    defaultD2dt2Scheme_.clear();
    interpolationSchemes_.clear();
    defaultInterpolationScheme_.clear();
    divSchemes_.clear(); // optional
    defaultDivScheme_.clear();
    gradSchemes_.clear(); // optional
    defaultGradScheme_.clear();
    lnGradSchemes_.clear();
    defaultLnGradScheme_.clear();
    laplacianSchemes_.clear(); // optional
    defaultLaplacianScheme_.clear();
    fluxRequired_.clear();
    defaultFluxRequired_ = false;
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
        defaultDdtScheme_ = ddtSchemes_.lookup("default");

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
        defaultD2dt2Scheme_ = d2dt2Schemes_.lookup("default");
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
        defaultInterpolationScheme_ =
            interpolationSchemes_.lookup("default");
    }


    divSchemes_ = dict.subDict("divSchemes");

    if
    (
        divSchemes_.found("default")
     && word(divSchemes_.lookup("default")) != "none"
    )
    {
        defaultDivScheme_ = divSchemes_.lookup("default");
    }


    gradSchemes_ = dict.subDict("gradSchemes");

    if
    (
        gradSchemes_.found("default")
     && word(gradSchemes_.lookup("default")) != "none"
    )
    {
        defaultGradScheme_ = gradSchemes_.lookup("default");
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
        defaultLnGradScheme_ = lnGradSchemes_.lookup("default");
    }


    laplacianSchemes_ = dict.subDict("laplacianSchemes");

    if
    (
        laplacianSchemes_.found("default")
     && word(laplacianSchemes_.lookup("default")) != "none"
    )
    {
        defaultLaplacianScheme_ = laplacianSchemes_.lookup("default");
    }


    if (dict.found("fluxRequired"))
    {
        fluxRequired_.merge(dict.subDict("fluxRequired"));

        if
        (
            fluxRequired_.found("default")
         && word(fluxRequired_.lookup("default")) != "none"
        )
        {
            defaultFluxRequired_ = Switch(fluxRequired_.lookup("default"));
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faSchemes::faSchemes(const objectRegistry& obr)
:
    IOdictionary
    (
        IOobject
        (
            "faSchemes",
            obr.time().system(),
            obr,
            (
                obr.readOpt() == IOobject::MUST_READ
             || obr.readOpt() == IOobject::READ_IF_PRESENT
              ? IOobject::MUST_READ_IF_MODIFIED
              : obr.readOpt()
            ),
            IOobject::NO_WRITE
        )
    ),
    ddtSchemes_
    (
        ITstream
        (
            objectPath() + ".ddtSchemes",
            tokenList()
        )()
    ),
    defaultDdtScheme_
    (
        ddtSchemes_.name() + ".default",
        tokenList()
    ),
    d2dt2Schemes_
    (
        ITstream
        (
            objectPath() + ".d2dt2Schemes",
            tokenList()
        )()
    ),
    defaultD2dt2Scheme_
    (
        d2dt2Schemes_.name() + ".default",
        tokenList()
    ),
    interpolationSchemes_
    (
        ITstream
        (
            objectPath() + ".interpolationSchemes",
            tokenList()
        )()
    ),
    defaultInterpolationScheme_
    (
        interpolationSchemes_.name() + ".default",
        tokenList()
    ),
    divSchemes_
    (
        ITstream
        (
            objectPath() + ".divSchemes",
            tokenList()
        )()
    ),
    defaultDivScheme_
    (
        divSchemes_.name() + ".default",
        tokenList()
    ),
    gradSchemes_
    (
        ITstream
        (
            objectPath() + ".gradSchemes",
            tokenList()
        )()
    ),
    defaultGradScheme_
    (
        gradSchemes_.name() + ".default",
        tokenList()
    ),
    lnGradSchemes_
    (
        ITstream
        (
            objectPath() + ".snGradSchemes",
            tokenList()
        )()
    ),
    defaultLnGradScheme_
    (
        lnGradSchemes_.name() + ".default",
        tokenList()
    ),
    laplacianSchemes_
    (
        ITstream
        (
            objectPath() + ".laplacianSchemes",
            tokenList()
        )()
    ),
    defaultLaplacianScheme_
    (
        laplacianSchemes_.name() + ".default",
        tokenList()
    ),
    fluxRequired_
    (
        ITstream
        (
            objectPath() + ".fluxRequired",
            tokenList()
        )()
    ),
    defaultFluxRequired_(false)
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
    else
    {
        return false;
    }
}


const Foam::dictionary& Foam::faSchemes::schemesDict() const
{
    if (found("select"))
    {
        return subDict(word(lookup("select")));
    }
    else
    {
        return *this;
    }
}


Foam::ITstream& Foam::faSchemes::ddtScheme(const word& name) const
{
    if (debug)
    {
        Info<< "Lookup ddtScheme for " << name << endl;
    }

    if (ddtSchemes_.found(name) || defaultDdtScheme_.empty())
    {
        return ddtSchemes_.lookup(name);
    }
    else
    {
        const_cast<ITstream&>(defaultDdtScheme_).rewind();
        return const_cast<ITstream&>(defaultDdtScheme_);
    }
}


Foam::ITstream& Foam::faSchemes::d2dt2Scheme(const word& name) const
{
    if (debug)
    {
        Info<< "Lookup d2dt2Scheme for " << name << endl;
    }

    if (d2dt2Schemes_.found(name) || defaultD2dt2Scheme_.empty())
    {
        return d2dt2Schemes_.lookup(name);
    }
    else
    {
        const_cast<ITstream&>(defaultD2dt2Scheme_).rewind();
        return const_cast<ITstream&>(defaultD2dt2Scheme_);
    }
}


Foam::ITstream& Foam::faSchemes::interpolationScheme(const word& name) const
{
    if (debug)
    {
        Info<< "Lookup interpolationScheme for " << name << endl;
    }
     
    if
    (
        interpolationSchemes_.found(name)
     || defaultInterpolationScheme_.empty()
    )
    {
        return interpolationSchemes_.lookup(name);
    }
    else
    {
        const_cast<ITstream&>(defaultInterpolationScheme_).rewind();
        return const_cast<ITstream&>(defaultInterpolationScheme_);
    }
}


Foam::ITstream& Foam::faSchemes::divScheme(const word& name) const
{
    if (debug)
    {
        Info<< "Lookup divScheme for " << name << endl;
    }

    if (divSchemes_.found(name) || defaultDivScheme_.empty())
    {
        return divSchemes_.lookup(name);
    }
    else
    {
        const_cast<ITstream&>(defaultDivScheme_).rewind();
        return const_cast<ITstream&>(defaultDivScheme_);
    }
}


Foam::ITstream& Foam::faSchemes::gradScheme(const word& name) const
{
    if (debug)
    {
        Info<< "Lookup gradScheme for " << name << endl;
    }

    if (gradSchemes_.found(name) || defaultGradScheme_.empty())
    {
        return gradSchemes_.lookup(name);
    }
    else
    {
        const_cast<ITstream&>(defaultGradScheme_).rewind();
        return const_cast<ITstream&>(defaultGradScheme_);
    }
}


Foam::ITstream& Foam::faSchemes::lnGradScheme(const word& name) const
{
    if (debug)
    {
        Info<< "Lookup snGradScheme for " << name << endl;
    }

    if (lnGradSchemes_.found(name) || defaultLnGradScheme_.empty())
    {
        return lnGradSchemes_.lookup(name);
    }
    else
    {
        const_cast<ITstream&>(defaultLnGradScheme_).rewind();
        return const_cast<ITstream&>(defaultLnGradScheme_);
    }
}


Foam::ITstream& Foam::faSchemes::laplacianScheme(const word& name) const
{
    if (debug)
    {
        Info<< "Lookup laplacianScheme for " << name << endl;
    }

    if (laplacianSchemes_.found(name) || defaultLaplacianScheme_.empty())
    {
        return laplacianSchemes_.lookup(name);
    }
    else
    {
        const_cast<ITstream&>(defaultLaplacianScheme_).rewind();
        return const_cast<ITstream&>(defaultLaplacianScheme_);
    }
}


void Foam::faSchemes::setFluxRequired(const word& name) const
{
    if (debug)
    {
        Info<< "Setting fluxRequired for " << name << endl;
    }

    fluxRequired_.add(name, true, true);
}


bool Foam::faSchemes::fluxRequired(const word& name) const
{
    if (debug)
    {
        Info<< "Lookup fluxRequired for " << name << endl;
    }

    if (fluxRequired_.found(name))
    {
        return true;
    }
    else
    {
        return defaultFluxRequired_;
    }
}


bool Foam::faSchemes::writeData(Ostream& os) const
{
    // Write dictionaries
    os << nl << "ddtSchemes";
    ddtSchemes_.write(os, true);

    os << nl << "d2dt2Schemes";
    d2dt2Schemes_.write(os, true);

    os << nl << "interpolationSchemes";
    interpolationSchemes_.write(os, true);

    os << nl << "divSchemes";
    divSchemes_.write(os, true);

    os << nl << "gradSchemes";
    gradSchemes_.write(os, true);

    os << nl << "lnGradSchemes";
    lnGradSchemes_.write(os, true);

    os << nl << "laplacianSchemes";
    laplacianSchemes_.write(os, true);

    os << nl << "fluxRequired";
    fluxRequired_.write(os, true);

    return true;
}


// ************************************************************************* //
