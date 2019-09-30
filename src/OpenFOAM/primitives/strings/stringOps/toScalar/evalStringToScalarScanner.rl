/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
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

Description
    Ragel lexer interface for lemon grammar of a simple string to
    scalar evaluation

\*---------------------------------------------------------------------------*/

#include "evalStringToScalarScanner.H"
#include "evalStringToScalarDriver.H"
#include "evalStringToScalarLemonParser.h"
#include "evalStringToScalarParser.H"
#include "error.H"

#pragma GCC diagnostic ignored "-Wimplicit-fallthrough"
#pragma GCC diagnostic ignored "-Wold-style-cast"

#ifndef FULLDEBUG
#define NDEBUG
#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::parsing::evalStringToScalar::scanner::debug = 0;


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Ragel lexer with lemon parser integration

// Ragel machine definition
// Ragel variables (p, pe, eof, cs, top, stack, ts, te, act) defined later...
//
// Can use 'variable p xxx;' etc to change these names

%%{
    machine evalScanner;
    write  data;
}%%

#define TOKEN_OF(T)         TOK_##T
#define EMIT_TOKEN(T)                                                         \
    /* Inform driver of last position */                                      \
    driver.parsePosition() = (p-buf);                                         \
    DebugInfo<< "TOKEN_" #T << " at " << driver.parsePosition() << nl;        \
    parser_->parse(TOKEN_OF(T), 0)


%%{
    machine evalScanner;

    action setPosition
    {
        // Inform driver of last position
        driver.parsePosition() = (p-buf);
    }

    action truncated
    {
        // Inform driver of last position
        driver.parsePosition() = 0;
        driver.reportFatal("Truncated input");
    }

    action emit_number {
        // Inform driver of last position
        driver.parsePosition() = (p-buf);

        DebugInfo
            << "NUMBER:" << std::string(ts, te-ts).c_str()
            << " at " << driver.parsePosition() << nl;

        scalar val;

        if (readScalar(std::string(ts, te-ts), val))
        {
            // Emit number
            parser_->parse(TOKEN_OF(NUMBER), val);
        }
        else
        {
            driver.reportFatal("Error reading scalar value");
        }
    }

    decimal = ((digit* '.' digit+) | (digit+ '.'?)) ;
    number  = (digit+ | decimal) ([Ee][\-+]? digit+)? ;
    dnl     = (any* -- '\n') '\n';  # Discard up to and including newline
    lfunc   = space* '(';           # Require functions to have '('

    functions = (
        'pi'         lfunc  @{ fhold; EMIT_TOKEN(PI); }
      | 'degToRad'   lfunc  @{ fhold; EMIT_TOKEN(DEG_TO_RAD); }
      | 'radToDeg'   lfunc  @{ fhold; EMIT_TOKEN(RAD_TO_DEG); }
      | 'exp'        lfunc  @{ fhold; EMIT_TOKEN(EXP); }
      | 'log'        lfunc  @{ fhold; EMIT_TOKEN(LOG); }
      | 'log10'      lfunc  @{ fhold; EMIT_TOKEN(LOG10); }
      | 'pow'        lfunc  @{ fhold; EMIT_TOKEN(POW); }
      | 'sqr'        lfunc  @{ fhold; EMIT_TOKEN(SQR); }
      | 'sqrt'       lfunc  @{ fhold; EMIT_TOKEN(SQRT); }
      | 'cbrt'       lfunc  @{ fhold; EMIT_TOKEN(CBRT); }
      | 'sin'        lfunc  @{ fhold; EMIT_TOKEN(SIN); }
      | 'cos'        lfunc  @{ fhold; EMIT_TOKEN(COS); }
      | 'tan'        lfunc  @{ fhold; EMIT_TOKEN(TAN); }
      | 'asin'       lfunc  @{ fhold; EMIT_TOKEN(ASIN); }
      | 'acos'       lfunc  @{ fhold; EMIT_TOKEN(ACOS); }
      | 'atan'       lfunc  @{ fhold; EMIT_TOKEN(ATAN); }
      | 'atan2'      lfunc  @{ fhold; EMIT_TOKEN(ATAN2); }
      | 'hypot'      lfunc  @{ fhold; EMIT_TOKEN(HYPOT); }
      | 'sinh'       lfunc  @{ fhold; EMIT_TOKEN(SINH); }
      | 'cosh'       lfunc  @{ fhold; EMIT_TOKEN(COSH); }
      | 'tanh'       lfunc  @{ fhold; EMIT_TOKEN(TANH); }
      | 'min'        lfunc  @{ fhold; EMIT_TOKEN(MIN); }
      | 'max'        lfunc  @{ fhold; EMIT_TOKEN(MAX); }
      | 'mag'        lfunc  @{ fhold; EMIT_TOKEN(MAG); }
      | 'magSqr'     lfunc  @{ fhold; EMIT_TOKEN(MAGSQR); }
      | 'rand'       lfunc  @{ fhold; EMIT_TOKEN(RAND); }
    );

    operators = (
        '('  @{ EMIT_TOKEN(LPAREN); }
      | ')'  @{ EMIT_TOKEN(RPAREN); }
      | '+'  @{ EMIT_TOKEN(PLUS); }
      | '-'  @{ EMIT_TOKEN(MINUS); }
      | '*'  @{ EMIT_TOKEN(TIMES); }
      | '/'  @{ EMIT_TOKEN(DIVIDE); }
      | ','  @{ EMIT_TOKEN(COMMA); }
    );


    main := |*
        space*;

        number => emit_number;
        functions;
        operators;

        '/*' any* :>> '*/' @setPosition;        # Multi-line comment
        '//' (any* -- '\n') '\n'* @setPosition; # (sloppy) 1-line comment
        space*;
    *|;
}%%


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::parsing::evalStringToScalar::scanner::~scanner()
{
    if (parser_)
    {
        delete parser_;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool Foam::parsing::evalStringToScalar::scanner::process
(
    const std::string& str,
    parseDriver& driver
)
{
    if (!parser_)
    {
        parser_ = new parser();
    }

    driver.content(str);
    parser_->start(driver);

    // Ragel token start/end (required naming)
    const char* ts;
    const char* te;

    // Local buffer data.
    // - p, pe, eof are required Ragel naming
    // - buf is our own naming

    const char* buf = &(str[0]);
    const char* eof = &(str[str.length()]);
    const char* p = buf;
    const char* pe = eof;

    // Initialize FSM variables
    %%{write init;}%%   /* ^^^ FSM initialization here ^^^ */;

    %%{write exec;}%%  /* ^^^ FSM execution here ^^^ */;

    if (%%{write error;}%% == cs)
    {
        driver.reportFatal("Parse error while scanning", (p-buf));
    }

    if (p != eof)
    {
        driver.reportFatal("Parsing failed with remaining content", (p-buf));
    }

    // Terminate parser execution
    parser_->parse(0, 0);
    parser_->stop();

    return true;
}


// ************************************************************************* //
