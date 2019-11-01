
#line 1 "evalStringToScalarScanner.rl"
/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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
#include "macros.H"

#pragma GCC diagnostic ignored "-Wimplicit-fallthrough"
#pragma GCC diagnostic ignored "-Wold-style-cast"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::parsing::evalStringToScalar::scanner::debug = 0;


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Ragel lexer with lemon parser integration

// Ragel machine definition
// Ragel variables (p, pe, eof, cs, top, stack, ts, te, act) defined later...
//
// Can use 'variable p xxx;' etc to change these names


#line 61 "evalStringToScalarScanner.cc"
static const int evalScanner_start = 4;
static const int evalScanner_first_final = 4;
static const int evalScanner_error = 0;

static const int evalScanner_en_main = 4;


#line 60 "evalStringToScalarScanner.rl"


#define TOKEN_OF(T)         TOK_##T
#define EMIT_TOKEN(T)                                                         \
    driver.parsePosition() = (ts-buf);                                        \
    DebugInfo<< STRINGIFY(T) << ": " << driver.parsePosition() << nl;         \
    parser_->parse(TOKEN_OF(T), 0);                                           \
    driver.parsePosition() = (p-buf);



#line 160 "evalStringToScalarScanner.rl"



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::parsing::evalStringToScalar::scanner::~scanner()
{
    if (parser_)
    {
        delete parser_;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::parsing::evalStringToScalar::scanner::process
(
    const std::string& str,
    size_t strBeg,
    size_t strLen,
    parseDriver& driver
)
{
    if (!parser_)
    {
        parser_ = new parser();
    }

    driver.content(str, strBeg, strLen);

    size_t strEnd = str.length();

    if (strBeg > str.length())
    {
        strBeg = str.length();
    }
    else if (strLen != std::string::npos)
    {
        strLen += strBeg;

        if (strLen < str.length())
        {
            strEnd = strLen;
        }
    }


    parser_->start(driver);

    // Ragel token start/end (required naming)
    const char* ts;
    const char* te;

    // Local buffer data.
    // - p, pe, eof are required Ragel naming
    // - buf is our own naming

    const char* buf = &(str[strBeg]);
    const char* eof = &(str[strEnd]);
    const char* p = buf;
    const char* pe = eof;

    // Initialize FSM variables
    
#line 147 "evalStringToScalarScanner.cc"
	{
	cs = evalScanner_start;
	ts = 0;
	te = 0;
	act = 0;
	}

#line 224 "evalStringToScalarScanner.rl"
   /* ^^^ FSM initialization here ^^^ */;

    
#line 159 "evalStringToScalarScanner.cc"
	{
	if ( p == pe )
		goto _test_eof;
	switch ( cs )
	{
tr2:
#line 73 "evalStringToScalarScanner.rl"
	{{p = ((te))-1;}{
        driver.parsePosition() = (ts-buf);

        DebugInfo
            << "Number:" << std::string(ts, te-ts).c_str()
            << " at " << driver.parsePosition() << nl;

        scalar val(0);

        if (readScalar(std::string(ts, te-ts), val))
        {
            // Emit number
            parser_->parse(TOKEN_OF(NUMBER), val);
        }
        else
        {
            // Range error
            driver.reportFatal("Error parsing number");
        }

        driver.parsePosition() = (p-buf);
    }}
	goto st4;
tr6:
#line 116 "evalStringToScalarScanner.rl"
	{te = p+1;{ EMIT_TOKEN(LPAREN); }}
	goto st4;
tr7:
#line 117 "evalStringToScalarScanner.rl"
	{te = p+1;{ EMIT_TOKEN(RPAREN); }}
	goto st4;
tr8:
#line 120 "evalStringToScalarScanner.rl"
	{te = p+1;{ EMIT_TOKEN(TIMES); }}
	goto st4;
tr9:
#line 118 "evalStringToScalarScanner.rl"
	{te = p+1;{ EMIT_TOKEN(PLUS); }}
	goto st4;
tr10:
#line 122 "evalStringToScalarScanner.rl"
	{te = p+1;{ EMIT_TOKEN(COMMA); }}
	goto st4;
tr11:
#line 119 "evalStringToScalarScanner.rl"
	{te = p+1;{ EMIT_TOKEN(MINUS); }}
	goto st4;
tr13:
#line 121 "evalStringToScalarScanner.rl"
	{te = p+1;{ EMIT_TOKEN(DIVIDE); }}
	goto st4;
tr28:
#line 111 "evalStringToScalarScanner.rl"
	{te = p;p--;}
	goto st4;
tr29:
#line 73 "evalStringToScalarScanner.rl"
	{te = p;p--;{
        driver.parsePosition() = (ts-buf);

        DebugInfo
            << "Number:" << std::string(ts, te-ts).c_str()
            << " at " << driver.parsePosition() << nl;

        scalar val(0);

        if (readScalar(std::string(ts, te-ts), val))
        {
            // Emit number
            parser_->parse(TOKEN_OF(NUMBER), val);
        }
        else
        {
            // Range error
            driver.reportFatal("Error parsing number");
        }

        driver.parsePosition() = (p-buf);
    }}
	goto st4;
tr31:
#line 96 "evalStringToScalarScanner.rl"
	{te = p;p--;{
        driver.parsePosition() = (ts-buf);
        const word ident = word::validate(ts, te);

        driver.reportFatal("Unknown function/type: " + ident);
        driver.parsePosition() = (p-buf);
    }}
	goto st4;
tr33:
#line 1 "NONE"
	{	switch( act ) {
	case 10:
	{{p = ((te))-1;} EMIT_TOKEN(PI); }
	break;
	case 11:
	{{p = ((te))-1;} EMIT_TOKEN(DEG_TO_RAD); }
	break;
	case 12:
	{{p = ((te))-1;} EMIT_TOKEN(RAD_TO_DEG); }
	break;
	case 13:
	{{p = ((te))-1;} EMIT_TOKEN(EXP); }
	break;
	case 15:
	{{p = ((te))-1;} EMIT_TOKEN(LOG10); }
	break;
	case 16:
	{{p = ((te))-1;} EMIT_TOKEN(POW); }
	break;
	case 18:
	{{p = ((te))-1;} EMIT_TOKEN(SQRT); }
	break;
	case 19:
	{{p = ((te))-1;} EMIT_TOKEN(CBRT); }
	break;
	case 23:
	{{p = ((te))-1;} EMIT_TOKEN(ASIN); }
	break;
	case 24:
	{{p = ((te))-1;} EMIT_TOKEN(ACOS); }
	break;
	case 26:
	{{p = ((te))-1;} EMIT_TOKEN(ATAN2); }
	break;
	case 27:
	{{p = ((te))-1;} EMIT_TOKEN(HYPOT); }
	break;
	case 28:
	{{p = ((te))-1;} EMIT_TOKEN(SINH); }
	break;
	case 29:
	{{p = ((te))-1;} EMIT_TOKEN(COSH); }
	break;
	case 30:
	{{p = ((te))-1;} EMIT_TOKEN(TANH); }
	break;
	case 31:
	{{p = ((te))-1;} EMIT_TOKEN(MIN); }
	break;
	case 32:
	{{p = ((te))-1;} EMIT_TOKEN(MAX); }
	break;
	case 34:
	{{p = ((te))-1;} EMIT_TOKEN(MAGSQR); }
	break;
	case 35:
	{{p = ((te))-1;} EMIT_TOKEN(FLOOR); }
	break;
	case 36:
	{{p = ((te))-1;} EMIT_TOKEN(CEIL); }
	break;
	case 37:
	{{p = ((te))-1;} EMIT_TOKEN(ROUND); }
	break;
	case 38:
	{{p = ((te))-1;} p--; EMIT_TOKEN(RAND); }
	break;
	case 39:
	{{p = ((te))-1;}
        driver.parsePosition() = (ts-buf);
        const word ident = word::validate(ts, te);

        driver.reportFatal("Unknown function/type: " + ident);
        driver.parsePosition() = (p-buf);
    }
	break;
	}
	}
	goto st4;
tr43:
#line 140 "evalStringToScalarScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(ATAN); }}
	goto st4;
tr53:
#line 136 "evalStringToScalarScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(COS); }}
	goto st4;
tr74:
#line 129 "evalStringToScalarScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(LOG); }}
	goto st4;
tr81:
#line 148 "evalStringToScalarScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(MAG); }}
	goto st4;
tr105:
#line 135 "evalStringToScalarScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(SIN); }}
	goto st4;
tr108:
#line 132 "evalStringToScalarScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(SQR); }}
	goto st4;
tr112:
#line 137 "evalStringToScalarScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(TAN); }}
	goto st4;
st4:
#line 1 "NONE"
	{ts = 0;}
	if ( ++p == pe )
		goto _test_eof4;
case 4:
#line 1 "NONE"
	{ts = p;}
#line 374 "evalStringToScalarScanner.cc"
	switch( (*p) ) {
		case 32: goto st5;
		case 40: goto tr6;
		case 41: goto tr7;
		case 42: goto tr8;
		case 43: goto tr9;
		case 44: goto tr10;
		case 45: goto tr11;
		case 46: goto st1;
		case 47: goto tr13;
		case 95: goto st9;
		case 97: goto st11;
		case 99: goto st19;
		case 100: goto st26;
		case 101: goto st33;
		case 102: goto st35;
		case 104: goto st39;
		case 108: goto st43;
		case 109: goto st47;
		case 112: goto st53;
		case 114: goto st55;
		case 115: goto st66;
		case 116: goto st71;
	}
	if ( (*p) < 48 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto st5;
	} else if ( (*p) > 57 ) {
		if ( (*p) > 90 ) {
			if ( 98 <= (*p) && (*p) <= 122 )
				goto st9;
		} else if ( (*p) >= 65 )
			goto st9;
	} else
		goto tr14;
	goto st0;
st0:
cs = 0;
	goto _out;
st5:
	if ( ++p == pe )
		goto _test_eof5;
case 5:
	if ( (*p) == 32 )
		goto st5;
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st5;
	goto tr28;
st1:
	if ( ++p == pe )
		goto _test_eof1;
case 1:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto tr0;
	goto st0;
tr0:
#line 1 "NONE"
	{te = p+1;}
	goto st6;
st6:
	if ( ++p == pe )
		goto _test_eof6;
case 6:
#line 438 "evalStringToScalarScanner.cc"
	switch( (*p) ) {
		case 69: goto st2;
		case 101: goto st2;
	}
	if ( 48 <= (*p) && (*p) <= 57 )
		goto tr0;
	goto tr29;
st2:
	if ( ++p == pe )
		goto _test_eof2;
case 2:
	switch( (*p) ) {
		case 43: goto st3;
		case 45: goto st3;
	}
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st7;
	goto tr2;
st3:
	if ( ++p == pe )
		goto _test_eof3;
case 3:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st7;
	goto tr2;
st7:
	if ( ++p == pe )
		goto _test_eof7;
case 7:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st7;
	goto tr29;
tr14:
#line 1 "NONE"
	{te = p+1;}
	goto st8;
st8:
	if ( ++p == pe )
		goto _test_eof8;
case 8:
#line 479 "evalStringToScalarScanner.cc"
	switch( (*p) ) {
		case 46: goto tr0;
		case 69: goto st2;
		case 101: goto st2;
	}
	if ( 48 <= (*p) && (*p) <= 57 )
		goto tr14;
	goto tr29;
st9:
	if ( ++p == pe )
		goto _test_eof9;
case 9:
	if ( (*p) == 95 )
		goto tr32;
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
tr32:
#line 1 "NONE"
	{te = p+1;}
#line 96 "evalStringToScalarScanner.rl"
	{act = 39;}
	goto st10;
tr38:
#line 1 "NONE"
	{te = p+1;}
#line 139 "evalStringToScalarScanner.rl"
	{act = 24;}
	goto st10;
tr40:
#line 1 "NONE"
	{te = p+1;}
#line 138 "evalStringToScalarScanner.rl"
	{act = 23;}
	goto st10;
tr44:
#line 1 "NONE"
	{te = p+1;}
#line 141 "evalStringToScalarScanner.rl"
	{act = 26;}
	goto st10;
tr49:
#line 1 "NONE"
	{te = p+1;}
#line 134 "evalStringToScalarScanner.rl"
	{act = 19;}
	goto st10;
tr51:
#line 1 "NONE"
	{te = p+1;}
#line 151 "evalStringToScalarScanner.rl"
	{act = 36;}
	goto st10;
tr54:
#line 1 "NONE"
	{te = p+1;}
#line 144 "evalStringToScalarScanner.rl"
	{act = 29;}
	goto st10;
tr61:
#line 1 "NONE"
	{te = p+1;}
#line 126 "evalStringToScalarScanner.rl"
	{act = 11;}
	goto st10;
tr63:
#line 1 "NONE"
	{te = p+1;}
#line 128 "evalStringToScalarScanner.rl"
	{act = 13;}
	goto st10;
tr67:
#line 1 "NONE"
	{te = p+1;}
#line 150 "evalStringToScalarScanner.rl"
	{act = 35;}
	goto st10;
tr71:
#line 1 "NONE"
	{te = p+1;}
#line 142 "evalStringToScalarScanner.rl"
	{act = 27;}
	goto st10;
tr76:
#line 1 "NONE"
	{te = p+1;}
#line 130 "evalStringToScalarScanner.rl"
	{act = 15;}
	goto st10;
tr80:
#line 1 "NONE"
	{te = p+1;}
#line 147 "evalStringToScalarScanner.rl"
	{act = 32;}
	goto st10;
tr84:
#line 1 "NONE"
	{te = p+1;}
#line 149 "evalStringToScalarScanner.rl"
	{act = 34;}
	goto st10;
tr85:
#line 1 "NONE"
	{te = p+1;}
#line 146 "evalStringToScalarScanner.rl"
	{act = 31;}
	goto st10;
tr86:
#line 1 "NONE"
	{te = p+1;}
#line 125 "evalStringToScalarScanner.rl"
	{act = 10;}
	goto st10;
tr88:
#line 1 "NONE"
	{te = p+1;}
#line 131 "evalStringToScalarScanner.rl"
	{act = 16;}
	goto st10;
tr97:
#line 1 "NONE"
	{te = p+1;}
#line 127 "evalStringToScalarScanner.rl"
	{act = 12;}
	goto st10;
tr98:
#line 1 "NONE"
	{te = p+1;}
#line 153 "evalStringToScalarScanner.rl"
	{act = 38;}
	goto st10;
tr101:
#line 1 "NONE"
	{te = p+1;}
#line 152 "evalStringToScalarScanner.rl"
	{act = 37;}
	goto st10;
tr106:
#line 1 "NONE"
	{te = p+1;}
#line 143 "evalStringToScalarScanner.rl"
	{act = 28;}
	goto st10;
tr109:
#line 1 "NONE"
	{te = p+1;}
#line 133 "evalStringToScalarScanner.rl"
	{act = 18;}
	goto st10;
tr113:
#line 1 "NONE"
	{te = p+1;}
#line 145 "evalStringToScalarScanner.rl"
	{act = 30;}
	goto st10;
st10:
	if ( ++p == pe )
		goto _test_eof10;
case 10:
#line 645 "evalStringToScalarScanner.cc"
	if ( (*p) == 95 )
		goto tr32;
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr33;
st11:
	if ( ++p == pe )
		goto _test_eof11;
case 11:
	switch( (*p) ) {
		case 95: goto tr32;
		case 99: goto st12;
		case 115: goto st14;
		case 116: goto st16;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st12:
	if ( ++p == pe )
		goto _test_eof12;
case 12:
	switch( (*p) ) {
		case 95: goto tr32;
		case 111: goto st13;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st13:
	if ( ++p == pe )
		goto _test_eof13;
case 13:
	switch( (*p) ) {
		case 95: goto tr32;
		case 115: goto tr38;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st14:
	if ( ++p == pe )
		goto _test_eof14;
case 14:
	switch( (*p) ) {
		case 95: goto tr32;
		case 105: goto st15;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st15:
	if ( ++p == pe )
		goto _test_eof15;
case 15:
	switch( (*p) ) {
		case 95: goto tr32;
		case 110: goto tr40;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st16:
	if ( ++p == pe )
		goto _test_eof16;
case 16:
	switch( (*p) ) {
		case 95: goto tr32;
		case 97: goto st17;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 98 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st17:
	if ( ++p == pe )
		goto _test_eof17;
case 17:
	switch( (*p) ) {
		case 95: goto tr32;
		case 110: goto st18;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st18:
	if ( ++p == pe )
		goto _test_eof18;
case 18:
	switch( (*p) ) {
		case 50: goto tr44;
		case 95: goto tr32;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr43;
st19:
	if ( ++p == pe )
		goto _test_eof19;
case 19:
	switch( (*p) ) {
		case 95: goto tr32;
		case 98: goto st20;
		case 101: goto st22;
		case 111: goto st24;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st20:
	if ( ++p == pe )
		goto _test_eof20;
case 20:
	switch( (*p) ) {
		case 95: goto tr32;
		case 114: goto st21;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st21:
	if ( ++p == pe )
		goto _test_eof21;
case 21:
	switch( (*p) ) {
		case 95: goto tr32;
		case 116: goto tr49;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st22:
	if ( ++p == pe )
		goto _test_eof22;
case 22:
	switch( (*p) ) {
		case 95: goto tr32;
		case 105: goto st23;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st23:
	if ( ++p == pe )
		goto _test_eof23;
case 23:
	switch( (*p) ) {
		case 95: goto tr32;
		case 108: goto tr51;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st24:
	if ( ++p == pe )
		goto _test_eof24;
case 24:
	switch( (*p) ) {
		case 95: goto tr32;
		case 115: goto st25;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st25:
	if ( ++p == pe )
		goto _test_eof25;
case 25:
	switch( (*p) ) {
		case 95: goto tr32;
		case 104: goto tr54;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr53;
st26:
	if ( ++p == pe )
		goto _test_eof26;
case 26:
	switch( (*p) ) {
		case 95: goto tr32;
		case 101: goto st27;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st27:
	if ( ++p == pe )
		goto _test_eof27;
case 27:
	switch( (*p) ) {
		case 95: goto tr32;
		case 103: goto st28;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st28:
	if ( ++p == pe )
		goto _test_eof28;
case 28:
	switch( (*p) ) {
		case 84: goto st29;
		case 95: goto tr32;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st29:
	if ( ++p == pe )
		goto _test_eof29;
case 29:
	switch( (*p) ) {
		case 95: goto tr32;
		case 111: goto st30;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st30:
	if ( ++p == pe )
		goto _test_eof30;
case 30:
	switch( (*p) ) {
		case 82: goto st31;
		case 95: goto tr32;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st31:
	if ( ++p == pe )
		goto _test_eof31;
case 31:
	switch( (*p) ) {
		case 95: goto tr32;
		case 97: goto st32;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 98 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st32:
	if ( ++p == pe )
		goto _test_eof32;
case 32:
	switch( (*p) ) {
		case 95: goto tr32;
		case 100: goto tr61;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st33:
	if ( ++p == pe )
		goto _test_eof33;
case 33:
	switch( (*p) ) {
		case 95: goto tr32;
		case 120: goto st34;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st34:
	if ( ++p == pe )
		goto _test_eof34;
case 34:
	switch( (*p) ) {
		case 95: goto tr32;
		case 112: goto tr63;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st35:
	if ( ++p == pe )
		goto _test_eof35;
case 35:
	switch( (*p) ) {
		case 95: goto tr32;
		case 108: goto st36;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st36:
	if ( ++p == pe )
		goto _test_eof36;
case 36:
	switch( (*p) ) {
		case 95: goto tr32;
		case 111: goto st37;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st37:
	if ( ++p == pe )
		goto _test_eof37;
case 37:
	switch( (*p) ) {
		case 95: goto tr32;
		case 111: goto st38;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st38:
	if ( ++p == pe )
		goto _test_eof38;
case 38:
	switch( (*p) ) {
		case 95: goto tr32;
		case 114: goto tr67;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st39:
	if ( ++p == pe )
		goto _test_eof39;
case 39:
	switch( (*p) ) {
		case 95: goto tr32;
		case 121: goto st40;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st40:
	if ( ++p == pe )
		goto _test_eof40;
case 40:
	switch( (*p) ) {
		case 95: goto tr32;
		case 112: goto st41;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st41:
	if ( ++p == pe )
		goto _test_eof41;
case 41:
	switch( (*p) ) {
		case 95: goto tr32;
		case 111: goto st42;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st42:
	if ( ++p == pe )
		goto _test_eof42;
case 42:
	switch( (*p) ) {
		case 95: goto tr32;
		case 116: goto tr71;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st43:
	if ( ++p == pe )
		goto _test_eof43;
case 43:
	switch( (*p) ) {
		case 95: goto tr32;
		case 111: goto st44;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st44:
	if ( ++p == pe )
		goto _test_eof44;
case 44:
	switch( (*p) ) {
		case 95: goto tr32;
		case 103: goto st45;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st45:
	if ( ++p == pe )
		goto _test_eof45;
case 45:
	switch( (*p) ) {
		case 49: goto st46;
		case 95: goto tr32;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr74;
st46:
	if ( ++p == pe )
		goto _test_eof46;
case 46:
	switch( (*p) ) {
		case 48: goto tr76;
		case 95: goto tr32;
	}
	if ( (*p) < 65 ) {
		if ( 49 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st47:
	if ( ++p == pe )
		goto _test_eof47;
case 47:
	switch( (*p) ) {
		case 95: goto tr32;
		case 97: goto st48;
		case 105: goto st52;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 98 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st48:
	if ( ++p == pe )
		goto _test_eof48;
case 48:
	switch( (*p) ) {
		case 95: goto tr32;
		case 103: goto st49;
		case 120: goto tr80;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st49:
	if ( ++p == pe )
		goto _test_eof49;
case 49:
	switch( (*p) ) {
		case 83: goto st50;
		case 95: goto tr32;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr81;
st50:
	if ( ++p == pe )
		goto _test_eof50;
case 50:
	switch( (*p) ) {
		case 95: goto tr32;
		case 113: goto st51;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st51:
	if ( ++p == pe )
		goto _test_eof51;
case 51:
	switch( (*p) ) {
		case 95: goto tr32;
		case 114: goto tr84;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st52:
	if ( ++p == pe )
		goto _test_eof52;
case 52:
	switch( (*p) ) {
		case 95: goto tr32;
		case 110: goto tr85;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st53:
	if ( ++p == pe )
		goto _test_eof53;
case 53:
	switch( (*p) ) {
		case 95: goto tr32;
		case 105: goto tr86;
		case 111: goto st54;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st54:
	if ( ++p == pe )
		goto _test_eof54;
case 54:
	switch( (*p) ) {
		case 95: goto tr32;
		case 119: goto tr88;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st55:
	if ( ++p == pe )
		goto _test_eof55;
case 55:
	switch( (*p) ) {
		case 95: goto tr32;
		case 97: goto st56;
		case 111: goto st63;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 98 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st56:
	if ( ++p == pe )
		goto _test_eof56;
case 56:
	switch( (*p) ) {
		case 95: goto tr32;
		case 100: goto st57;
		case 110: goto st62;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st57:
	if ( ++p == pe )
		goto _test_eof57;
case 57:
	switch( (*p) ) {
		case 84: goto st58;
		case 95: goto tr32;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st58:
	if ( ++p == pe )
		goto _test_eof58;
case 58:
	switch( (*p) ) {
		case 95: goto tr32;
		case 111: goto st59;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st59:
	if ( ++p == pe )
		goto _test_eof59;
case 59:
	switch( (*p) ) {
		case 68: goto st60;
		case 95: goto tr32;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st60:
	if ( ++p == pe )
		goto _test_eof60;
case 60:
	switch( (*p) ) {
		case 95: goto tr32;
		case 101: goto st61;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st61:
	if ( ++p == pe )
		goto _test_eof61;
case 61:
	switch( (*p) ) {
		case 95: goto tr32;
		case 103: goto tr97;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st62:
	if ( ++p == pe )
		goto _test_eof62;
case 62:
	switch( (*p) ) {
		case 95: goto tr32;
		case 100: goto tr98;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st63:
	if ( ++p == pe )
		goto _test_eof63;
case 63:
	switch( (*p) ) {
		case 95: goto tr32;
		case 117: goto st64;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st64:
	if ( ++p == pe )
		goto _test_eof64;
case 64:
	switch( (*p) ) {
		case 95: goto tr32;
		case 110: goto st65;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st65:
	if ( ++p == pe )
		goto _test_eof65;
case 65:
	switch( (*p) ) {
		case 95: goto tr32;
		case 100: goto tr101;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st66:
	if ( ++p == pe )
		goto _test_eof66;
case 66:
	switch( (*p) ) {
		case 95: goto tr32;
		case 105: goto st67;
		case 113: goto st69;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st67:
	if ( ++p == pe )
		goto _test_eof67;
case 67:
	switch( (*p) ) {
		case 95: goto tr32;
		case 110: goto st68;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st68:
	if ( ++p == pe )
		goto _test_eof68;
case 68:
	switch( (*p) ) {
		case 95: goto tr32;
		case 104: goto tr106;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr105;
st69:
	if ( ++p == pe )
		goto _test_eof69;
case 69:
	switch( (*p) ) {
		case 95: goto tr32;
		case 114: goto st70;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st70:
	if ( ++p == pe )
		goto _test_eof70;
case 70:
	switch( (*p) ) {
		case 95: goto tr32;
		case 116: goto tr109;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr108;
st71:
	if ( ++p == pe )
		goto _test_eof71;
case 71:
	switch( (*p) ) {
		case 95: goto tr32;
		case 97: goto st72;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 98 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st72:
	if ( ++p == pe )
		goto _test_eof72;
case 72:
	switch( (*p) ) {
		case 95: goto tr32;
		case 110: goto st73;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr31;
st73:
	if ( ++p == pe )
		goto _test_eof73;
case 73:
	switch( (*p) ) {
		case 95: goto tr32;
		case 104: goto tr113;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr32;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr32;
	} else
		goto tr32;
	goto tr112;
	}
	_test_eof4: cs = 4; goto _test_eof; 
	_test_eof5: cs = 5; goto _test_eof; 
	_test_eof1: cs = 1; goto _test_eof; 
	_test_eof6: cs = 6; goto _test_eof; 
	_test_eof2: cs = 2; goto _test_eof; 
	_test_eof3: cs = 3; goto _test_eof; 
	_test_eof7: cs = 7; goto _test_eof; 
	_test_eof8: cs = 8; goto _test_eof; 
	_test_eof9: cs = 9; goto _test_eof; 
	_test_eof10: cs = 10; goto _test_eof; 
	_test_eof11: cs = 11; goto _test_eof; 
	_test_eof12: cs = 12; goto _test_eof; 
	_test_eof13: cs = 13; goto _test_eof; 
	_test_eof14: cs = 14; goto _test_eof; 
	_test_eof15: cs = 15; goto _test_eof; 
	_test_eof16: cs = 16; goto _test_eof; 
	_test_eof17: cs = 17; goto _test_eof; 
	_test_eof18: cs = 18; goto _test_eof; 
	_test_eof19: cs = 19; goto _test_eof; 
	_test_eof20: cs = 20; goto _test_eof; 
	_test_eof21: cs = 21; goto _test_eof; 
	_test_eof22: cs = 22; goto _test_eof; 
	_test_eof23: cs = 23; goto _test_eof; 
	_test_eof24: cs = 24; goto _test_eof; 
	_test_eof25: cs = 25; goto _test_eof; 
	_test_eof26: cs = 26; goto _test_eof; 
	_test_eof27: cs = 27; goto _test_eof; 
	_test_eof28: cs = 28; goto _test_eof; 
	_test_eof29: cs = 29; goto _test_eof; 
	_test_eof30: cs = 30; goto _test_eof; 
	_test_eof31: cs = 31; goto _test_eof; 
	_test_eof32: cs = 32; goto _test_eof; 
	_test_eof33: cs = 33; goto _test_eof; 
	_test_eof34: cs = 34; goto _test_eof; 
	_test_eof35: cs = 35; goto _test_eof; 
	_test_eof36: cs = 36; goto _test_eof; 
	_test_eof37: cs = 37; goto _test_eof; 
	_test_eof38: cs = 38; goto _test_eof; 
	_test_eof39: cs = 39; goto _test_eof; 
	_test_eof40: cs = 40; goto _test_eof; 
	_test_eof41: cs = 41; goto _test_eof; 
	_test_eof42: cs = 42; goto _test_eof; 
	_test_eof43: cs = 43; goto _test_eof; 
	_test_eof44: cs = 44; goto _test_eof; 
	_test_eof45: cs = 45; goto _test_eof; 
	_test_eof46: cs = 46; goto _test_eof; 
	_test_eof47: cs = 47; goto _test_eof; 
	_test_eof48: cs = 48; goto _test_eof; 
	_test_eof49: cs = 49; goto _test_eof; 
	_test_eof50: cs = 50; goto _test_eof; 
	_test_eof51: cs = 51; goto _test_eof; 
	_test_eof52: cs = 52; goto _test_eof; 
	_test_eof53: cs = 53; goto _test_eof; 
	_test_eof54: cs = 54; goto _test_eof; 
	_test_eof55: cs = 55; goto _test_eof; 
	_test_eof56: cs = 56; goto _test_eof; 
	_test_eof57: cs = 57; goto _test_eof; 
	_test_eof58: cs = 58; goto _test_eof; 
	_test_eof59: cs = 59; goto _test_eof; 
	_test_eof60: cs = 60; goto _test_eof; 
	_test_eof61: cs = 61; goto _test_eof; 
	_test_eof62: cs = 62; goto _test_eof; 
	_test_eof63: cs = 63; goto _test_eof; 
	_test_eof64: cs = 64; goto _test_eof; 
	_test_eof65: cs = 65; goto _test_eof; 
	_test_eof66: cs = 66; goto _test_eof; 
	_test_eof67: cs = 67; goto _test_eof; 
	_test_eof68: cs = 68; goto _test_eof; 
	_test_eof69: cs = 69; goto _test_eof; 
	_test_eof70: cs = 70; goto _test_eof; 
	_test_eof71: cs = 71; goto _test_eof; 
	_test_eof72: cs = 72; goto _test_eof; 
	_test_eof73: cs = 73; goto _test_eof; 

	_test_eof: {}
	if ( p == eof )
	{
	switch ( cs ) {
	case 5: goto tr28;
	case 6: goto tr29;
	case 2: goto tr2;
	case 3: goto tr2;
	case 7: goto tr29;
	case 8: goto tr29;
	case 9: goto tr31;
	case 10: goto tr33;
	case 11: goto tr31;
	case 12: goto tr31;
	case 13: goto tr31;
	case 14: goto tr31;
	case 15: goto tr31;
	case 16: goto tr31;
	case 17: goto tr31;
	case 18: goto tr43;
	case 19: goto tr31;
	case 20: goto tr31;
	case 21: goto tr31;
	case 22: goto tr31;
	case 23: goto tr31;
	case 24: goto tr31;
	case 25: goto tr53;
	case 26: goto tr31;
	case 27: goto tr31;
	case 28: goto tr31;
	case 29: goto tr31;
	case 30: goto tr31;
	case 31: goto tr31;
	case 32: goto tr31;
	case 33: goto tr31;
	case 34: goto tr31;
	case 35: goto tr31;
	case 36: goto tr31;
	case 37: goto tr31;
	case 38: goto tr31;
	case 39: goto tr31;
	case 40: goto tr31;
	case 41: goto tr31;
	case 42: goto tr31;
	case 43: goto tr31;
	case 44: goto tr31;
	case 45: goto tr74;
	case 46: goto tr31;
	case 47: goto tr31;
	case 48: goto tr31;
	case 49: goto tr81;
	case 50: goto tr31;
	case 51: goto tr31;
	case 52: goto tr31;
	case 53: goto tr31;
	case 54: goto tr31;
	case 55: goto tr31;
	case 56: goto tr31;
	case 57: goto tr31;
	case 58: goto tr31;
	case 59: goto tr31;
	case 60: goto tr31;
	case 61: goto tr31;
	case 62: goto tr31;
	case 63: goto tr31;
	case 64: goto tr31;
	case 65: goto tr31;
	case 66: goto tr31;
	case 67: goto tr31;
	case 68: goto tr105;
	case 69: goto tr31;
	case 70: goto tr108;
	case 71: goto tr31;
	case 72: goto tr31;
	case 73: goto tr112;
	}
	}

	_out: {}
	}

#line 226 "evalStringToScalarScanner.rl"
  /* ^^^ FSM execution here ^^^ */;

    if (0 == cs)
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
