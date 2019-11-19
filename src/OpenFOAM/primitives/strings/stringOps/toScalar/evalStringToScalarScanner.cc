
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


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Ragel lexer with lemon parser integration

// Ragel machine definition
// Ragel variables (p, pe, eof, cs, top, stack, ts, te, act) defined later...
//
// Can use 'variable p xxx;' etc to change these names


#line 56 "evalStringToScalarScanner.cc"
static const int evalScanner_start = 7;
static const int evalScanner_first_final = 7;
static const int evalScanner_error = 0;

static const int evalScanner_en_main = 7;


#line 55 "evalStringToScalarScanner.rl"



#undef DebugScannerInfo
#define DebugScannerInfo if (debug & 4) InfoErr


#define TOKEN_OF(T)         TOK_##T
#define EMIT_TOKEN(T)                                                         \
    driver.parsePosition() = (ts-buf);                                        \
    DebugScannerInfo << STRINGIFY(T) << ": "<< driver.parsePosition() << nl;  \
    parser_->parse(TOKEN_OF(T), 0);                                           \
    driver.parsePosition() = (p-buf);



#line 181 "evalStringToScalarScanner.rl"



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

#line 245 "evalStringToScalarScanner.rl"
   /* ^^^ FSM initialization here ^^^ */;

    
#line 159 "evalStringToScalarScanner.cc"
	{
	if ( p == pe )
		goto _test_eof;
	switch ( cs )
	{
tr0:
#line 134 "evalStringToScalarScanner.rl"
	{te = p+1;{ EMIT_TOKEN(LAND); }}
	goto st7;
tr3:
#line 73 "evalStringToScalarScanner.rl"
	{{p = ((te))-1;}{
        driver.parsePosition() = (ts-buf);

        DebugScannerInfo
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
	goto st7;
tr6:
#line 132 "evalStringToScalarScanner.rl"
	{te = p+1;{ EMIT_TOKEN(EQUAL); }}
	goto st7;
tr7:
#line 135 "evalStringToScalarScanner.rl"
	{te = p+1;{ EMIT_TOKEN(LOR); }}
	goto st7;
tr10:
#line 117 "evalStringToScalarScanner.rl"
	{te = p+1;{ EMIT_TOKEN(PERCENT); }}
	goto st7;
tr12:
#line 118 "evalStringToScalarScanner.rl"
	{te = p+1;{ EMIT_TOKEN(LPAREN); }}
	goto st7;
tr13:
#line 119 "evalStringToScalarScanner.rl"
	{te = p+1;{ EMIT_TOKEN(RPAREN); }}
	goto st7;
tr14:
#line 120 "evalStringToScalarScanner.rl"
	{te = p+1;{ EMIT_TOKEN(TIMES); }}
	goto st7;
tr15:
#line 121 "evalStringToScalarScanner.rl"
	{te = p+1;{ EMIT_TOKEN(PLUS); }}
	goto st7;
tr16:
#line 123 "evalStringToScalarScanner.rl"
	{te = p+1;{ EMIT_TOKEN(COMMA); }}
	goto st7;
tr17:
#line 122 "evalStringToScalarScanner.rl"
	{te = p+1;{ EMIT_TOKEN(MINUS); }}
	goto st7;
tr19:
#line 124 "evalStringToScalarScanner.rl"
	{te = p+1;{ EMIT_TOKEN(DIVIDE); }}
	goto st7;
tr21:
#line 127 "evalStringToScalarScanner.rl"
	{te = p+1;{ EMIT_TOKEN(COLON); }}
	goto st7;
tr25:
#line 126 "evalStringToScalarScanner.rl"
	{te = p+1;{ EMIT_TOKEN(QUESTION); }}
	goto st7;
tr41:
#line 111 "evalStringToScalarScanner.rl"
	{te = p;p--;}
	goto st7;
tr42:
#line 116 "evalStringToScalarScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(NOT); }}
	goto st7;
tr43:
#line 133 "evalStringToScalarScanner.rl"
	{te = p+1;{ EMIT_TOKEN(NOT_EQUAL); }}
	goto st7;
tr44:
#line 73 "evalStringToScalarScanner.rl"
	{te = p;p--;{
        driver.parsePosition() = (ts-buf);

        DebugScannerInfo
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
	goto st7;
tr46:
#line 128 "evalStringToScalarScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(LESS); }}
	goto st7;
tr47:
#line 129 "evalStringToScalarScanner.rl"
	{te = p+1;{ EMIT_TOKEN(LESS_EQ); }}
	goto st7;
tr48:
#line 130 "evalStringToScalarScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(GREATER); }}
	goto st7;
tr49:
#line 131 "evalStringToScalarScanner.rl"
	{te = p+1;{ EMIT_TOKEN(GREATER_EQ); }}
	goto st7;
tr50:
#line 96 "evalStringToScalarScanner.rl"
	{te = p;p--;{
        driver.parsePosition() = (ts-buf);
        const word ident = word::validate(ts, te);

        driver.reportFatal("Unknown function/type: " + ident);
        driver.parsePosition() = (p-buf);
    }}
	goto st7;
tr52:
#line 1 "NONE"
	{	switch( act ) {
	case 23:
	{{p = ((te))-1;} EMIT_TOKEN(PI); }
	break;
	case 24:
	{{p = ((te))-1;} EMIT_TOKEN(DEG_TO_RAD); }
	break;
	case 25:
	{{p = ((te))-1;} EMIT_TOKEN(RAD_TO_DEG); }
	break;
	case 26:
	{{p = ((te))-1;} EMIT_TOKEN(EXP); }
	break;
	case 28:
	{{p = ((te))-1;} EMIT_TOKEN(LOG10); }
	break;
	case 29:
	{{p = ((te))-1;} EMIT_TOKEN(POW); }
	break;
	case 31:
	{{p = ((te))-1;} EMIT_TOKEN(SQRT); }
	break;
	case 32:
	{{p = ((te))-1;} EMIT_TOKEN(CBRT); }
	break;
	case 36:
	{{p = ((te))-1;} EMIT_TOKEN(ASIN); }
	break;
	case 37:
	{{p = ((te))-1;} EMIT_TOKEN(ACOS); }
	break;
	case 39:
	{{p = ((te))-1;} EMIT_TOKEN(ATAN2); }
	break;
	case 40:
	{{p = ((te))-1;} EMIT_TOKEN(HYPOT); }
	break;
	case 41:
	{{p = ((te))-1;} EMIT_TOKEN(SINH); }
	break;
	case 42:
	{{p = ((te))-1;} EMIT_TOKEN(COSH); }
	break;
	case 43:
	{{p = ((te))-1;} EMIT_TOKEN(TANH); }
	break;
	case 44:
	{{p = ((te))-1;} EMIT_TOKEN(MIN); }
	break;
	case 45:
	{{p = ((te))-1;} EMIT_TOKEN(MAX); }
	break;
	case 47:
	{{p = ((te))-1;} EMIT_TOKEN(MAGSQR); }
	break;
	case 48:
	{{p = ((te))-1;} EMIT_TOKEN(FLOOR); }
	break;
	case 49:
	{{p = ((te))-1;} EMIT_TOKEN(CEIL); }
	break;
	case 50:
	{{p = ((te))-1;} EMIT_TOKEN(ROUND); }
	break;
	case 51:
	{{p = ((te))-1;} EMIT_TOKEN(RAND); }
	break;
	case 52:
	{{p = ((te))-1;} EMIT_TOKEN(BOOL); }
	break;
	case 53:
	{{p = ((te))-1;} EMIT_TOKEN(BOOL_FALSE); }
	break;
	case 54:
	{{p = ((te))-1;} EMIT_TOKEN(BOOL_TRUE); }
	break;
	case 55:
	{{p = ((te))-1;}
        driver.parsePosition() = (ts-buf);
        const word ident = word::validate(ts, te);

        driver.reportFatal("Unknown function/type: " + ident);
        driver.parsePosition() = (p-buf);
    }
	break;
	}
	}
	goto st7;
tr62:
#line 156 "evalStringToScalarScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(ATAN); }}
	goto st7;
tr75:
#line 152 "evalStringToScalarScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(COS); }}
	goto st7;
tr100:
#line 145 "evalStringToScalarScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(LOG); }}
	goto st7;
tr107:
#line 164 "evalStringToScalarScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(MAG); }}
	goto st7;
tr131:
#line 151 "evalStringToScalarScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(SIN); }}
	goto st7;
tr134:
#line 148 "evalStringToScalarScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(SQR); }}
	goto st7;
tr139:
#line 153 "evalStringToScalarScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(TAN); }}
	goto st7;
st7:
#line 1 "NONE"
	{ts = 0;}
	if ( ++p == pe )
		goto _test_eof7;
case 7:
#line 1 "NONE"
	{ts = p;}
#line 431 "evalStringToScalarScanner.cc"
	switch( (*p) ) {
		case 32: goto st8;
		case 33: goto st9;
		case 37: goto tr10;
		case 38: goto st1;
		case 40: goto tr12;
		case 41: goto tr13;
		case 42: goto tr14;
		case 43: goto tr15;
		case 44: goto tr16;
		case 45: goto tr17;
		case 46: goto st2;
		case 47: goto tr19;
		case 58: goto tr21;
		case 60: goto st13;
		case 61: goto st5;
		case 62: goto st14;
		case 63: goto tr25;
		case 95: goto st15;
		case 97: goto st17;
		case 98: goto st25;
		case 99: goto st28;
		case 100: goto st35;
		case 101: goto st42;
		case 102: goto st44;
		case 104: goto st51;
		case 108: goto st55;
		case 109: goto st59;
		case 112: goto st65;
		case 114: goto st67;
		case 115: goto st78;
		case 116: goto st83;
		case 124: goto st6;
	}
	if ( (*p) < 48 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto st8;
	} else if ( (*p) > 57 ) {
		if ( (*p) > 90 ) {
			if ( 103 <= (*p) && (*p) <= 122 )
				goto st15;
		} else if ( (*p) >= 65 )
			goto st15;
	} else
		goto tr20;
	goto st0;
st0:
cs = 0;
	goto _out;
st8:
	if ( ++p == pe )
		goto _test_eof8;
case 8:
	if ( (*p) == 32 )
		goto st8;
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st8;
	goto tr41;
st9:
	if ( ++p == pe )
		goto _test_eof9;
case 9:
	if ( (*p) == 61 )
		goto tr43;
	goto tr42;
st1:
	if ( ++p == pe )
		goto _test_eof1;
case 1:
	if ( (*p) == 38 )
		goto tr0;
	goto st0;
st2:
	if ( ++p == pe )
		goto _test_eof2;
case 2:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto tr2;
	goto st0;
tr2:
#line 1 "NONE"
	{te = p+1;}
	goto st10;
st10:
	if ( ++p == pe )
		goto _test_eof10;
case 10:
#line 519 "evalStringToScalarScanner.cc"
	switch( (*p) ) {
		case 69: goto st3;
		case 101: goto st3;
	}
	if ( 48 <= (*p) && (*p) <= 57 )
		goto tr2;
	goto tr44;
st3:
	if ( ++p == pe )
		goto _test_eof3;
case 3:
	switch( (*p) ) {
		case 43: goto st4;
		case 45: goto st4;
	}
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st11;
	goto tr3;
st4:
	if ( ++p == pe )
		goto _test_eof4;
case 4:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st11;
	goto tr3;
st11:
	if ( ++p == pe )
		goto _test_eof11;
case 11:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st11;
	goto tr44;
tr20:
#line 1 "NONE"
	{te = p+1;}
	goto st12;
st12:
	if ( ++p == pe )
		goto _test_eof12;
case 12:
#line 560 "evalStringToScalarScanner.cc"
	switch( (*p) ) {
		case 46: goto tr2;
		case 69: goto st3;
		case 101: goto st3;
	}
	if ( 48 <= (*p) && (*p) <= 57 )
		goto tr20;
	goto tr44;
st13:
	if ( ++p == pe )
		goto _test_eof13;
case 13:
	if ( (*p) == 61 )
		goto tr47;
	goto tr46;
st5:
	if ( ++p == pe )
		goto _test_eof5;
case 5:
	if ( (*p) == 61 )
		goto tr6;
	goto st0;
st14:
	if ( ++p == pe )
		goto _test_eof14;
case 14:
	if ( (*p) == 61 )
		goto tr49;
	goto tr48;
st15:
	if ( ++p == pe )
		goto _test_eof15;
case 15:
	if ( (*p) == 95 )
		goto tr51;
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
tr51:
#line 1 "NONE"
	{te = p+1;}
#line 96 "evalStringToScalarScanner.rl"
	{act = 55;}
	goto st16;
tr57:
#line 1 "NONE"
	{te = p+1;}
#line 155 "evalStringToScalarScanner.rl"
	{act = 37;}
	goto st16;
tr59:
#line 1 "NONE"
	{te = p+1;}
#line 154 "evalStringToScalarScanner.rl"
	{act = 36;}
	goto st16;
tr63:
#line 1 "NONE"
	{te = p+1;}
#line 157 "evalStringToScalarScanner.rl"
	{act = 39;}
	goto st16;
tr66:
#line 1 "NONE"
	{te = p+1;}
#line 170 "evalStringToScalarScanner.rl"
	{act = 52;}
	goto st16;
tr71:
#line 1 "NONE"
	{te = p+1;}
#line 150 "evalStringToScalarScanner.rl"
	{act = 32;}
	goto st16;
tr73:
#line 1 "NONE"
	{te = p+1;}
#line 167 "evalStringToScalarScanner.rl"
	{act = 49;}
	goto st16;
tr76:
#line 1 "NONE"
	{te = p+1;}
#line 160 "evalStringToScalarScanner.rl"
	{act = 42;}
	goto st16;
tr83:
#line 1 "NONE"
	{te = p+1;}
#line 142 "evalStringToScalarScanner.rl"
	{act = 24;}
	goto st16;
tr85:
#line 1 "NONE"
	{te = p+1;}
#line 144 "evalStringToScalarScanner.rl"
	{act = 26;}
	goto st16;
tr90:
#line 1 "NONE"
	{te = p+1;}
#line 173 "evalStringToScalarScanner.rl"
	{act = 53;}
	goto st16;
tr93:
#line 1 "NONE"
	{te = p+1;}
#line 166 "evalStringToScalarScanner.rl"
	{act = 48;}
	goto st16;
tr97:
#line 1 "NONE"
	{te = p+1;}
#line 158 "evalStringToScalarScanner.rl"
	{act = 40;}
	goto st16;
tr102:
#line 1 "NONE"
	{te = p+1;}
#line 146 "evalStringToScalarScanner.rl"
	{act = 28;}
	goto st16;
tr106:
#line 1 "NONE"
	{te = p+1;}
#line 163 "evalStringToScalarScanner.rl"
	{act = 45;}
	goto st16;
tr110:
#line 1 "NONE"
	{te = p+1;}
#line 165 "evalStringToScalarScanner.rl"
	{act = 47;}
	goto st16;
tr111:
#line 1 "NONE"
	{te = p+1;}
#line 162 "evalStringToScalarScanner.rl"
	{act = 44;}
	goto st16;
tr112:
#line 1 "NONE"
	{te = p+1;}
#line 141 "evalStringToScalarScanner.rl"
	{act = 23;}
	goto st16;
tr114:
#line 1 "NONE"
	{te = p+1;}
#line 147 "evalStringToScalarScanner.rl"
	{act = 29;}
	goto st16;
tr123:
#line 1 "NONE"
	{te = p+1;}
#line 143 "evalStringToScalarScanner.rl"
	{act = 25;}
	goto st16;
tr124:
#line 1 "NONE"
	{te = p+1;}
#line 169 "evalStringToScalarScanner.rl"
	{act = 51;}
	goto st16;
tr127:
#line 1 "NONE"
	{te = p+1;}
#line 168 "evalStringToScalarScanner.rl"
	{act = 50;}
	goto st16;
tr132:
#line 1 "NONE"
	{te = p+1;}
#line 159 "evalStringToScalarScanner.rl"
	{act = 41;}
	goto st16;
tr135:
#line 1 "NONE"
	{te = p+1;}
#line 149 "evalStringToScalarScanner.rl"
	{act = 31;}
	goto st16;
tr140:
#line 1 "NONE"
	{te = p+1;}
#line 161 "evalStringToScalarScanner.rl"
	{act = 43;}
	goto st16;
tr142:
#line 1 "NONE"
	{te = p+1;}
#line 174 "evalStringToScalarScanner.rl"
	{act = 54;}
	goto st16;
st16:
	if ( ++p == pe )
		goto _test_eof16;
case 16:
#line 765 "evalStringToScalarScanner.cc"
	if ( (*p) == 95 )
		goto tr51;
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr52;
st17:
	if ( ++p == pe )
		goto _test_eof17;
case 17:
	switch( (*p) ) {
		case 95: goto tr51;
		case 99: goto st18;
		case 115: goto st20;
		case 116: goto st22;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st18:
	if ( ++p == pe )
		goto _test_eof18;
case 18:
	switch( (*p) ) {
		case 95: goto tr51;
		case 111: goto st19;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st19:
	if ( ++p == pe )
		goto _test_eof19;
case 19:
	switch( (*p) ) {
		case 95: goto tr51;
		case 115: goto tr57;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st20:
	if ( ++p == pe )
		goto _test_eof20;
case 20:
	switch( (*p) ) {
		case 95: goto tr51;
		case 105: goto st21;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st21:
	if ( ++p == pe )
		goto _test_eof21;
case 21:
	switch( (*p) ) {
		case 95: goto tr51;
		case 110: goto tr59;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st22:
	if ( ++p == pe )
		goto _test_eof22;
case 22:
	switch( (*p) ) {
		case 95: goto tr51;
		case 97: goto st23;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 98 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st23:
	if ( ++p == pe )
		goto _test_eof23;
case 23:
	switch( (*p) ) {
		case 95: goto tr51;
		case 110: goto st24;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st24:
	if ( ++p == pe )
		goto _test_eof24;
case 24:
	switch( (*p) ) {
		case 50: goto tr63;
		case 95: goto tr51;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr62;
st25:
	if ( ++p == pe )
		goto _test_eof25;
case 25:
	switch( (*p) ) {
		case 95: goto tr51;
		case 111: goto st26;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st26:
	if ( ++p == pe )
		goto _test_eof26;
case 26:
	switch( (*p) ) {
		case 95: goto tr51;
		case 111: goto st27;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st27:
	if ( ++p == pe )
		goto _test_eof27;
case 27:
	switch( (*p) ) {
		case 95: goto tr51;
		case 108: goto tr66;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st28:
	if ( ++p == pe )
		goto _test_eof28;
case 28:
	switch( (*p) ) {
		case 95: goto tr51;
		case 98: goto st29;
		case 101: goto st31;
		case 111: goto st33;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st29:
	if ( ++p == pe )
		goto _test_eof29;
case 29:
	switch( (*p) ) {
		case 95: goto tr51;
		case 114: goto st30;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st30:
	if ( ++p == pe )
		goto _test_eof30;
case 30:
	switch( (*p) ) {
		case 95: goto tr51;
		case 116: goto tr71;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st31:
	if ( ++p == pe )
		goto _test_eof31;
case 31:
	switch( (*p) ) {
		case 95: goto tr51;
		case 105: goto st32;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st32:
	if ( ++p == pe )
		goto _test_eof32;
case 32:
	switch( (*p) ) {
		case 95: goto tr51;
		case 108: goto tr73;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st33:
	if ( ++p == pe )
		goto _test_eof33;
case 33:
	switch( (*p) ) {
		case 95: goto tr51;
		case 115: goto st34;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st34:
	if ( ++p == pe )
		goto _test_eof34;
case 34:
	switch( (*p) ) {
		case 95: goto tr51;
		case 104: goto tr76;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr75;
st35:
	if ( ++p == pe )
		goto _test_eof35;
case 35:
	switch( (*p) ) {
		case 95: goto tr51;
		case 101: goto st36;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st36:
	if ( ++p == pe )
		goto _test_eof36;
case 36:
	switch( (*p) ) {
		case 95: goto tr51;
		case 103: goto st37;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st37:
	if ( ++p == pe )
		goto _test_eof37;
case 37:
	switch( (*p) ) {
		case 84: goto st38;
		case 95: goto tr51;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st38:
	if ( ++p == pe )
		goto _test_eof38;
case 38:
	switch( (*p) ) {
		case 95: goto tr51;
		case 111: goto st39;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st39:
	if ( ++p == pe )
		goto _test_eof39;
case 39:
	switch( (*p) ) {
		case 82: goto st40;
		case 95: goto tr51;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st40:
	if ( ++p == pe )
		goto _test_eof40;
case 40:
	switch( (*p) ) {
		case 95: goto tr51;
		case 97: goto st41;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 98 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st41:
	if ( ++p == pe )
		goto _test_eof41;
case 41:
	switch( (*p) ) {
		case 95: goto tr51;
		case 100: goto tr83;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st42:
	if ( ++p == pe )
		goto _test_eof42;
case 42:
	switch( (*p) ) {
		case 95: goto tr51;
		case 120: goto st43;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st43:
	if ( ++p == pe )
		goto _test_eof43;
case 43:
	switch( (*p) ) {
		case 95: goto tr51;
		case 112: goto tr85;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st44:
	if ( ++p == pe )
		goto _test_eof44;
case 44:
	switch( (*p) ) {
		case 95: goto tr51;
		case 97: goto st45;
		case 108: goto st48;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 98 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st45:
	if ( ++p == pe )
		goto _test_eof45;
case 45:
	switch( (*p) ) {
		case 95: goto tr51;
		case 108: goto st46;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st46:
	if ( ++p == pe )
		goto _test_eof46;
case 46:
	switch( (*p) ) {
		case 95: goto tr51;
		case 115: goto st47;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st47:
	if ( ++p == pe )
		goto _test_eof47;
case 47:
	switch( (*p) ) {
		case 95: goto tr51;
		case 101: goto tr90;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st48:
	if ( ++p == pe )
		goto _test_eof48;
case 48:
	switch( (*p) ) {
		case 95: goto tr51;
		case 111: goto st49;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st49:
	if ( ++p == pe )
		goto _test_eof49;
case 49:
	switch( (*p) ) {
		case 95: goto tr51;
		case 111: goto st50;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st50:
	if ( ++p == pe )
		goto _test_eof50;
case 50:
	switch( (*p) ) {
		case 95: goto tr51;
		case 114: goto tr93;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st51:
	if ( ++p == pe )
		goto _test_eof51;
case 51:
	switch( (*p) ) {
		case 95: goto tr51;
		case 121: goto st52;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st52:
	if ( ++p == pe )
		goto _test_eof52;
case 52:
	switch( (*p) ) {
		case 95: goto tr51;
		case 112: goto st53;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st53:
	if ( ++p == pe )
		goto _test_eof53;
case 53:
	switch( (*p) ) {
		case 95: goto tr51;
		case 111: goto st54;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st54:
	if ( ++p == pe )
		goto _test_eof54;
case 54:
	switch( (*p) ) {
		case 95: goto tr51;
		case 116: goto tr97;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st55:
	if ( ++p == pe )
		goto _test_eof55;
case 55:
	switch( (*p) ) {
		case 95: goto tr51;
		case 111: goto st56;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st56:
	if ( ++p == pe )
		goto _test_eof56;
case 56:
	switch( (*p) ) {
		case 95: goto tr51;
		case 103: goto st57;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st57:
	if ( ++p == pe )
		goto _test_eof57;
case 57:
	switch( (*p) ) {
		case 49: goto st58;
		case 95: goto tr51;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr100;
st58:
	if ( ++p == pe )
		goto _test_eof58;
case 58:
	switch( (*p) ) {
		case 48: goto tr102;
		case 95: goto tr51;
	}
	if ( (*p) < 65 ) {
		if ( 49 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st59:
	if ( ++p == pe )
		goto _test_eof59;
case 59:
	switch( (*p) ) {
		case 95: goto tr51;
		case 97: goto st60;
		case 105: goto st64;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 98 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st60:
	if ( ++p == pe )
		goto _test_eof60;
case 60:
	switch( (*p) ) {
		case 95: goto tr51;
		case 103: goto st61;
		case 120: goto tr106;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st61:
	if ( ++p == pe )
		goto _test_eof61;
case 61:
	switch( (*p) ) {
		case 83: goto st62;
		case 95: goto tr51;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr107;
st62:
	if ( ++p == pe )
		goto _test_eof62;
case 62:
	switch( (*p) ) {
		case 95: goto tr51;
		case 113: goto st63;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st63:
	if ( ++p == pe )
		goto _test_eof63;
case 63:
	switch( (*p) ) {
		case 95: goto tr51;
		case 114: goto tr110;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st64:
	if ( ++p == pe )
		goto _test_eof64;
case 64:
	switch( (*p) ) {
		case 95: goto tr51;
		case 110: goto tr111;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st65:
	if ( ++p == pe )
		goto _test_eof65;
case 65:
	switch( (*p) ) {
		case 95: goto tr51;
		case 105: goto tr112;
		case 111: goto st66;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st66:
	if ( ++p == pe )
		goto _test_eof66;
case 66:
	switch( (*p) ) {
		case 95: goto tr51;
		case 119: goto tr114;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st67:
	if ( ++p == pe )
		goto _test_eof67;
case 67:
	switch( (*p) ) {
		case 95: goto tr51;
		case 97: goto st68;
		case 111: goto st75;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 98 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st68:
	if ( ++p == pe )
		goto _test_eof68;
case 68:
	switch( (*p) ) {
		case 95: goto tr51;
		case 100: goto st69;
		case 110: goto st74;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st69:
	if ( ++p == pe )
		goto _test_eof69;
case 69:
	switch( (*p) ) {
		case 84: goto st70;
		case 95: goto tr51;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st70:
	if ( ++p == pe )
		goto _test_eof70;
case 70:
	switch( (*p) ) {
		case 95: goto tr51;
		case 111: goto st71;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st71:
	if ( ++p == pe )
		goto _test_eof71;
case 71:
	switch( (*p) ) {
		case 68: goto st72;
		case 95: goto tr51;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st72:
	if ( ++p == pe )
		goto _test_eof72;
case 72:
	switch( (*p) ) {
		case 95: goto tr51;
		case 101: goto st73;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st73:
	if ( ++p == pe )
		goto _test_eof73;
case 73:
	switch( (*p) ) {
		case 95: goto tr51;
		case 103: goto tr123;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st74:
	if ( ++p == pe )
		goto _test_eof74;
case 74:
	switch( (*p) ) {
		case 95: goto tr51;
		case 100: goto tr124;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st75:
	if ( ++p == pe )
		goto _test_eof75;
case 75:
	switch( (*p) ) {
		case 95: goto tr51;
		case 117: goto st76;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st76:
	if ( ++p == pe )
		goto _test_eof76;
case 76:
	switch( (*p) ) {
		case 95: goto tr51;
		case 110: goto st77;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st77:
	if ( ++p == pe )
		goto _test_eof77;
case 77:
	switch( (*p) ) {
		case 95: goto tr51;
		case 100: goto tr127;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st78:
	if ( ++p == pe )
		goto _test_eof78;
case 78:
	switch( (*p) ) {
		case 95: goto tr51;
		case 105: goto st79;
		case 113: goto st81;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st79:
	if ( ++p == pe )
		goto _test_eof79;
case 79:
	switch( (*p) ) {
		case 95: goto tr51;
		case 110: goto st80;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st80:
	if ( ++p == pe )
		goto _test_eof80;
case 80:
	switch( (*p) ) {
		case 95: goto tr51;
		case 104: goto tr132;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr131;
st81:
	if ( ++p == pe )
		goto _test_eof81;
case 81:
	switch( (*p) ) {
		case 95: goto tr51;
		case 114: goto st82;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st82:
	if ( ++p == pe )
		goto _test_eof82;
case 82:
	switch( (*p) ) {
		case 95: goto tr51;
		case 116: goto tr135;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr134;
st83:
	if ( ++p == pe )
		goto _test_eof83;
case 83:
	switch( (*p) ) {
		case 95: goto tr51;
		case 97: goto st84;
		case 114: goto st86;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 98 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st84:
	if ( ++p == pe )
		goto _test_eof84;
case 84:
	switch( (*p) ) {
		case 95: goto tr51;
		case 110: goto st85;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st85:
	if ( ++p == pe )
		goto _test_eof85;
case 85:
	switch( (*p) ) {
		case 95: goto tr51;
		case 104: goto tr140;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr139;
st86:
	if ( ++p == pe )
		goto _test_eof86;
case 86:
	switch( (*p) ) {
		case 95: goto tr51;
		case 117: goto st87;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st87:
	if ( ++p == pe )
		goto _test_eof87;
case 87:
	switch( (*p) ) {
		case 95: goto tr51;
		case 101: goto tr142;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr51;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr51;
	} else
		goto tr51;
	goto tr50;
st6:
	if ( ++p == pe )
		goto _test_eof6;
case 6:
	if ( (*p) == 124 )
		goto tr7;
	goto st0;
	}
	_test_eof7: cs = 7; goto _test_eof; 
	_test_eof8: cs = 8; goto _test_eof; 
	_test_eof9: cs = 9; goto _test_eof; 
	_test_eof1: cs = 1; goto _test_eof; 
	_test_eof2: cs = 2; goto _test_eof; 
	_test_eof10: cs = 10; goto _test_eof; 
	_test_eof3: cs = 3; goto _test_eof; 
	_test_eof4: cs = 4; goto _test_eof; 
	_test_eof11: cs = 11; goto _test_eof; 
	_test_eof12: cs = 12; goto _test_eof; 
	_test_eof13: cs = 13; goto _test_eof; 
	_test_eof5: cs = 5; goto _test_eof; 
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
	_test_eof74: cs = 74; goto _test_eof; 
	_test_eof75: cs = 75; goto _test_eof; 
	_test_eof76: cs = 76; goto _test_eof; 
	_test_eof77: cs = 77; goto _test_eof; 
	_test_eof78: cs = 78; goto _test_eof; 
	_test_eof79: cs = 79; goto _test_eof; 
	_test_eof80: cs = 80; goto _test_eof; 
	_test_eof81: cs = 81; goto _test_eof; 
	_test_eof82: cs = 82; goto _test_eof; 
	_test_eof83: cs = 83; goto _test_eof; 
	_test_eof84: cs = 84; goto _test_eof; 
	_test_eof85: cs = 85; goto _test_eof; 
	_test_eof86: cs = 86; goto _test_eof; 
	_test_eof87: cs = 87; goto _test_eof; 
	_test_eof6: cs = 6; goto _test_eof; 

	_test_eof: {}
	if ( p == eof )
	{
	switch ( cs ) {
	case 8: goto tr41;
	case 9: goto tr42;
	case 10: goto tr44;
	case 3: goto tr3;
	case 4: goto tr3;
	case 11: goto tr44;
	case 12: goto tr44;
	case 13: goto tr46;
	case 14: goto tr48;
	case 15: goto tr50;
	case 16: goto tr52;
	case 17: goto tr50;
	case 18: goto tr50;
	case 19: goto tr50;
	case 20: goto tr50;
	case 21: goto tr50;
	case 22: goto tr50;
	case 23: goto tr50;
	case 24: goto tr62;
	case 25: goto tr50;
	case 26: goto tr50;
	case 27: goto tr50;
	case 28: goto tr50;
	case 29: goto tr50;
	case 30: goto tr50;
	case 31: goto tr50;
	case 32: goto tr50;
	case 33: goto tr50;
	case 34: goto tr75;
	case 35: goto tr50;
	case 36: goto tr50;
	case 37: goto tr50;
	case 38: goto tr50;
	case 39: goto tr50;
	case 40: goto tr50;
	case 41: goto tr50;
	case 42: goto tr50;
	case 43: goto tr50;
	case 44: goto tr50;
	case 45: goto tr50;
	case 46: goto tr50;
	case 47: goto tr50;
	case 48: goto tr50;
	case 49: goto tr50;
	case 50: goto tr50;
	case 51: goto tr50;
	case 52: goto tr50;
	case 53: goto tr50;
	case 54: goto tr50;
	case 55: goto tr50;
	case 56: goto tr50;
	case 57: goto tr100;
	case 58: goto tr50;
	case 59: goto tr50;
	case 60: goto tr50;
	case 61: goto tr107;
	case 62: goto tr50;
	case 63: goto tr50;
	case 64: goto tr50;
	case 65: goto tr50;
	case 66: goto tr50;
	case 67: goto tr50;
	case 68: goto tr50;
	case 69: goto tr50;
	case 70: goto tr50;
	case 71: goto tr50;
	case 72: goto tr50;
	case 73: goto tr50;
	case 74: goto tr50;
	case 75: goto tr50;
	case 76: goto tr50;
	case 77: goto tr50;
	case 78: goto tr50;
	case 79: goto tr50;
	case 80: goto tr131;
	case 81: goto tr50;
	case 82: goto tr134;
	case 83: goto tr50;
	case 84: goto tr50;
	case 85: goto tr139;
	case 86: goto tr50;
	case 87: goto tr50;
	}
	}

	_out: {}
	}

#line 247 "evalStringToScalarScanner.rl"
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
