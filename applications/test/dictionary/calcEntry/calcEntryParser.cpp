

#include <wchar.h>
#include "calcEntryParser.h"


namespace Foam {
namespace functionEntries {
namespace calcEntryInternal {


// ----------------------------------------------------------------------------
// Parser Implementation
// ----------------------------------------------------------------------------

void Parser::SynErr(int n) {
	if (errDist >= minErrDist) errors->SynErr(la->line, la->col, n);
	errDist = 0;
}


void Parser::SemErr(const wchar_t* msg) {
	if (errDist >= minErrDist) errors->Error(t->line, t->col, msg);
	errDist = 0;
}


void Parser::Get() {
	for (;;) {
		t = la;
		la = scanner->Scan();
		if (la->kind <= maxT) {
			++errDist;
			break;
		}

		if (dummyToken != t) {
			dummyToken->kind = t->kind;
			dummyToken->pos = t->pos;
			dummyToken->col = t->col;
			dummyToken->line = t->line;
			dummyToken->next = NULL;
			coco_string_delete(dummyToken->val);
			dummyToken->val = coco_string_create(t->val);
			t = dummyToken;
		}
		la = t;
	}
}


void Parser::Expect(int n) {
	if (la->kind == n) {
		Get();
	}
	else {
		SynErr(n);
	}
}


void Parser::ExpectWeak(int n, int follow) {
	if (la->kind == n) {
		Get();
	}
	else {
		SynErr(n);
		while (!StartOf(follow)) {
			Get();
		}
	}
}


bool Parser::WeakSeparator(int n, int syFol, int repFol) {
	if (la->kind == n) {
		Get();
		return true;
	}
	else if (StartOf(repFol)) {
		return false;
	}
	else {
		SynErr(n);
		while (!(StartOf(syFol) || StartOf(repFol) || StartOf(0))) {
			Get();
		}
		return StartOf(syFol);
	}
}


void Parser::calcEntry() {
		val = 0; 
		if (la->kind == 5) {
			Get();
			Expr(val);
			Expect(6);
			scanner->buffer->SetPos(t->pos + 1);
			
		} else if (StartOf(1)) {
			Expr(val);
			Expect(0);
		} else SynErr(15);
}

void Parser::Expr(scalar& val) {
		scalar val2 = 0; 
		Term(val);
		while (la->kind == 7 || la->kind == 8) {
			if (la->kind == 7) {
				Get();
				Term(val2);
				val += val2; 
			} else {
				Get();
				Term(val2);
				val -= val2; 
			}
		}
}

void Parser::Term(scalar& val) {
		scalar val2 = 0; 
		Factor(val);
		while (la->kind == 9 || la->kind == 10) {
			if (la->kind == 9) {
				Get();
				Factor(val2);
				val *= val2; 
			} else {
				Get();
				Factor(val2);
				val /= val2; 
			}
		}
}

void Parser::Factor(scalar& val) {
		bool negative = false; 
		if (la->kind == 7 || la->kind == 8) {
			if (la->kind == 7) {
				Get();
			} else {
				Get();
				negative = true; 
			}
		}
		if (la->kind == 1) {
			Func(val);
		} else if (la->kind == 11) {
			Get();
			Expr(val);
			Expect(12);
		} else if (la->kind == 3) {
			Get();
			val = getDictLookup(); 
		} else if (la->kind == 4) {
			Get();
			val = coco_string_toDouble(t->val); 
		} else SynErr(16);
		if (negative) { val = -val; } 
}

void Parser::Func(scalar& val) {
		Expect(1);
		char* str = coco_string_create_char(t->val);
		word funcName(str);
		coco_string_delete(str);
		DynamicList<scalar> param(4);  // hold parameter values
		
		Expect(11);
		if (StartOf(1)) {
			scalar x; 
			Expr(x);
			param.append(x); 
			while (la->kind == 13) {
				Get();
				Expr(x);
				param.append(x); 
			}
		}
		Expect(12);
		val = dispatch(funcName, param); 
}



void Parser::Parse() {
	t = NULL;
	if (dummyToken) {    // safety: someone might call Parse() twice
		delete dummyToken;
	}
	la = dummyToken = new Token();
	la->val = coco_string_create(L"Dummy Token");
	Get();
	calcEntry();
	// let grammar deal with end-of-file expectations

}


Parser::Parser(Scanner* scan, Errors* err)
:
	dummyToken(NULL),
	deleteErrorsDestruct_(!err),
	errDist(minErrDist),
	scanner(scan),
	errors(err),
	t(NULL),
	la(NULL)
{
	if (!errors) {   // add in default error handling
		errors = new Errors();
	}
	// user-defined initialization:
dict_ = 0;
    val = 0;

/*---------------------------------------------------------------------------*/


}


bool Parser::StartOf(int s) {
	const bool T = true;
	const bool x = false;

	static const bool set[2][16] = {
		{T,x,x,x, x,x,x,x, x,x,x,x, x,x,x,x},
		{x,T,x,T, T,x,x,T, T,x,x,T, x,x,x,x}
	};



	return set[s][la->kind];
}


Parser::~Parser() {
	if (deleteErrorsDestruct_) {    // delete default error handling
		delete errors;
	}
	delete dummyToken;
	// user-defined destruction:

}


// ----------------------------------------------------------------------------
// Errors Implementation
// ----------------------------------------------------------------------------

Errors::Errors()
:
	count(0)
{}


Errors::~Errors()
{}


void Errors::clear() {
	count = 0;
}


wchar_t* Errors::strerror(int n)
{
	wchar_t* s;
	switch (n) {
			case 0: s = coco_string_create(L"EOF expected"); break;
			case 1: s = coco_string_create(L"ident expected"); break;
			case 2: s = coco_string_create(L"string expected"); break;
			case 3: s = coco_string_create(L"variable expected"); break;
			case 4: s = coco_string_create(L"number expected"); break;
			case 5: s = coco_string_create(L"\"{\" expected"); break;
			case 6: s = coco_string_create(L"\"}\" expected"); break;
			case 7: s = coco_string_create(L"\"+\" expected"); break;
			case 8: s = coco_string_create(L"\"-\" expected"); break;
			case 9: s = coco_string_create(L"\"*\" expected"); break;
			case 10: s = coco_string_create(L"\"/\" expected"); break;
			case 11: s = coco_string_create(L"\"(\" expected"); break;
			case 12: s = coco_string_create(L"\")\" expected"); break;
			case 13: s = coco_string_create(L"\",\" expected"); break;
			case 14: s = coco_string_create(L"??? expected"); break;
			case 15: s = coco_string_create(L"invalid calcEntry"); break;
			case 16: s = coco_string_create(L"invalid Factor"); break;

		default:
		{
			wchar_t format[20];
			coco_swprintf(format, 20, L"error %d", n);
			s = coco_string_create(format);
		}
		break;
	}
	return s;
}


void Errors::Warning(const wchar_t* msg) {
	wprintf(L"%ls\n", msg);
}


void Errors::Warning(int line, int col, const wchar_t* msg) {
	wprintf(L"-- line %d col %d: %ls\n", line, col, msg);
}


void Errors::Error(int line, int col, const wchar_t* msg) {
	wprintf(L"-- line %d col %d: %ls\n", line, col, msg);
	count++;
}


void Errors::SynErr(int line, int col, int n) {
	wchar_t* msg = this->strerror(n);
	this->Error(line, col, msg);
	coco_string_delete(msg);
}


void Errors::Exception(const wchar_t* msg) {
	wprintf(L"%ls", msg);
	::exit(1);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // namespace
} // namespace
} // namespace


// ************************************************************************* //
