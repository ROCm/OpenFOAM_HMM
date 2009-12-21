

#ifndef COCO_calcEntryPARSER_H__
#define COCO_calcEntryPARSER_H__

#include "dictionary.H"
#include "scalar.H"
#include "error.H"
#include "wchar.H"


#include "calcEntryScanner.h"

namespace Foam {
namespace functionEntries {
namespace calcEntryInternal {


//! Parser error handing
class Errors {
public:
	int count;      //!< The number of errors detected

	//! Allocate and return a string describing the given error code.
	/** It is the responsibility of the caller to free this string,
	 *  eg, with coco_string_delete()
	 */
	static wchar_t* strerror(int n);

	Errors();               //!< Construct null - start with no errors
	virtual ~Errors();      //!< Destructor
	virtual void clear();   //!< Clear the error count

	//! Handle a general warning 'msg'
	virtual void Warning(const wchar_t* msg);
	//! Handle a general warning 'msg'
	virtual void Warning(int line, int col, const wchar_t* msg);
	//! Handle general error 'msg' (eg, a semantic error)
	virtual void Error(int line, int col, const wchar_t* msg);
	//! Handle syntax error 'n', uses strerror for the message, calls Error()
	virtual void SynErr(int line, int col, int n);
	//! Handle a general exception 'msg'
	virtual void Exception(const wchar_t* msg);

}; // Errors


//! A Coco/R Parser
class Parser {
private:
	enum {
		_EOF=0,
		_ident=1,
		_string=2,
		_variable=3,
		_number=4,
	};
	static const int maxT = 13;

	static const int minErrDist = 2; //!< min. distance before reporting errors

	Token *dummyToken;
	bool deleteErrorsDestruct_; //!< delete the 'errors' member in destructor
	int  errDist;

	void SynErr(int n);         //!< Handle syntax error 'n'
	void Get();
	void Expect(int n);
	bool StartOf(int s);
	void ExpectWeak(int n, int follow);
	bool WeakSeparator(int n, int syFol, int repFol);

public:
	Scanner *scanner;
	Errors  *errors;

	Token *t;                   //!< last recognized token
	Token *la;                  //!< lookahead token

static const int debug = 0;

    //! The parent dictionary
    mutable dictionary* dict_;

    //! The calculation result
    scalar val;

    //! token -> scalar
    scalar getScalar() const
    {
        return coco_string_toDouble(t->val);
    }

    //! token -> string
    std::string getString() const
    {
        char* str = coco_string_create_char(t->val);
        std::string s(str);
        coco_string_delete(str);
        return s;
    }


    //! attach a dictionary
    void dict(const dictionary& dict) const
    {
        dict_ = const_cast<dictionary*>(&dict);
    }


    //! lookup dictionary entry
    scalar getDictLookup() const
    {
        scalar dictValue = 0;

        if (!dict_)
        {
            FatalErrorIn
            (
                "SimpleCalc::getDictEntry() const"
            )   << "No dictionary attached!"
                << exit(FatalError);

            return 0;
        }

        char* chars = coco_string_create_char
        (
            t->val,
            1,
            (coco_string_length(t->val) - 1)
        );
        word keyword(chars);
        coco_string_delete(chars);

        if (debug)
        {
            Info<<"lookup: " << keyword << nl;
        }

        entry* entryPtr = dict_->lookupEntryPtr(keyword, true, false);
        if (entryPtr && !entryPtr->isDict())
        {
            entryPtr->stream() >> dictValue;
        }
        else
        {
            FatalErrorIn
            (
                "SimpleCalc::getDictEntry() const"
            )   << "keyword " << keyword << " is undefined in dictionary "
                << exit(FatalError);
        }


        return dictValue;
    }

    scalar Result() const
    {
        return val;
    }


// * * * * * * * * * * * * * * *  CHARACTERS * * * * * * * * * * * * * * * * //



	//! Construct for the specified scanner
	/**
	 *  Use the default error handling, or optionally provide an error
	 *  handler, which will not be deleted upon destruction.
	 */
	Parser(Scanner* scan, Errors* err = 0);
	~Parser();      //!< Destructor - cleanup errors and dummyToken
	void SemErr(const wchar_t* msg);    //!< Handle semantic error

	void SimpleCalc();
	void Expr(scalar& val);
	void Term(scalar& val);
	void Factor(scalar& val);

	void Parse();                       //!< Execute the parse operation

}; // end Parser

} // namespace
} // namespace
} // namespace


#endif // COCO_calcEntryPARSER_H__

