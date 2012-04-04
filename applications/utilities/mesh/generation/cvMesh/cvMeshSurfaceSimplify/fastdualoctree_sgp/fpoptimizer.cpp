#line 1 "fpoptimizer/fpoptimizer_header.txt"
/***************************************************************************\
|* Function Parser for C++ v4.0.3                                          *|
|*-------------------------------------------------------------------------*|
|* Function optimizer                                                      *|
|*-------------------------------------------------------------------------*|
|* Copyright: Joel Yliluoma                                                *|
\***************************************************************************/

/* NOTE:
   This is a concatenation of all the header and source files of the
   original optimizer source code. All the code has been concatenated
   into this single file for convenience of usage (in other words, to
   simply use the optimizer, it's enough to add this file to the project
   rather than a multitude of files which the original optimizer source
   code is composed of).

   Thus this file exists for the usage of the Function parser library
   only, and is not suitable for developing it further. If you want to
   develop the library further, you should download the development
   version of the library, which has all the original source files.
 */

#include "fpconfig.h"
#ifdef FP_SUPPORT_OPTIMIZER

#define FPOPTIMIZER_MERGED_FILE

#line 1 "fpoptimizer/fpoptimizer_hash.h"
#ifndef FPoptimizerHashH
#define FPoptimizerHashH

#ifdef _MSC_VER

typedef unsigned long long fphash_value_t;
#define FPHASH_CONST(x) x##ULL

#else

#include <stdint.h>
typedef uint_fast64_t fphash_value_t;
#define FPHASH_CONST(x) x##ULL

#endif

namespace FUNCTIONPARSERTYPES
{
    struct fphash_t
    {
        fphash_value_t hash1, hash2;

        bool operator==(const fphash_t& rhs) const
        { return hash1 == rhs.hash1 && hash2 == rhs.hash2; }

        bool operator!=(const fphash_t& rhs) const
        { return hash1 != rhs.hash1 || hash2 != rhs.hash2; }

        bool operator<(const fphash_t& rhs) const
        { return hash1 != rhs.hash1 ? hash1 < rhs.hash1 : hash2 < rhs.hash2; }
    };
}

#endif

#line 1 "fpoptimizer/fpoptimizer_autoptr.h"
#ifndef FPOptimizerAutoPtrH
#define FPOptimizerAutoPtrH

template<typename Ref>
class FPOPT_autoptr
{
public:
    FPOPT_autoptr()                   : p(0)   { }
    FPOPT_autoptr(Ref*        b) : p(b)   { Birth(); }
    FPOPT_autoptr(const FPOPT_autoptr& b) : p(b.p) { Birth(); }

    inline Ref& operator* () const { return *p; }
    inline Ref* operator->() const { return p; }

    FPOPT_autoptr& operator= (Ref*        b) { Set(b); return *this; }
    FPOPT_autoptr& operator= (const FPOPT_autoptr& b) { Set(b.p); return *this; }
#ifdef __GXX_EXPERIMENTAL_CXX0X__
    FPOPT_autoptr(FPOPT_autoptr&& b)      : p(b.p) { b.p = 0; }
    FPOPT_autoptr& operator= (FPOPT_autoptr&& b) { if(p != b.p) { Forget(); p=b.p; b.p=0; }
                                                   return *this; }
#endif

    ~FPOPT_autoptr() { Forget(); }

    void UnsafeSetP(Ref* newp) { p = newp; }
    void swap(FPOPT_autoptr<Ref>& b) { Ref* tmp=p; p=b.p; b.p=tmp; }

private:
    inline static void Have(Ref* p2);
    inline void Forget();
    inline void Birth();
    inline void Set(Ref* p2);
private:
    Ref* p;
};

//
template<typename Ref>
inline void FPOPT_autoptr<Ref>::Forget()
{
    if(!p) return;
    p->RefCount -= 1;
    if(!p->RefCount) delete p;
    //assert(p->RefCount >= 0);
}
template<typename Ref>
inline void FPOPT_autoptr<Ref>::Have(Ref* p2)
{
    if(p2) ++(p2->RefCount);
}
template<typename Ref>
inline void FPOPT_autoptr<Ref>::Birth()
{
    Have(p);
}
template<typename Ref>
inline void FPOPT_autoptr<Ref>::Set(Ref* p2)
{
    Have(p2);
    Forget();
    p = p2;
}

#endif

#line 1 "fpoptimizer/fpoptimizer_codetree.h"
#ifndef FPOptimizer_CodeTreeH
#define FPOptimizer_CodeTreeH

#include "fpconfig.h"
#include "fparser.h"
#include "fptypes.h"

#ifdef FP_SUPPORT_OPTIMIZER

#include <vector>
#include <utility>

// line removed for fpoptimizer.c: #include "fpoptimizer_hash.h"
// line removed for fpoptimizer.c: #include "fpoptimizer_autoptr.h"

#ifdef FP_EPSILON
 #define NEGATIVE_MAXIMUM (-FP_EPSILON)
#else
 #define NEGATIVE_MAXIMUM (-1e-14)
#endif

namespace FPoptimizer_Grammar
{
    struct Grammar;
}

namespace FPoptimizer_ByteCode
{
    class ByteCodeSynth;
}

namespace FPoptimizer_CodeTree
{
    class CodeTree;

    struct MinMaxTree
    {
        double min,max;
        bool has_min, has_max;
        MinMaxTree() : min(),max(),has_min(false),has_max(false) { }
        MinMaxTree(double mi,double ma): min(mi),max(ma),has_min(true),has_max(true) { }
        MinMaxTree(bool,double ma): min(),max(ma),has_min(false),has_max(true) { }
        MinMaxTree(double mi,bool): min(mi),max(),has_min(true),has_max(false) { }
    };

    struct CodeTreeData;
    class CodeTree
    {
        typedef FPOPT_autoptr<CodeTreeData> DataP;
        DataP data;

    public:
    public:
        CodeTree();
        ~CodeTree();

        explicit CodeTree(double v); // produce an immed
        struct VarTag { };
        explicit CodeTree(unsigned varno, VarTag); // produce a var reference
        struct CloneTag { };
        explicit CodeTree(const CodeTree& b, CloneTag);

        /* Generates a CodeTree from the given bytecode */
        void GenerateFrom(
            const std::vector<unsigned>& byteCode,
            const std::vector<double>& immed,
            const FunctionParser::Data& data,
            bool keep_powi = false);

        void GenerateFrom(
            const std::vector<unsigned>& byteCode,
            const std::vector<double>& immed,
            const FunctionParser::Data& data,
            const std::vector<CodeTree>& var_trees,
            bool keep_powi = false);

        void SynthesizeByteCode(
            std::vector<unsigned>& byteCode,
            std::vector<double>&   immed,
            size_t& stacktop_max);

        void SynthesizeByteCode(
            FPoptimizer_ByteCode::ByteCodeSynth& synth,
            bool MustPopTemps=true) const;

        size_t SynthCommonSubExpressions(
            FPoptimizer_ByteCode::ByteCodeSynth& synth) const;

        void SetParams(const std::vector<CodeTree>& RefParams);
        void SetParamsMove(std::vector<CodeTree>& RefParams);

        CodeTree GetUniqueRef();
        // ^use this when CodeTree tmp=x; tmp.CopyOnWrite(); does not do exactly what you want

#ifdef __GXX_EXPERIMENTAL_CXX0X__
        void SetParams(std::vector<CodeTree>&& RefParams);
#endif
        void SetParam(size_t which, const CodeTree& b);
        void SetParamMove(size_t which, CodeTree& b);
        void AddParam(const CodeTree& param);
        void AddParamMove(CodeTree& param);
        void AddParams(const std::vector<CodeTree>& RefParams);
        void AddParamsMove(std::vector<CodeTree>& RefParams);
        void AddParamsMove(std::vector<CodeTree>& RefParams, size_t replacing_slot);
        void DelParam(size_t index);
        void DelParams();

        void Become(const CodeTree& b);

        inline size_t GetParamCount() const { return GetParams().size(); }
        inline CodeTree& GetParam(size_t n) { return GetParams()[n]; }
        inline const CodeTree& GetParam(size_t n) const { return GetParams()[n]; }
        inline void SetOpcode(FUNCTIONPARSERTYPES::OPCODE o);
        inline void SetFuncOpcode(FUNCTIONPARSERTYPES::OPCODE o, unsigned f);
        inline void SetVar(unsigned v);
        inline void SetImmed(double v);
        inline FUNCTIONPARSERTYPES::OPCODE GetOpcode() const;
        inline FUNCTIONPARSERTYPES::fphash_t GetHash() const;
        inline const std::vector<CodeTree>& GetParams() const;
        inline std::vector<CodeTree>& GetParams();
        inline size_t GetDepth() const;
        inline const double& GetImmed() const;
        inline unsigned GetVar() const;
        inline unsigned GetFuncNo() const;
        inline bool IsDefined() const { return GetOpcode() != FUNCTIONPARSERTYPES::cNop; }

        inline bool    IsImmed() const { return GetOpcode() == FUNCTIONPARSERTYPES::cImmed; }
        inline bool      IsVar() const { return GetOpcode() == FUNCTIONPARSERTYPES::VarBegin; }
        bool    IsLongIntegerImmed() const { return IsImmed() && GetImmed() == (double)GetLongIntegerImmed(); }
        long   GetLongIntegerImmed() const { return (long)GetImmed(); }
        bool    IsLogicalValue() const;
        inline unsigned GetRefCount() const;
        /* This function calculates the minimum and maximum values
         * of the tree's result. If an estimate cannot be made,
         * has_min/has_max are indicated as false.
         */
        MinMaxTree CalculateResultBoundaries_do() const;
        MinMaxTree CalculateResultBoundaries() const;

        enum TriTruthValue { IsAlways, IsNever, Unknown };
        TriTruthValue GetEvennessInfo() const;

        bool IsAlwaysSigned(bool positive) const;
        bool IsAlwaysParity(bool odd) const
            { return GetEvennessInfo() == (odd?IsNever:IsAlways); }
        bool IsAlwaysInteger(bool integer) const;

        void ConstantFolding();
        bool ConstantFolding_AndLogic();
        bool ConstantFolding_OrLogic();
        bool ConstantFolding_MulLogicItems();
        bool ConstantFolding_AddLogicItems();
        bool ConstantFolding_IfOperations();
        bool ConstantFolding_PowOperations();
        bool ConstantFolding_ComparisonOperations();
        template<typename CondType>
        bool ConstantFolding_LogicCommon(CondType cond_type, bool is_logical);
        bool ConstantFolding_MulGrouping();
        bool ConstantFolding_AddGrouping();
        bool ConstantFolding_Assimilate();

        void Rehash(bool constantfolding = true);
        inline void Mark_Incompletely_Hashed();
        inline bool Is_Incompletely_Hashed() const;

        inline const FPoptimizer_Grammar::Grammar* GetOptimizedUsing() const;
        inline void SetOptimizedUsing(const FPoptimizer_Grammar::Grammar* g);

        bool RecreateInversionsAndNegations(bool prefer_base2 = false);
        void FixIncompleteHashes();

        void swap(CodeTree& b) { data.swap(b.data); }
        bool IsIdenticalTo(const CodeTree& b) const;
        void CopyOnWrite();
    };

    struct CodeTreeData
    {
        int RefCount;

        /* Describing the codetree node */
        FUNCTIONPARSERTYPES::OPCODE Opcode;
        union
        {
            double   Value;   // In case of cImmed:   value of the immed
            unsigned Var;     // In case of VarBegin: variable number
            unsigned Funcno;  // In case of cFCall or cPCall
        };

        // Parameters for the function
        //  These use the sign:
        //   For cAdd: operands to add together (0 to n)
        //             sign indicates that the value is negated before adding (0-x)
        //   For cMul: operands to multiply together (0 to n)
        //             sign indicates that the value is inverted before multiplying (1/x)
        //   For cAnd: operands to bitwise-and together (0 to n)
        //             sign indicates that the value is inverted before anding (!x)
        //   For cOr:  operands to bitwise-or together (0 to n)
        //             sign indicates that the value is inverted before orring (!x)
        //  These don't use the sign (sign is always false):
        //   For cMin: operands to select the minimum of
        //   For cMax: operands to select the maximum of
        //   For cImmed, not used
        //   For VarBegin, not used
        //   For cIf:  operand 1 = condition, operand 2 = yes-branch, operand 3 = no-branch
        //   For anything else: the parameters required by the operation/function
        std::vector<CodeTree> Params;

        /* Internal operation */
        FUNCTIONPARSERTYPES::fphash_t      Hash;
        size_t        Depth;
        const FPoptimizer_Grammar::Grammar* OptimizedUsing;

        CodeTreeData();
        CodeTreeData(const CodeTreeData& b);
        explicit CodeTreeData(double i);
#ifdef __GXX_EXPERIMENTAL_CXX0X__
        CodeTreeData(CodeTreeData&& b);
#endif

        bool IsIdenticalTo(const CodeTreeData& b) const;
        void Sort();
        void Recalculate_Hash_NoRecursion();
    };

    inline void CodeTree::SetOpcode(FUNCTIONPARSERTYPES::OPCODE o)
        { data->Opcode = o; }
    inline void CodeTree::SetFuncOpcode(FUNCTIONPARSERTYPES::OPCODE o, unsigned f)
        { SetOpcode(o); data->Funcno = f; }
    inline void CodeTree::SetVar(unsigned v)
        { SetOpcode(FUNCTIONPARSERTYPES::VarBegin); data->Var = v; }
    inline void CodeTree::SetImmed(double v)
        { SetOpcode(FUNCTIONPARSERTYPES::cImmed); data->Value = v; }
    inline FUNCTIONPARSERTYPES::OPCODE CodeTree::GetOpcode() const { return data->Opcode; }
    inline FUNCTIONPARSERTYPES::fphash_t CodeTree::GetHash() const { return data->Hash; }
    inline const std::vector<CodeTree>& CodeTree::GetParams() const { return data->Params; }
    inline std::vector<CodeTree>& CodeTree::GetParams() { return data->Params; }
    inline size_t CodeTree::GetDepth() const { return data->Depth; }
    inline const double& CodeTree::GetImmed() const { return data->Value; }
    inline unsigned CodeTree::GetVar() const { return data->Var; }
    inline unsigned CodeTree::GetFuncNo() const { return data->Funcno; }

    inline const FPoptimizer_Grammar::Grammar* CodeTree::GetOptimizedUsing() const
        { return data->OptimizedUsing; }
    inline void CodeTree::SetOptimizedUsing(const FPoptimizer_Grammar::Grammar* g)
        { data->OptimizedUsing = g; }
    inline unsigned CodeTree::GetRefCount() const { return data->RefCount; }

    inline void CodeTree::Mark_Incompletely_Hashed() { data->Depth = 0; }
    inline bool CodeTree::Is_Incompletely_Hashed() const { return data->Depth == 0; }

#ifdef FUNCTIONPARSER_SUPPORT_DEBUG_OUTPUT
    void DumpHashes(const FPoptimizer_CodeTree::CodeTree& tree, std::ostream& o = std::cout);
    void DumpTree(const FPoptimizer_CodeTree::CodeTree& tree, std::ostream& o = std::cout);
    void DumpTreeWithIndent(const FPoptimizer_CodeTree::CodeTree& tree, std::ostream& o = std::cout, const std::string& indent = "\\");
#endif
}

#endif

#endif

#line 1 "fpoptimizer/fpoptimizer_grammar.h"
#ifndef FPOPT_NAN_CONST

#include <iostream>

#include "fpconfig.h"
#include "fparser.h"
#include "fptypes.h"

#define FPOPT_NAN_CONST (-1712345.25) /* Would use 0.0 / 0.0 here, but some compilers don't accept it. */

namespace FPoptimizer_CodeTree
{
    class CodeTree;
}

namespace FPoptimizer_Grammar
{
    enum ImmedConstraint_Value
    {
        ValueMask = 0x07,
        Value_AnyNum     = 0x0, // any value
        Value_EvenInt    = 0x1, // any even integer (0,2,4, etc)
        Value_OddInt     = 0x2, // any odd integer (1,3, etc)
        Value_IsInteger  = 0x3, // any integer-value (excludes e.g. 0.2)
        Value_NonInteger = 0x4, // any non-integer (excludes e.g. 1 or 5)
        Value_Logical    = 0x5  // a result of cNot,cNotNot,cAnd,cOr or comparators
    };
    enum ImmedConstraint_Sign
    {
        SignMask  = 0x18,
        Sign_AnySign     = 0x00, // - or +
        Sign_Positive    = 0x08, // positive only
        Sign_Negative    = 0x10, // negative only
        Sign_NoIdea      = 0x18  // where sign cannot be guessed
    };
    enum ImmedConstraint_Oneness
    {
        OnenessMask   = 0x60,
        Oneness_Any      = 0x00,
        Oneness_One      = 0x20, // +1 or -1
        Oneness_NotOne   = 0x40  // anything but +1 or -1
    };
    enum ImmedConstraint_Constness
    {
        ConstnessMask = 0x80,
        Constness_Any    = 0x00,
        Constness_Const  = 0x80
    };

    /* The param_opcode field of the ParamSpec has the following
     * possible values (from enum SpecialOpcode):
     *   NumConstant:
     *      this describes a specific constant value (constvalue)
     *      that must be matched / synthesized.
     *   ParamHolder:
     *      this describes any node
     *      that must be matched / synthesized.
     *      "index" is the ID of the NamedHolder:
     *      In matching, all NamedHolders having the same ID
     *      must match the identical node.
     *      In synthesizing, the node matched by
     *      a NamedHolder with this ID must be synthesized.
     *    SubFunction:
     *      this describes a subtree
     *      that must be matched / synthesized.
     *      The subtree is described in subfunc_opcode,param_begin..+param_count.
     *      If the type is GroupFunction, the tree is expected
     *      to yield a constant value which is tested.
     */
    enum SpecialOpcode
    {
        NumConstant,        // Holds a particular value (syntax-time constant)
        ParamHolder,        // Holds a particular named param
        SubFunction         // Holds an opcode and the params
    };

    enum ParamMatchingType
    {
        PositionalParams, // this set of params in this order
        SelectedParams,   // this set of params in any order
        AnyParams,        // these params are included
        GroupFunction     // this function represents a constant value
    };

    enum RuleType
    {
        ProduceNewTree, // replace self with the first (and only) from replaced_param
        ReplaceParams   // replace indicate params with replaced_params
    };

#ifdef __GNUC__
# define PACKED_GRAMMAR_ATTRIBUTE __attribute__((packed))
#else
# define PACKED_GRAMMAR_ATTRIBUTE
#endif

    enum { PARAM_INDEX_BITS = 9 };

    /* A ParamSpec object describes
     * either a parameter (leaf, node) that must be matched,
     * or a parameter (leaf, node) that must be synthesized.
     */
    typedef std::pair<SpecialOpcode, const void*> ParamSpec;
    ParamSpec ParamSpec_Extract(unsigned paramlist, unsigned index);
    bool ParamSpec_Compare(const void* a, const void* b, SpecialOpcode type);
    unsigned ParamSpec_GetDepCode(const ParamSpec& b);

    struct ParamSpec_ParamHolder
    {
        unsigned index       : 8; // holder ID
        unsigned constraints : 8; // constraints
        unsigned depcode     :16;
    } PACKED_GRAMMAR_ATTRIBUTE;

    struct ParamSpec_NumConstant
    {
        double constvalue;        // the value
    } PACKED_GRAMMAR_ATTRIBUTE;

    struct ParamSpec_SubFunctionData
    {
        /* Expected parameters (leaves) of the tree: */
        unsigned param_count         : 2;
        unsigned param_list          : 30;
        /* The opcode that the tree must have when SubFunction */
        FUNCTIONPARSERTYPES::OPCODE subfunc_opcode : 8;

        /* When matching, type describes the method of matching.
         *
         *               Sample input tree:      (cOr 2 3)  (cOr 2 4) (cOr 3 2) (cOr 4 2 3) (cOr 2)
         * Possible methods:
         *    PositionalParams, e.g. (cOr [2 3]):  match     no match  no match  no match   no match
         *      The nodes described here are
         *      to be matched, in this order.
         *    SelectedParams,   e.g. (cOr {2 3}):  match     no match   match    no match   no match
         *      The nodes described here are
         *      to be matched, in any order.
         *    AnyParams,        e.g. (cOr 2 3  ):  match     no match   match     match     no match
         *      At least the nodes described here
         *      are to be matched, in any order.
         * When synthesizing, the type is ignored.
         */
        ParamMatchingType match_type : 3; /* When SubFunction */
        /* Note: match_type needs 2, but we specify 3 because
         * otherwise Microsoft VC++ borks things up
         * as it interprets the value as signed.
         */
        /* Optional restholder index for capturing the rest of parameters (0=not used)
         * Only valid when match_type = AnyParams
         */
        unsigned restholder_index : 5;
    } PACKED_GRAMMAR_ATTRIBUTE; // size: 2+30+6+2+8=48 bits=6 bytes

    struct ParamSpec_SubFunction
    {
        ParamSpec_SubFunctionData data;
        unsigned constraints : 8; // constraints
        unsigned depcode     : 8;
    } PACKED_GRAMMAR_ATTRIBUTE; // 8 bytes

    /* Theoretical minimal sizes in each param_opcode cases:
     * Assume param_opcode needs 3 bits.
     *    NumConstant:   3 + 64              (or 3+4 if just use index to clist[])
     *    ParamHolder:   3 + 7 + 2           (7 for constraints, 2 for immed index)
     *    SubFunction:   3 + 7 + 2 + 2 + 3*9 = 41
     */

    /* A rule describes a pattern for matching
     * and the method how to reconstruct the
     * matched node(s) in the tree.
     */
    struct Rule
    {
        /* If the rule matched, this field describes how to perform
         * the replacement.
         *   When type==ProduceNewTree,
         *       the source tree is replaced entirely with
         *       the new tree described at repl_param_begin[0].
         *   When type==ReplaceParams,
         *       the matching leaves in the source tree are removed
         *       and new leaves are constructedfrom the trees
         *       described at repl_param_begin[0..repl_param_count].
         *       Other leaves remain intact.
         */
        RuleType  ruletype         : 2;
        bool      logical_context  : 1;

        /* The replacement parameters (if NewTree, begin[0] represents the new tree) */
        unsigned  repl_param_count : 2; /* Assumed to be 1 when type == ProduceNewTree */
        unsigned  repl_param_list  : 27;

        /* The function that we must match. Always a SubFunction. */
        ParamSpec_SubFunctionData match_tree;
    } PACKED_GRAMMAR_ATTRIBUTE; // size: 2+1+2+27 + 48 = 80 bits = 10 bytes

    /* Grammar is a set of rules for tree substitutions. */
    struct Grammar
    {
        /* The rules of this grammar */
        unsigned rule_count;
        unsigned char rule_list[999]; // maximum limit...
        /* Note: Changing the limit has no effect to performance of
         * fparser. The limit is only actually used within grammar_parser.
         * A too low limit causes a memory corruption during the parse.
         * A too high limit just may cause inconvenience.
         * The actual grammar items linked to fparser are optimized for size,
         * and the size of the Grammar object may be considerably smaller
         * than what is indicated by this prototype.
         */
    };

    extern "C" {
        extern const Rule      grammar_rules[];
    #ifndef FPOPTIMIZER_MERGED_FILE
        extern const Grammar   grammar_optimize_round1;
        extern const Grammar   grammar_optimize_round2;
        extern const Grammar   grammar_optimize_round3;
        extern const Grammar   grammar_optimize_round4;
        extern const Grammar   grammar_optimize_shortcut_logical_evaluation;
        extern const Grammar   grammar_optimize_nonshortcut_logical_evaluation;
        extern const Grammar   grammar_optimize_ignore_if_sideeffects;
        extern const Grammar   grammar_optimize_abslogical;
        extern const Grammar   grammar_optimize_base2_expand;
    #endif
    }

    void DumpParam(const ParamSpec& p, std::ostream& o = std::cout);
    void DumpParams(unsigned paramlist, unsigned count, std::ostream& o = std::cout);
}

#endif

#line 1 "fpoptimizer/fpoptimizer_consts.h"
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

#define CONSTANT_E     2.71828182845904523536  // exp(1)
#define CONSTANT_EI    0.3678794411714423215955 // exp(-1)
#define CONSTANT_2E    7.3890560989306502272304 // exp(2)
#define CONSTANT_2EI   0.135335283236612691894 // exp(-2)
#define CONSTANT_PI    M_PI                    // atan2(0,-1)
#define CONSTANT_L10   2.30258509299404590109  // log(10)
#define CONSTANT_L2    0.69314718055994530943  // log(2)
#define CONSTANT_L10I  0.43429448190325176116  // 1/log(10)
#define CONSTANT_L2I   1.4426950408889634074   // 1/log(2)
#define CONSTANT_L10E  CONSTANT_L10I           // log10(e)
#define CONSTANT_L10EI CONSTANT_L10            // 1/log10(e)
#define CONSTANT_L10B  0.3010299956639811952137 // log10(2)
#define CONSTANT_L10BI 3.3219280948873623478703 // 1/log10(2)
#define CONSTANT_LB10  CONSTANT_L10BI          // log2(10)
#define CONSTANT_LB10I CONSTANT_L10B           // 1/log2(10)
#define CONSTANT_L2E   CONSTANT_L2I            // log2(e)
#define CONSTANT_L2EI  CONSTANT_L2             // 1/log2(e)
#define CONSTANT_DR    (180.0 / M_PI)          // 180/pi
#define CONSTANT_RD    (M_PI / 180.0)          // pi/180

#define CONSTANT_POS_INF     HUGE_VAL  // positive infinity, from math.h
#define CONSTANT_NEG_INF   (-HUGE_VAL) // negative infinity
#define CONSTANT_PIHALF (M_PI / 2)

#line 1 "fpoptimizer/fpoptimizer_optimize.h"
// line removed for fpoptimizer.c: #include "fpoptimizer_codetree.h"
// line removed for fpoptimizer.c: #include "fpoptimizer_grammar.h"

#ifdef FP_SUPPORT_OPTIMIZER

#include <vector>
#include <utility>
#include <iostream>

//#define DEBUG_SUBSTITUTIONS

namespace FPoptimizer_Optimize
{
    using namespace FPoptimizer_Grammar;
    using namespace FPoptimizer_CodeTree;
    using namespace FUNCTIONPARSERTYPES;

    /* This struct collects information regarding the matching process so far */
    class MatchInfo
    {
    public:
        std::vector<std::pair<bool,std::vector<CodeTree>
                             > > restholder_matches;
        std::vector<CodeTree> paramholder_matches;
        std::vector<unsigned> matched_params;
    public:
        /* These functions save data from matching */
        bool SaveOrTestRestHolder(
            unsigned restholder_index,
            const std::vector<CodeTree>& treelist)
        {
            if(restholder_matches.size() <= restholder_index)
            {
                restholder_matches.resize(restholder_index+1);
                restholder_matches[restholder_index].first  = true;
                restholder_matches[restholder_index].second = treelist;
                return true;
            }
            if(restholder_matches[restholder_index].first == false)
            {
                restholder_matches[restholder_index].first  = true;
                restholder_matches[restholder_index].second = treelist;
                return true;
            }
            const std::vector<CodeTree>& found =
                restholder_matches[restholder_index].second;
            if(treelist.size() != found.size())
                return false;
            for(size_t a=0; a<treelist.size(); ++a)
                if(!treelist[a].IsIdenticalTo(found[a]))
                    return false;
            return true;
        }

        void SaveRestHolder(
            unsigned restholder_index,
            std::vector<CodeTree>& treelist)
        {
            if(restholder_matches.size() <= restholder_index)
                restholder_matches.resize(restholder_index+1);
            restholder_matches[restholder_index].first = true;
            restholder_matches[restholder_index].second.swap(treelist);
        }

        bool SaveOrTestParamHolder(
            unsigned paramholder_index,
            const CodeTree& treeptr)
        {
            if(paramholder_matches.size() <= paramholder_index)
            {
                paramholder_matches.reserve(paramholder_index+1);
                paramholder_matches.resize(paramholder_index);
                paramholder_matches.push_back(treeptr);
                return true;
            }
            if(!paramholder_matches[paramholder_index].IsDefined())
            {
                paramholder_matches[paramholder_index] = treeptr;
                return true;
            }
            return treeptr.IsIdenticalTo(paramholder_matches[paramholder_index]);
        }

        void SaveMatchedParamIndex(unsigned index)
        {
            matched_params.push_back(index);
        }

        /* These functions retrieve the data from matching
         * for use when synthesizing the resulting tree.
         */
        const CodeTree& GetParamHolderValueIfFound( unsigned paramholder_index ) const
        {
            static const CodeTree dummytree;
            if(paramholder_matches.size() <= paramholder_index)
                return dummytree;
            return paramholder_matches[paramholder_index];
        }

        const CodeTree& GetParamHolderValue( unsigned paramholder_index ) const
            { return paramholder_matches[paramholder_index]; }

        bool HasRestHolder(unsigned restholder_index) const
            { return restholder_matches.size() > restholder_index
                  && restholder_matches[restholder_index].first == true; }

        const std::vector<CodeTree>& GetRestHolderValues( unsigned restholder_index ) const
        {
            static const std::vector<CodeTree> empty_result;
            if(restholder_matches.size() <= restholder_index)
                return empty_result;
            return restholder_matches[restholder_index].second;
        }

        const std::vector<unsigned>& GetMatchedParamIndexes() const
            { return matched_params; }

        /* */
        void swap(MatchInfo& b)
        {
            restholder_matches.swap(b.restholder_matches);
            paramholder_matches.swap(b.paramholder_matches);
            matched_params.swap(b.matched_params);
        }
        MatchInfo& operator=(const MatchInfo& b)
        {
            restholder_matches = b.restholder_matches;
            paramholder_matches = b.paramholder_matches;
            matched_params = b.matched_params;
            return *this;
        }
    };

    class MatchPositionSpecBase;

    // For iterating through match candidates
    typedef FPOPT_autoptr<MatchPositionSpecBase> MatchPositionSpecBaseP;

    class MatchPositionSpecBase
    {
    public:
        int RefCount;
    public:
        MatchPositionSpecBase() : RefCount(0) { }
        virtual ~MatchPositionSpecBase() { }
    };
    struct MatchResultType
    {
        bool found;
        MatchPositionSpecBaseP specs;

        MatchResultType(bool f) : found(f), specs() { }
        MatchResultType(bool f,
                        const MatchPositionSpecBaseP& s) : found(f), specs(s) { }
    };

    /* Synthesize the given grammatic rule's replacement into the codetree */
    void SynthesizeRule(
        const Rule& rule,
        CodeTree& tree,
        MatchInfo& info);

    /* Test the given parameter to a given CodeTree */
    MatchResultType TestParam(
        const ParamSpec& parampair,
        const CodeTree& tree,
        const MatchPositionSpecBaseP& start_at,
        MatchInfo& info);

    /* Test the list of parameters to a given CodeTree */
    MatchResultType TestParams(
        const ParamSpec_SubFunctionData& model_tree,
        const CodeTree& tree,
        const MatchPositionSpecBaseP& start_at,
        MatchInfo& info,
        bool TopLevel);

    bool ApplyGrammar(const Grammar& grammar,
                      FPoptimizer_CodeTree::CodeTree& tree,
                      bool from_logical_context = false);
    void ApplyGrammars(FPoptimizer_CodeTree::CodeTree& tree);

    bool IsLogisticallyPlausibleParamsMatch(
        const ParamSpec_SubFunctionData& params,
        const CodeTree& tree);
}

namespace FPoptimizer_Grammar
{
    void DumpMatch(const Rule& rule,
                   const FPoptimizer_CodeTree::CodeTree& tree,
                   const FPoptimizer_Optimize::MatchInfo& info,
                   bool DidMatch,
                   std::ostream& o = std::cout);
    void DumpMatch(const Rule& rule,
                   const FPoptimizer_CodeTree::CodeTree& tree,
                   const FPoptimizer_Optimize::MatchInfo& info,
                   const char* whydump,
                   std::ostream& o = std::cout);
}

#endif

#line 1 "fpoptimizer/crc32.h"
/* crc32 */

#ifdef _MSC_VER

 typedef unsigned int crc32_t;

#else

 #include <stdint.h>
 typedef uint_least32_t crc32_t;

#endif

namespace crc32
{
    enum { startvalue = 0xFFFFFFFFUL, poly = 0xEDB88320UL };

    /* This code constructs the CRC32 table at compile-time,
     * avoiding the need for a huge explicitly written table of magical numbers. */
    template<crc32_t crc> // One byte of a CRC32 (eight bits):
    struct b8
    {
        enum { b1 = (crc & 1) ? (poly ^ (crc >> 1)) : (crc >> 1),
               b2 = (b1  & 1) ? (poly ^ (b1  >> 1)) : (b1  >> 1),
               b3 = (b2  & 1) ? (poly ^ (b2  >> 1)) : (b2  >> 1),
               b4 = (b3  & 1) ? (poly ^ (b3  >> 1)) : (b3  >> 1),
               b5 = (b4  & 1) ? (poly ^ (b4  >> 1)) : (b4  >> 1),
               b6 = (b5  & 1) ? (poly ^ (b5  >> 1)) : (b5  >> 1),
               b7 = (b6  & 1) ? (poly ^ (b6  >> 1)) : (b6  >> 1),
               res= (b7  & 1) ? (poly ^ (b7  >> 1)) : (b7  >> 1) };
    };
    inline crc32_t update(crc32_t crc, unsigned/* char */b) // __attribute__((pure))
    {
        // Four values of the table
        #define B4(n) b8<n>::res,b8<n+1>::res,b8<n+2>::res,b8<n+3>::res
        // Sixteen values of the table
        #define R(n) B4(n),B4(n+4),B4(n+8),B4(n+12)
        // The whole table, index by steps of 16
        static const crc32_t table[256] =
        { R(0x00),R(0x10),R(0x20),R(0x30), R(0x40),R(0x50),R(0x60),R(0x70),
          R(0x80),R(0x90),R(0xA0),R(0xB0), R(0xC0),R(0xD0),R(0xE0),R(0xF0) };
        #undef R
        #undef B4
        return ((crc >> 8) /* & 0x00FFFFFF*/) ^ table[/*(unsigned char)*/(crc^b)&0xFF];
    }
    inline crc32_t calc_upd(crc32_t c, const unsigned char* buf, size_t size)
    {
        crc32_t value = c;
        for(size_t p=0; p<size; ++p) value = update(value, buf[p]);
        return value;
    }
    inline crc32_t calc(const unsigned char* buf, size_t size)
    {
        return calc_upd(startvalue, buf, size);
    }
}

#line 1 "fpoptimizer/fpoptimizer_opcodename.h"
// line removed for fpoptimizer.c: #include "fpoptimizer_grammar.h"
#include <string>

const std::string FP_GetOpcodeName(FPoptimizer_Grammar::SpecialOpcode opcode, bool pad=false);
const std::string FP_GetOpcodeName(FUNCTIONPARSERTYPES::OPCODE opcode,        bool pad=false);

#line 1 "fpoptimizer/fpoptimizer_opcodename.c"
#include <string>
#include <sstream>
#include <assert.h>

#include <iostream>

#include "fpconfig.h"
#include "fptypes.h"

// line removed for fpoptimizer.c: #include "fpoptimizer_grammar.h"
// line removed for fpoptimizer.c: #include "fpoptimizer_opcodename.h"

using namespace FPoptimizer_Grammar;
using namespace FUNCTIONPARSERTYPES;

const std::string FP_GetOpcodeName(FPoptimizer_Grammar::SpecialOpcode opcode, bool pad)
{
#if 1
    /* Symbolic meanings for the opcodes? */
    const char* p = 0;
    switch( opcode )
    {
        case NumConstant:   p = "NumConstant"; break;
        case ParamHolder:   p = "ParamHolder"; break;
        case SubFunction:   p = "SubFunction"; break;
    }
    std::ostringstream tmp;
    //if(!p) std::cerr << "o=" << opcode << "\n";
    assert(p);
    tmp << p;
    if(pad) while(tmp.str().size() < 12) tmp << ' ';
    return tmp.str();
#else
    /* Just numeric meanings */
    std::ostringstream tmp;
    tmp << opcode;
    if(pad) while(tmp.str().size() < 5) tmp << ' ';
    return tmp.str();
#endif
}

const std::string FP_GetOpcodeName(FUNCTIONPARSERTYPES::OPCODE opcode,        bool pad)
{
#if 1
    /* Symbolic meanings for the opcodes? */
    const char* p = 0;
    switch(opcode)
    {
        case cAbs: p = "cAbs"; break;
        case cAcos: p = "cAcos"; break;
        case cAcosh: p = "cAcosh"; break;
        case cAsin: p = "cAsin"; break;
        case cAsinh: p = "cAsinh"; break;
        case cAtan: p = "cAtan"; break;
        case cAtan2: p = "cAtan2"; break;
        case cAtanh: p = "cAtanh"; break;
        case cCbrt: p = "cCbrt"; break;
        case cCeil: p = "cCeil"; break;
        case cCos: p = "cCos"; break;
        case cCosh: p = "cCosh"; break;
        case cCot: p = "cCot"; break;
        case cCsc: p = "cCsc"; break;
        case cEval: p = "cEval"; break;
        case cExp: p = "cExp"; break;
        case cExp2: p = "cExp2"; break;
        case cFloor: p = "cFloor"; break;
        case cIf: p = "cIf"; break;
        case cInt: p = "cInt"; break;
        case cLog: p = "cLog"; break;
        case cLog2: p = "cLog2"; break;
        case cLog10: p = "cLog10"; break;
        case cMax: p = "cMax"; break;
        case cMin: p = "cMin"; break;
        case cPow: p = "cPow"; break;
        case cSec: p = "cSec"; break;
        case cSin: p = "cSin"; break;
        case cSinh: p = "cSinh"; break;
        case cSqrt: p = "cSqrt"; break;
        case cTan: p = "cTan"; break;
        case cTanh: p = "cTanh"; break;
        case cTrunc: p = "cTrunc"; break;
        case cImmed: p = "cImmed"; break;
        case cJump: p = "cJump"; break;
        case cNeg: p = "cNeg"; break;
        case cAdd: p = "cAdd"; break;
        case cSub: p = "cSub"; break;
        case cMul: p = "cMul"; break;
        case cDiv: p = "cDiv"; break;
        case cMod: p = "cMod"; break;
        case cEqual: p = "cEqual"; break;
        case cNEqual: p = "cNEqual"; break;
        case cLess: p = "cLess"; break;
        case cLessOrEq: p = "cLessOrEq"; break;
        case cGreater: p = "cGreater"; break;
        case cGreaterOrEq: p = "cGreaterOrEq"; break;
        case cNot: p = "cNot"; break;
        case cAnd: p = "cAnd"; break;
        case cOr: p = "cOr"; break;
        case cDeg: p = "cDeg"; break;
        case cRad: p = "cRad"; break;
        case cFCall: p = "cFCall"; break;
        case cPCall: p = "cPCall"; break;
#ifdef FP_SUPPORT_OPTIMIZER
        case cFetch: p = "cFetch"; break;
        case cPopNMov: p = "cPopNMov"; break;
        case cLog2by: p = "cLog2by"; break;
#endif
        case cAbsNot: p = "cAbsNot"; break;
        case cAbsNotNot: p = "cAbsNotNot"; break;
        case cAbsAnd: p = "cAbsAnd"; break;
        case cAbsOr: p = "cAbsOr"; break;
        case cAbsIf: p = "cAbsIf"; break;
        case cDup: p = "cDup"; break;
        case cInv: p = "cInv"; break;
        case cSqr: p = "cSqr"; break;
        case cRDiv: p = "cRDiv"; break;
        case cRSub: p = "cRSub"; break;
        case cNotNot: p = "cNotNot"; break;
        case cRSqrt: p = "cRSqrt"; break;
        case cNop: p = "cNop"; break;
        case VarBegin: p = "VarBegin"; break;
    }
    std::ostringstream tmp;
    //if(!p) std::cerr << "o=" << opcode << "\n";
    assert(p);
    tmp << p;
    if(pad) while(tmp.str().size() < 12) tmp << ' ';
    return tmp.str();
#else
    /* Just numeric meanings */
    std::ostringstream tmp;
    tmp << opcode;
    if(pad) while(tmp.str().size() < 5) tmp << ' ';
    return tmp.str();
#endif
}

#line 1 "fpoptimizer/fpoptimizer_bytecodesynth.h"
#include "fpconfig.h"
#include "fparser.h"
#include "fptypes.h"

#ifdef FP_SUPPORT_OPTIMIZER

#include <vector>
#include <utility>

// line removed for fpoptimizer.c: #include "fpoptimizer_codetree.h"

#ifndef FP_GENERATING_POWI_TABLE
enum { MAX_POWI_BYTECODE_LENGTH = 20 };
#else
enum { MAX_POWI_BYTECODE_LENGTH = 999 };
#endif
enum { MAX_MULI_BYTECODE_LENGTH = 3 };

namespace FPoptimizer_ByteCode
{
    class ByteCodeSynth
    {
    public:
        ByteCodeSynth()
            : ByteCode(), Immed(), StackTop(0), StackMax(0)
        {
            /* estimate the initial requirements as such */
            ByteCode.reserve(64);
            Immed.reserve(8);
            StackState.reserve(16);
        }

        void Pull(std::vector<unsigned>& bc,
                  std::vector<double>&   imm,
                  size_t& StackTop_max)
        {
            ByteCode.swap(bc);
            Immed.swap(imm);
            StackTop_max = StackMax;
        }

        size_t GetByteCodeSize() const { return ByteCode.size(); }
        size_t GetStackTop()     const { return StackTop; }

        void PushVar(unsigned varno)
        {
            ByteCode.push_back(varno);
            SetStackTop(StackTop+1);
        }

        void PushImmed(double immed)
        {
            using namespace FUNCTIONPARSERTYPES;
            ByteCode.push_back(cImmed);
            Immed.push_back(immed);
            SetStackTop(StackTop+1);
        }

        void StackTopIs(const FPoptimizer_CodeTree::CodeTree& tree)
        {
            if(StackTop > 0)
            {
                StackState[StackTop-1].first = true;
                StackState[StackTop-1].second = tree;
            }
        }

        void EatNParams(unsigned eat_count)
        {
            SetStackTop(StackTop - eat_count);
        }

        void ProducedNParams(unsigned produce_count)
        {
            SetStackTop(StackTop + produce_count);
        }

        void AddOperation(unsigned opcode, unsigned eat_count, unsigned produce_count = 1)
        {
            EatNParams(eat_count);

            using namespace FUNCTIONPARSERTYPES;

            if(!ByteCode.empty() && opcode == cMul && ByteCode.back() == cDup)
                ByteCode.back() = cSqr;
            else
                ByteCode.push_back(opcode);

            ProducedNParams(produce_count);
        }

        void DoPopNMov(size_t targetpos, size_t srcpos)
        {
            using namespace FUNCTIONPARSERTYPES;
            ByteCode.push_back(cPopNMov);
            ByteCode.push_back( (unsigned) targetpos);
            ByteCode.push_back( (unsigned) srcpos);

            SetStackTop(srcpos+1);
            StackState[targetpos] = StackState[srcpos];
            SetStackTop(targetpos+1);
        }

        void DoDup(size_t src_pos)
        {
            using namespace FUNCTIONPARSERTYPES;
            if(src_pos == StackTop-1)
            {
                ByteCode.push_back(cDup);
            }
            else
            {
                ByteCode.push_back(cFetch);
                ByteCode.push_back( (unsigned) src_pos);
            }
            SetStackTop(StackTop + 1);
            StackState[StackTop-1] = StackState[src_pos];
        }

        size_t FindPos(const FPoptimizer_CodeTree::CodeTree& tree) const
        {
            /*std::cout << "Stack state now(" << StackTop << "):\n";
            for(size_t a=0; a<StackTop; ++a)
            {
                std::cout << a << ": ";
                if(StackState[a].first)
                    DumpTree(StackState[a].second);
                else
                    std::cout << "?";
                std::cout << "\n";
            }*/
            for(size_t a=StackTop; a-->0; )
                if(StackState[a].first && StackState[a].second.IsIdenticalTo(tree))
                    return a;
            return ~size_t(0);
        }

        bool Find(const FPoptimizer_CodeTree::CodeTree& tree) const
        {
            return FindPos(tree) != ~size_t(0);
        }

        bool FindAndDup(const FPoptimizer_CodeTree::CodeTree& tree)
        {
            size_t pos = FindPos(tree);
            if(pos != ~size_t(0))
            {
                DoDup(pos);
                return true;
            }
            return false;
        }

        struct IfData
        {
            size_t ofs;
        };

        void SynthIfStep1(IfData& ifdata, FUNCTIONPARSERTYPES::OPCODE op)
        {
            using namespace FUNCTIONPARSERTYPES;
            SetStackTop(StackTop-1); // the If condition was popped.

            ifdata.ofs = ByteCode.size();
            ByteCode.push_back(op);
            ByteCode.push_back(0); // code index
            ByteCode.push_back(0); // Immed index
        }
        void SynthIfStep2(IfData& ifdata)
        {
            using namespace FUNCTIONPARSERTYPES;
            SetStackTop(StackTop-1); // ignore the pushed then-branch result.

            ByteCode[ifdata.ofs+1] = unsigned( ByteCode.size()+2 );
            ByteCode[ifdata.ofs+2] = unsigned( Immed.size()      );

            ifdata.ofs = ByteCode.size();
            ByteCode.push_back(cJump);
            ByteCode.push_back(0); // code index
            ByteCode.push_back(0); // Immed index
        }
        void SynthIfStep3(IfData& ifdata)
        {
            using namespace FUNCTIONPARSERTYPES;
            SetStackTop(StackTop-1); // ignore the pushed else-branch result.

            ByteCode[ifdata.ofs+1] = unsigned( ByteCode.size()-1 );
            ByteCode[ifdata.ofs+2] = unsigned( Immed.size()      );

            SetStackTop(StackTop+1); // one or the other was pushed.

            /* Threading jumps:
             * If there are any cJumps that point
             * to the cJump instruction we just changed,
             * change them to point to this target as well.
             * This screws up PrintByteCode() majorly.
             */
            for(size_t a=0; a<ifdata.ofs; ++a)
            {
                if(ByteCode[a]   == cJump
                && ByteCode[a+1] == ifdata.ofs-1)
                {
                    ByteCode[a+1] = unsigned( ByteCode.size()-1 );
                    ByteCode[a+2] = unsigned( Immed.size()      );
                }
                switch(ByteCode[a])
                {
                    case cAbsIf:
                    case cIf:
                    case cJump:
                    case cPopNMov: a += 2; break;
                    case cFCall:
                    case cPCall:
                    case cFetch: a += 1; break;
                    default: break;
                }
            }
        }

    protected:
        void SetStackTop(size_t value)
        {
            StackTop = value;
            if(StackTop > StackMax)
            {
                StackMax = StackTop;
                StackState.resize(StackMax);
            }
        }

    private:
        std::vector<unsigned> ByteCode;
        std::vector<double>   Immed;

        std::vector<
            std::pair<bool/*known*/, FPoptimizer_CodeTree::CodeTree/*tree*/>
                   > StackState;
        size_t StackTop;
        size_t StackMax;
    };

    struct SequenceOpCode;
    extern const SequenceOpCode AddSequence; /* Multiplication implemented with adds */
    extern const SequenceOpCode MulSequence; /* Exponentiation implemented with muls */

    /* Generate a sequence that multiplies or exponentifies the
     * last operand in the stack by the given constant integer
     * amount (positive or negative).
     */
    void AssembleSequence(
        long count,
        const SequenceOpCode& sequencing,
        ByteCodeSynth& synth);
}

#endif

#line 1 "fpoptimizer/fpoptimizer_bytecodesynth.c"
// line removed for fpoptimizer.c: #include "fpoptimizer_bytecodesynth.h"

#ifdef FP_SUPPORT_OPTIMIZER

using namespace FUNCTIONPARSERTYPES;

namespace FPoptimizer_ByteCode
{
    const struct SequenceOpCode
    {
        double basevalue;
        unsigned op_flip;
        unsigned op_normal, op_normal_flip;
        unsigned op_inverse, op_inverse_flip;
    } AddSequence = {0.0, cNeg, cAdd, cAdd, cSub, cRSub },
      MulSequence = {1.0, cInv, cMul, cMul, cDiv, cRDiv };
}

using namespace FPoptimizer_ByteCode;

#define POWI_TABLE_SIZE 256
#define POWI_WINDOW_SIZE 3
namespace FPoptimizer_ByteCode
{
    #ifndef FP_GENERATING_POWI_TABLE
    extern const
    unsigned char powi_table[POWI_TABLE_SIZE];
    const
    #endif
    unsigned char powi_table[POWI_TABLE_SIZE] =
    {
          0,   1,   1,   1,   2,   1,   2,   1, /*   0 -   7 */
          4,   1,   2,   1,   4,   1,   2, 131, /*   8 -  15 */
          8,   1,   2,   1,   4,   1,   2,   1, /*  16 -  23 */
          8, 133,   2, 131,   4,   1,  15,   1, /*  24 -  31 */
         16,   1,   2,   1,   4,   1,   2, 131, /*  32 -  39 */
          8,   1,   2,   1,   4, 133,   2,   1, /*  40 -  47 */
         16,   1,  25, 131,   4,   1,  27,   5, /*  48 -  55 */
          8,   3,   2,   1,  30,   1,  31,   3, /*  56 -  63 */
         32,   1,   2,   1,   4,   1,   2,   1, /*  64 -  71 */
          8,   1,   2, 131,   4,   1,  39,   1, /*  72 -  79 */
         16, 137,   2,   1,   4, 133,   2, 131, /*  80 -  87 */
          8,   1,  45, 135,   4,  31,   2,   5, /*  88 -  95 */
         32,   1,   2, 131,  50,   1,  51,   1, /*  96 - 103 */
          8,   3,   2,   1,  54,   1,  55,   3, /* 104 - 111 */
         16,   1,  57, 133,   4, 137,   2, 135, /* 112 - 119 */
         60,   1,  61,   3,  62, 133,  63,   1, /* 120 - 127 */
        130,   1,   2,   1, 130,   1,   2, 131, /* 128 - 135 */
        130,   1,   2,   1, 130,   1,   2, 139, /* 136 - 143 */
        130,   1,   2, 131, 130,   1,  30,   1, /* 144 - 151 */
        130, 137,   2,  31, 130,   1,   2, 131, /* 152 - 159 */
        130,   1, 130,   1, 130, 133,   2,   1, /* 160 - 167 */
        130,   1, 130,   1,   2,   1, 130, 133, /* 168 - 175 */
        130,   1,   2,   1, 130,   1,   2,  61, /* 176 - 183 */
        130, 133,  62, 139, 130, 137, 130,   1, /* 184 - 191 */
        130,   1,   2, 131, 130,   1, 130,   1, /* 192 - 199 */
        130,   1,   2,   1, 130,   1,   2, 131, /* 200 - 207 */
        130,   1, 130,   1, 130, 131,   2, 133, /* 208 - 215 */
        130,   1,   2, 131, 130, 141, 130,   1, /* 216 - 223 */
        130, 133,   2,   1, 130,   1,   5, 135, /* 224 - 231 */
        130,   1, 130,   1,   2, 131, 130,   1, /* 232 - 239 */
        130,   1,   2, 131, 130, 133, 130, 141, /* 240 - 247 */
        130, 131, 130,   1, 130,   1,   2, 131  /* 248 - 255 */
    }; /* as in gcc, but custom-optimized for stack calculation */
}
static const int POWI_CACHE_SIZE = 256;

#define FPO(x) /**/
//#define FPO(x) x


namespace
{
    class PowiCache
    {
    private:
        int cache[POWI_CACHE_SIZE];
        int cache_needed[POWI_CACHE_SIZE];

    public:
        PowiCache()
            : cache(), cache_needed() /* Assume we have no factors in the cache */
        {
            /* Decide which factors we would need multiple times.
             * Output:
             *   cache[]        = these factors were generated
             *   cache_needed[] = number of times these factors were desired
             */
            cache[1] = 1; // We have this value already.
        }

        bool Plan_Add(long value, int count)
        {
            if(value >= POWI_CACHE_SIZE) return false;
            //FPO(fprintf(stderr, "%ld will be needed %d times more\n", count, need_count));
            cache_needed[value] += count;
            return cache[value] != 0;
        }

        void Plan_Has(long value)
        {
            if(value < POWI_CACHE_SIZE)
                cache[value] = 1; // This value has been generated
        }

        void Start(size_t value1_pos)
        {
            for(int n=2; n<POWI_CACHE_SIZE; ++n)
                cache[n] = -1; /* Stack location for each component */

            Remember(1, value1_pos);

            DumpContents();
        }

        int Find(long value) const
        {
            if(value < POWI_CACHE_SIZE)
            {
                if(cache[value] >= 0)
                {
                    // found from the cache
                    FPO(fprintf(stderr, "* I found %ld from cache (%u,%d)\n",
                        value, (unsigned)cache[value], cache_needed[value]));
                    return cache[value];
                }
            }
            return -1;
        }

        void Remember(long value, size_t stackpos)
        {
            if(value >= POWI_CACHE_SIZE) return;

            FPO(fprintf(stderr, "* Remembering that %ld can be found at %u (%d uses remain)\n",
                value, (unsigned)stackpos, cache_needed[value]));
            cache[value] = (int) stackpos;
        }

        void DumpContents() const
        {
            FPO(for(int a=1; a<POWI_CACHE_SIZE; ++a)
                if(cache[a] >= 0 || cache_needed[a] > 0)
                {
                    fprintf(stderr, "== cache: sp=%d, val=%d, needs=%d\n",
                        cache[a], a, cache_needed[a]);
                })
        }

        int UseGetNeeded(long value)
        {
            if(value >= 0 && value < POWI_CACHE_SIZE)
                return --cache_needed[value];
            return 0;
        }
    };


    size_t AssembleSequence_Subdivide(
        long count,
        PowiCache& cache,
        const SequenceOpCode& sequencing,
        ByteCodeSynth& synth);

    void Subdivide_Combine(
        size_t apos, long aval,
        size_t bpos, long bval,
        PowiCache& cache,

        unsigned cumulation_opcode,
        unsigned cimulation_opcode_flip,

        ByteCodeSynth& synth);

    void PlanNtimesCache
        (long value,
         PowiCache& cache,
         int need_count,
         int recursioncount=0)
    {
        if(value < 1) return;

    #ifdef FP_GENERATING_POWI_TABLE
        if(recursioncount > 32) throw false;
    #endif

        if(cache.Plan_Add(value, need_count)) return;

        long half = 1;
        if(value < POWI_TABLE_SIZE)
        {
            half = powi_table[value];
            if(half & 128)
            {
                half &= 127;
                if(half & 64)
                    half = -(half & 63) - 1;

                FPO(fprintf(stderr, "value=%ld, half=%ld, otherhalf=%ld\n", value,half,value/half));

                PlanNtimesCache(half,      cache, 1, recursioncount+1);
                cache.Plan_Has(half);
                return;
            }
            else if(half & 64)
            {
                half = -(half & 63) - 1;
            }
        }
        else if(value & 1)
            half = value & ((1 << POWI_WINDOW_SIZE) - 1); // that is, value & 7
        else
            half = value / 2;

        long otherhalf = value-half;
        if(half > otherhalf || half<0) std::swap(half,otherhalf);

        FPO(fprintf(stderr, "value=%ld, half=%ld, otherhalf=%ld\n", value,half,otherhalf));

        if(half == otherhalf)
        {
            PlanNtimesCache(half,      cache, 2, recursioncount+1);
        }
        else
        {
            PlanNtimesCache(half,      cache, 1, recursioncount+1);
            PlanNtimesCache(otherhalf>0?otherhalf:-otherhalf,
                                       cache, 1, recursioncount+1);
        }
        cache.Plan_Has(value);
    }

    size_t AssembleSequence_Subdivide(
        long value,
        PowiCache& cache,
        const SequenceOpCode& sequencing,
        ByteCodeSynth& synth)
    {
        int cachepos = cache.Find(value);
        if(cachepos >= 0)
        {
            // found from the cache
            return cachepos;
        }

        long half = 1;
        if(value < POWI_TABLE_SIZE)
        {
            half = powi_table[value];
            if(half & 128)
            {
                half &= 127;
                if(half & 64)
                    half = -(half & 63) - 1;

                FPO(fprintf(stderr, "* I want %ld, my plan is %ld * %ld\n", value, half, value/half));
                size_t half_pos = AssembleSequence_Subdivide(half, cache, sequencing, synth);
                if(cache.UseGetNeeded(half) > 0
                || half_pos != synth.GetStackTop()-1)
                {
                    synth.DoDup(half_pos);
                    cache.Remember(half, synth.GetStackTop()-1);
                }
                AssembleSequence(value/half, sequencing, synth);
                size_t stackpos = synth.GetStackTop()-1;
                cache.Remember(value, stackpos);
                cache.DumpContents();
                return stackpos;
            }
            else if(half & 64)
            {
                half = -(half & 63) - 1;
            }
        }
        else if(value & 1)
            half = value & ((1 << POWI_WINDOW_SIZE) - 1); // that is, value & 7
        else
            half = value / 2;

        long otherhalf = value-half;
        if(half > otherhalf || half<0) std::swap(half,otherhalf);

        FPO(fprintf(stderr, "* I want %ld, my plan is %ld + %ld\n", value, half, value-half));

        if(half == otherhalf)
        {
            size_t half_pos = AssembleSequence_Subdivide(half, cache, sequencing, synth);

            // self-cumulate the subdivide result
            Subdivide_Combine(half_pos,half, half_pos,half, cache,
                sequencing.op_normal, sequencing.op_normal_flip,
                synth);
        }
        else
        {
            long part1 = half;
            long part2 = otherhalf>0?otherhalf:-otherhalf;

            size_t part1_pos = AssembleSequence_Subdivide(part1, cache, sequencing, synth);
            size_t part2_pos = AssembleSequence_Subdivide(part2, cache, sequencing, synth);

            FPO(fprintf(stderr, "Subdivide(%ld: %ld, %ld)\n", value, half, otherhalf));

            Subdivide_Combine(part1_pos,part1, part2_pos,part2, cache,
                otherhalf>0 ? sequencing.op_normal      : sequencing.op_inverse,
                otherhalf>0 ? sequencing.op_normal_flip : sequencing.op_inverse_flip,
                synth);
        }

        size_t stackpos = synth.GetStackTop()-1;
        cache.Remember(value, stackpos);
        cache.DumpContents();
        return stackpos;
    }

    void Subdivide_Combine(
        size_t apos, long aval,
        size_t bpos, long bval,
        PowiCache& cache,
        unsigned cumulation_opcode,
        unsigned cumulation_opcode_flip,
        ByteCodeSynth& synth)
    {
        /*FPO(fprintf(stderr, "== making result for (sp=%u, val=%d, needs=%d) and (sp=%u, val=%d, needs=%d), stacktop=%u\n",
            (unsigned)apos, aval, aval>=0 ? cache_needed[aval] : -1,
            (unsigned)bpos, bval, bval>=0 ? cache_needed[bval] : -1,
            (unsigned)synth.GetStackTop()));*/

        // Figure out whether we can trample a and b
        int a_needed = cache.UseGetNeeded(aval);
        int b_needed = cache.UseGetNeeded(bval);

        bool flipped = false;

        #define DUP_BOTH() do { \
            if(apos < bpos) { size_t tmp=apos; apos=bpos; bpos=tmp; flipped=!flipped; } \
            FPO(fprintf(stderr, "-> dup(%u) dup(%u) op\n", (unsigned)apos, (unsigned)bpos)); \
            synth.DoDup(apos); \
            synth.DoDup(apos==bpos ? synth.GetStackTop()-1 : bpos); } while(0)
        #define DUP_ONE(p) do { \
            FPO(fprintf(stderr, "-> dup(%u) op\n", (unsigned)p)); \
            synth.DoDup(p); \
        } while(0)

        if(a_needed > 0)
        {
            if(b_needed > 0)
            {
                // If they must both be preserved, make duplicates
                // First push the one that is at the larger stack
                // address. This increases the odds of possibly using cDup.
                DUP_BOTH();

                //SCENARIO 1:
                // Input:  x B A x x
                // Temp:   x B A x x A B
                // Output: x B A x x R
                //SCENARIO 2:
                // Input:  x A B x x
                // Temp:   x A B x x B A
                // Output: x A B x x R
            }
            else
            {
                // A must be preserved, but B can be trampled over

                // SCENARIO 1:
                //  Input:  x B x x A
                //   Temp:  x B x x A A B   (dup both, later first)
                //  Output: x B x x A R
                // SCENARIO 2:
                //  Input:  x A x x B
                //   Temp:  x A x x B A
                //  Output: x A x x R       -- only commutative cases
                // SCENARIO 3:
                //  Input:  x x x B A
                //   Temp:  x x x B A A B   (dup both, later first)
                //  Output: x x x B A R
                // SCENARIO 4:
                //  Input:  x x x A B
                //   Temp:  x x x A B A     -- only commutative cases
                //  Output: x x x A R
                // SCENARIO 5:
                //  Input:  x A B x x
                //   Temp:  x A B x x A B   (dup both, later first)
                //  Output: x A B x x R

                // if B is not at the top, dup both.
                if(bpos != synth.GetStackTop()-1)
                    DUP_BOTH();    // dup both
                else
                {
                    DUP_ONE(apos); // just dup A
                    flipped=!flipped;
                }
            }
        }
        else if(b_needed > 0)
        {
            // B must be preserved, but A can be trampled over
            // This is a mirror image of the a_needed>0 case, so I'll cut the chase
            if(apos != synth.GetStackTop()-1)
                DUP_BOTH();
            else
                DUP_ONE(bpos);
        }
        else
        {
            // Both can be trampled over.
            // SCENARIO 1:
            //  Input:  x B x x A
            //   Temp:  x B x x A B
            //  Output: x B x x R
            // SCENARIO 2:
            //  Input:  x A x x B
            //   Temp:  x A x x B A
            //  Output: x A x x R       -- only commutative cases
            // SCENARIO 3:
            //  Input:  x x x B A
            //  Output: x x x R         -- only commutative cases
            // SCENARIO 4:
            //  Input:  x x x A B
            //  Output: x x x R
            // SCENARIO 5:
            //  Input:  x A B x x
            //   Temp:  x A B x x A B   (dup both, later first)
            //  Output: x A B x x R
            // SCENARIO 6:
            //  Input:  x x x C
            //   Temp:  x x x C C   (c is both A and B)
            //  Output: x x x R

            if(apos == bpos && apos == synth.GetStackTop()-1)
                DUP_ONE(apos); // scenario 6
            else if(apos == synth.GetStackTop()-1 && bpos == synth.GetStackTop()-2)
            {
                FPO(fprintf(stderr, "-> op\n")); // scenario 3
                flipped=!flipped;
            }
            else if(apos == synth.GetStackTop()-2 && bpos == synth.GetStackTop()-1)
                FPO(fprintf(stderr, "-> op\n")); // scenario 4
            else if(apos == synth.GetStackTop()-1)
                DUP_ONE(bpos); // scenario 1
            else if(bpos == synth.GetStackTop()-1)
            {
                DUP_ONE(apos); // scenario 2
                flipped=!flipped;
            }
            else
                DUP_BOTH(); // scenario 5
        }
        // Add them together.
        synth.AddOperation(flipped ? cumulation_opcode_flip : cumulation_opcode, 2);
    }

    void LightWeight(
        long count,
        const SequenceOpCode& sequencing,
        ByteCodeSynth& synth)
    {
        while(count < 256)
        {
            int half = FPoptimizer_ByteCode::powi_table[count];
            if(half & 128)
            {
                half &= 127;
                LightWeight(half,       sequencing, synth);
                count /= half;
            }
            else break;
        }
        if(count == 1) return;
        if(!(count & 1))
        {
            synth.AddOperation(cSqr, 1);
            LightWeight(count/2, sequencing, synth);
        }
        else
        {
            synth.DoDup(synth.GetStackTop()-1);
            LightWeight(count-1, sequencing, synth);
            synth.AddOperation(cMul, 2);
        }
    }
}

namespace FPoptimizer_ByteCode
{
    void AssembleSequence(
        long count,
        const SequenceOpCode& sequencing,
        ByteCodeSynth& synth)
    {
        if(count == 0)
            synth.PushImmed(sequencing.basevalue);
        else
        {
            bool needs_flip = false;
            if(count < 0)
            {
                needs_flip = true;
                count = -count;
            }

            if(false)
                LightWeight(count,sequencing,synth);
            else if(count > 1)
            {
                /* To prevent calculating the same factors over and over again,
                 * we use a cache. */
                PowiCache cache;
                PlanNtimesCache(count, cache, 1);

                size_t stacktop_desired = synth.GetStackTop();

                cache.Start( synth.GetStackTop()-1 );

                FPO(fprintf(stderr, "Calculating result for %ld...\n", count));
                size_t res_stackpos = AssembleSequence_Subdivide(
                    count, cache, sequencing,
                    synth);

                size_t n_excess = synth.GetStackTop() - stacktop_desired;
                if(n_excess > 0 || res_stackpos != stacktop_desired-1)
                {
                    // Remove the cache values
                    synth.DoPopNMov(stacktop_desired-1, res_stackpos);
                }
            }

            if(needs_flip)
                synth.AddOperation(sequencing.op_flip, 1);
        }
    }
}

#endif

#line 1 "fpoptimizer/fpoptimizer_codetree.c"
#include <list>
#include <algorithm>

// line removed for fpoptimizer.c: #include "fpoptimizer_codetree.h"
#include "fptypes.h"
// line removed for fpoptimizer.c: #include "fpoptimizer_consts.h"

#ifdef FP_SUPPORT_OPTIMIZER

using namespace FUNCTIONPARSERTYPES;
//using namespace FPoptimizer_Grammar;

namespace FPoptimizer_CodeTree
{
    CodeTree::CodeTree()
        : data(new CodeTreeData)
    {
        data->Opcode = cNop;
    }

    CodeTree::CodeTree(double i)
        : data(new CodeTreeData(i))
    {
        data->Recalculate_Hash_NoRecursion();
    }

    CodeTree::CodeTree(unsigned v, CodeTree::VarTag)
        : data(new CodeTreeData)
    {
        data->Opcode = VarBegin;
        data->Var    = v;
        data->Recalculate_Hash_NoRecursion();
    }

    CodeTree::CodeTree(const CodeTree& b, CodeTree::CloneTag)
        : data(new CodeTreeData(*b.data))
    {
    }

    CodeTree::~CodeTree()
    {
    }

    struct ParamComparer
    {
        bool operator() (const CodeTree& a, const CodeTree& b) const
        {
            if(a.GetDepth() != b.GetDepth())
                return a.GetDepth() > b.GetDepth();
            return a.GetHash() < b.GetHash();
        }
    };

    void CodeTreeData::Sort()
    {
        /* If the tree is commutative, order the parameters
         * in a set order in order to make equality tests
         * efficient in the optimizer
         */
        switch(Opcode)
        {
            case cAdd:
            case cMul:
            case cMin:
            case cMax:
            case cAnd:
            case cOr:
            case cEqual:
            case cNEqual:
                std::sort(Params.begin(), Params.end(), ParamComparer());
                break;
            case cLess:
                if(ParamComparer() (Params[1], Params[0]))
                    { std::swap(Params[0], Params[1]); Opcode = cGreater; }
                break;
            case cLessOrEq:
                if(ParamComparer() (Params[1], Params[0]))
                    { std::swap(Params[0], Params[1]); Opcode = cGreaterOrEq; }
                break;
            case cGreater:
                if(ParamComparer() (Params[1], Params[0]))
                    { std::swap(Params[0], Params[1]); Opcode = cLess; }
                break;
            case cGreaterOrEq:
                if(ParamComparer() (Params[1], Params[0]))
                    { std::swap(Params[0], Params[1]); Opcode = cLessOrEq; }
                break;
            default:
                break;
        }
    }

    void CodeTree::AddParam(const CodeTree& param)
    {
        //std::cout << "AddParam called\n";
        data->Params.push_back(param);
    }
    void CodeTree::AddParamMove(CodeTree& param)
    {
        data->Params.push_back(CodeTree());
        data->Params.back().swap(param);
    }
    void CodeTree::SetParam(size_t which, const CodeTree& b)
    {
        DataP slot_holder ( data->Params[which].data );
        data->Params[which] = b;
    }
    void CodeTree::SetParamMove(size_t which, CodeTree& b)
    {
        DataP slot_holder ( data->Params[which].data );
        data->Params[which].swap(b);
    }

    void CodeTree::AddParams(const std::vector<CodeTree>& RefParams)
    {
        data->Params.insert(data->Params.end(), RefParams.begin(), RefParams.end());
    }
    void CodeTree::AddParamsMove(std::vector<CodeTree>& RefParams)
    {
        size_t endpos = data->Params.size(), added = RefParams.size();
        data->Params.resize(endpos + added, CodeTree());
        for(size_t p=0; p<added; ++p)
            data->Params[endpos+p].swap( RefParams[p] );
    }
    void CodeTree::AddParamsMove(std::vector<CodeTree>& RefParams, size_t replacing_slot)
    {
        DataP slot_holder ( data->Params[replacing_slot].data );
        DelParam(replacing_slot);
        AddParamsMove(RefParams);
    /*
        const size_t n_added = RefParams.size();
        const size_t oldsize = data->Params.size();
        const size_t newsize = oldsize + n_added - 1;
        if(RefParams.empty())
            DelParam(replacing_slot);
        else
        {
            //    0 1 2 3 4 5 6 7 8 9 10 11
            //    a a a a X b b b b b
            //    a a a a Y Y Y b b b b  b
            //
            //   replacing_slot = 4
            //   n_added = 3
            //   oldsize = 10
            //   newsize = 12
            //   tail_length = 5

            data->Params.resize(newsize);
            data->Params[replacing_slot].data = 0;
            const size_t tail_length = oldsize - replacing_slot -1;
            for(size_t tail=0; tail<tail_length; ++tail)
                data->Params[newsize-1-tail].data.UnsafeSetP(
                &*data->Params[newsize-1-tail-(n_added-1)].data);
            for(size_t head=1; head<n_added; ++head)
                data->Params[replacing_slot+head].data.UnsafeSetP( 0 );
            for(size_t p=0; p<n_added; ++p)
                data->Params[replacing_slot+p].swap( RefParams[p] );
        }
    */
    }

    void CodeTree::SetParams(const std::vector<CodeTree>& RefParams)
    {
        //std::cout << "SetParams called" << (do_clone ? ", clone" : ", no clone") << "\n";
        std::vector<CodeTree> tmp(RefParams);
        data->Params.swap(tmp);
    }

    void CodeTree::SetParamsMove(std::vector<CodeTree>& RefParams)
    {
        data->Params.swap(RefParams);
        RefParams.clear();
    }

#ifdef __GXX_EXPERIMENTAL_CXX0X__
    void CodeTree::SetParams(std::vector<CodeTree>&& RefParams)
    {
        //std::cout << "SetParams&& called\n";
        SetParamsMove(RefParams);
    }
#endif

    void CodeTree::DelParam(size_t index)
    {
        std::vector<CodeTree>& Params = data->Params;
        //std::cout << "DelParam(" << index << ") called\n";
    #ifdef __GXX_EXPERIMENTAL_CXX0X__
        /* rvalue reference semantics makes this optimal */
        Params.erase( Params.begin() + index );
    #else
        /* This labor evades the need for refcount +1/-1 shuffling */
        Params[index].data = 0;
        for(size_t p=index; p+1<Params.size(); ++p)
            Params[p].data.UnsafeSetP( &*Params[p+1].data );
        Params[Params.size()-1].data.UnsafeSetP( 0 );
        Params.resize(Params.size()-1);
    #endif
    }

    void CodeTree::DelParams()
    {
        data->Params.clear();
    }

    /* Is the value of this tree definitely odd(true) or even(false)? */
    CodeTree::TriTruthValue CodeTree::GetEvennessInfo() const
    {
        if(!IsImmed()) return Unknown;
        if(!IsLongIntegerImmed()) return Unknown;
        return (GetLongIntegerImmed() & 1) ? IsNever : IsAlways;
    }

    bool CodeTree::IsLogicalValue() const
    {
        switch(data->Opcode)
        {
            case cImmed:
                return FloatEqual(data->Value, 0.0)
                    || FloatEqual(data->Value, 1.0);
            case cAnd:
            case cOr:
            case cNot:
            case cNotNot:
            case cAbsAnd:
            case cAbsOr:
            case cAbsNot:
            case cAbsNotNot:
            case cEqual:
            case cNEqual:
            case cLess:
            case cLessOrEq:
            case cGreater:
            case cGreaterOrEq:
                /* These operations always produce truth values (0 or 1) */
                return true;
            case cMul:
            {
                std::vector<CodeTree>& Params = data->Params;
                for(size_t a=0; a<Params.size(); ++a)
                    if(!Params[a].IsLogicalValue())
                        return false;
                return true;
            }
            case cIf:
            case cAbsIf:
            {
                std::vector<CodeTree>& Params = data->Params;
                return Params[1].IsLogicalValue()
                    && Params[2].IsLogicalValue();
            }
            default:
                break;
        }
        return false; // Not a logical value.
    }

    bool CodeTree::IsAlwaysInteger(bool integer) const
    {
        switch(data->Opcode)
        {
            case cImmed:
                return IsLongIntegerImmed() ? integer==true : integer==false;
            case cFloor:
            case cInt:
                return integer==true;
            case cAnd:
            case cOr:
            case cNot:
            case cNotNot:
            case cEqual:
            case cNEqual:
            case cLess:
            case cLessOrEq:
            case cGreater:
            case cGreaterOrEq:
                /* These operations always produce truth values (0 or 1) */
                return integer==true; /* 0 and 1 are both integers */
            case cIf:
            {
                std::vector<CodeTree>& Params = data->Params;
                return Params[1].IsAlwaysInteger(integer)
                    && Params[2].IsAlwaysInteger(integer);
                return true; /* 0 and 1 are both integers */
            }
            case cAdd:
            case cMul:
            {
                for(size_t a=GetParamCount(); a-- > 0; )
                    if(!GetParam(a).IsAlwaysInteger(integer))
                        return false;
                return true;
            }
            default:
                break;
        }
        return false; /* Don't know whether it's integer. */
    }

    bool CodeTree::IsAlwaysSigned(bool positive) const
    {
        MinMaxTree tmp = CalculateResultBoundaries();

        if(positive)
            return tmp.has_min && tmp.min >= 0.0
              && (!tmp.has_max || tmp.max >= 0.0);
        else
            return tmp.has_max && tmp.max < 0.0
              && (!tmp.has_min || tmp.min < 0.0);
    }

    bool CodeTree::IsIdenticalTo(const CodeTree& b) const
    {
        //if((!&*data) != (!&*b.data)) return false;
        if(&*data == &*b.data) return true;
        return data->IsIdenticalTo(*b.data);
    }

    bool CodeTreeData::IsIdenticalTo(const CodeTreeData& b) const
    {
        if(Hash   != b.Hash) return false; // a quick catch-all
        if(Opcode != b.Opcode) return false;
        switch(Opcode)
        {
            case cImmed:   return FloatEqual(Value, b.Value);
            case VarBegin: return Var == b.Var;
            case cFCall:
            case cPCall:   if(Funcno != b.Funcno) return false; break;
            default: break;
        }
        if(Params.size() != b.Params.size()) return false;
        for(size_t a=0; a<Params.size(); ++a)
        {
            if(!Params[a].IsIdenticalTo(b.Params[a])) return false;
        }
        return true;
    }

    void CodeTree::Become(const CodeTree& b)
    {
        if(&b != this && &*data != &*b.data)
        {
            DataP tmp = b.data;
            CopyOnWrite();
            data.swap(tmp);
        }
    }

    void CodeTree::CopyOnWrite()
    {
        if(data->RefCount > 1)
            data = new CodeTreeData(*data);
    }

    CodeTree CodeTree::GetUniqueRef()
    {
        if(data->RefCount > 1)
            return CodeTree(*this, CloneTag());
        return *this;
    }

    CodeTreeData::CodeTreeData()
        : RefCount(0),
          Opcode(cNop), Params(), Hash(), Depth(1), OptimizedUsing(0)
    {
    }

    CodeTreeData::CodeTreeData(const CodeTreeData& b)
        : RefCount(0),
          Opcode(b.Opcode),
          Params(b.Params),
          Hash(b.Hash),
          Depth(b.Depth),
          OptimizedUsing(b.OptimizedUsing)
    {
        switch(Opcode)
        {
            case VarBegin: Var   = b.Var; break;
            case cImmed:   Value = b.Value; break;
            case cPCall:
            case cFCall:   Funcno = b.Funcno; break;
            default: break;
        }
    }

#ifdef __GXX_EXPERIMENTAL_CXX0X__
    CodeTreeData::CodeTreeData(CodeTreeData&& b)
        : RefCount(0),
          Opcode(b.Opcode),
          Params(b.Params),
          Hash(b.Hash),
          Depth(b.Depth),
          OptimizedUsing(b.OptimizedUsing)
    {
        switch(Opcode)
        {
            case VarBegin: Var   = b.Var; break;
            case cImmed:   Value = b.Value; break;
            case cPCall:
            case cFCall:   Funcno = b.Funcno; break;
            default: break;
        }
    }
#endif

    CodeTreeData::CodeTreeData(double i)
        : RefCount(0), Opcode(cImmed), Params(), Hash(), Depth(1), OptimizedUsing(0)
    {
        Value = i;
    }
}

#endif

#line 1 "fpoptimizer/fpoptimizer_debug.c"
// line removed for fpoptimizer.c: #include "fpoptimizer_codetree.h"
// line removed for fpoptimizer.c: #include "fpoptimizer_opcodename.h"

#ifdef FP_SUPPORT_OPTIMIZER

#include <sstream>
#include <string>
#include <map>
#include <set>
#include <iostream>

using namespace FUNCTIONPARSERTYPES;

#ifdef FUNCTIONPARSER_SUPPORT_DEBUG_OUTPUT
namespace
{
    void DumpHashesFrom(
        const FPoptimizer_CodeTree::CodeTree& tree,
        std::map<fphash_t, std::set<std::string> >& done,
        std::ostream& o)
    {
        for(size_t a=0; a<tree.GetParamCount(); ++a)
            DumpHashesFrom(tree.GetParam(a), done, o);

        std::ostringstream buf;
        DumpTree(tree, buf);
        done[tree.GetHash()].insert(buf.str());
    }
}
#endif

namespace FPoptimizer_CodeTree
{
#ifdef FUNCTIONPARSER_SUPPORT_DEBUG_OUTPUT
    void DumpHashes(const CodeTree& tree, std::ostream& o)
    {
        std::map<fphash_t, std::set<std::string> > done;
        DumpHashesFrom(tree, done, o);

        for(std::map<fphash_t, std::set<std::string> >::const_iterator
            i = done.begin();
            i != done.end();
            ++i)
        {
            const std::set<std::string>& flist = i->second;
            if(flist.size() != 1) o << "ERROR - HASH COLLISION?\n";
            for(std::set<std::string>::const_iterator
                j = flist.begin();
                j != flist.end();
                ++j)
            {
                o << '[' << std::hex << i->first.hash1
                              << ',' << i->first.hash2
                              << ']' << std::dec;
                o << ": " << *j << "\n";
            }
        }
    }

    void DumpTree(const CodeTree& tree, std::ostream& o)
    {
        //o << "/*" << tree.Depth << "*/";
        const char* sep2 = " ";
        /*
        o << '[' << std::hex << tree.Hash.hash1
                      << ',' << tree.Hash.hash2
                      << ']' << std::dec;
        */
        switch(tree.GetOpcode())
        {
            case cImmed:
                o << tree.GetImmed();
                /*
                o << "(" << std::hex
                  << *(const uint_least64_t*)&tree.GetImmed()
                  << std::dec << ")";
                */
                return;
            case VarBegin: o << "Var" << (tree.GetVar() - VarBegin); return;
            case cAdd: sep2 = " +"; break;
            case cMul: sep2 = " *"; break;
            case cAnd: sep2 = " &"; break;
            case cOr: sep2 = " |"; break;
            case cPow: sep2 = " ^"; break;
            default:
                sep2 = " ";
                o << FP_GetOpcodeName(tree.GetOpcode());
                if(tree.GetOpcode() == cFCall || tree.GetOpcode() == cPCall)
                    o << ':' << tree.GetFuncNo();
        }
        o << '(';
        if(tree.GetParamCount() <= 1 && sep2[1]) o << (sep2+1) << ' ';
        for(size_t a=0; a<tree.GetParamCount(); ++a)
        {
            if(a > 0) o << ' ';

            DumpTree(tree.GetParam(a), o);

            if(a+1 < tree.GetParamCount()) o << sep2;
        }
        o << ')';
    }

    void DumpTreeWithIndent(const CodeTree& tree, std::ostream& o, const std::string& indent)
    {
        o << '[' << std::hex << (void*)(&tree.GetParams())
                 << std::dec
                 << ',' << tree.GetRefCount()
                 << ']';
        o << indent << '_';

        switch(tree.GetOpcode())
        {
            case cImmed:   o << "cImmed " << tree.GetImmed(); o << '\n'; return;
            case VarBegin: o << "VarBegin " << (tree.GetVar() - VarBegin); o << '\n'; return;
            default:
                o << FP_GetOpcodeName(tree.GetOpcode());
                if(tree.GetOpcode() == cFCall || tree.GetOpcode() == cPCall)
                    o << ':' << tree.GetFuncNo();
                o << '\n';
        }
        for(size_t a=0; a<tree.GetParamCount(); ++a)
        {
            std::string ind = indent;
            for(size_t p=0; p < ind.size(); p+=2)
                if(ind[p] == '\\')
                    ind[p] = ' ';
            ind += (a+1 < tree.GetParamCount()) ? " |" : " \\";
            DumpTreeWithIndent(tree.GetParam(a), o, ind);
        }
        o << std::flush;
    }
#endif
}

#endif

#line 1 "fpoptimizer/fpoptimizer_grammar.c"
#include "fparser.h"
#include "fptypes.h"
// line removed for fpoptimizer.c: #include "fpoptimizer_grammar.h"
// line removed for fpoptimizer.c: #include "fpoptimizer_optimize.h"

// line removed for fpoptimizer.c: #include "fpoptimizer_opcodename.h"
using namespace FPoptimizer_Grammar;
using namespace FUNCTIONPARSERTYPES;

#include <cctype>

namespace FPoptimizer_Grammar
{
    bool ParamSpec_Compare(const void* aa, const void* bb, SpecialOpcode type)
    {
        switch(type)
        {
            case ParamHolder:
            {
                ParamSpec_ParamHolder& a = *(ParamSpec_ParamHolder*) aa;
                ParamSpec_ParamHolder& b = *(ParamSpec_ParamHolder*) bb;
                return a.constraints == b.constraints
                    && a.index       == b.index
                    && a.depcode     == b.depcode;
            }
            case NumConstant:
            {
                ParamSpec_NumConstant& a = *(ParamSpec_NumConstant*) aa;
                ParamSpec_NumConstant& b = *(ParamSpec_NumConstant*) bb;
                return FloatEqual(a.constvalue, b.constvalue);
            }
            case SubFunction:
            {
                ParamSpec_SubFunction& a = *(ParamSpec_SubFunction*) aa;
                ParamSpec_SubFunction& b = *(ParamSpec_SubFunction*) bb;
                return a.constraints    == b.constraints
                    && a.data.subfunc_opcode   == b.data.subfunc_opcode
                    && a.data.match_type       == b.data.match_type
                    && a.data.param_count      == b.data.param_count
                    && a.data.param_list       == b.data.param_list
                    && a.data.restholder_index == b.data.restholder_index
                    && a.depcode               == b.depcode;
            }
        }
        return true;
    }

    unsigned ParamSpec_GetDepCode(const ParamSpec& b)
    {
        switch(b.first)
        {
            case ParamHolder:
            {
                const ParamSpec_ParamHolder* s = (const ParamSpec_ParamHolder*) b.second;
                return s->depcode;
            }
            case SubFunction:
            {
                const ParamSpec_SubFunction* s = (const ParamSpec_SubFunction*) b.second;
                return s->depcode;
            }
            default: break;
        }
        return 0;
    }

    void DumpParam(const ParamSpec& parampair, std::ostream& o)
    {
        static const char ParamHolderNames[][2]  = {"%","&", "x","y","z","a","b","c"};

        //o << "/*p" << (&p-pack.plist) << "*/";
        unsigned constraints = 0;
        switch(parampair.first)
        {
            case NumConstant:
              { const ParamSpec_NumConstant& param = *(const ParamSpec_NumConstant*) parampair.second;
                o.precision(12);
                o << param.constvalue; break; }
            case ParamHolder:
              { const ParamSpec_ParamHolder& param = *(const ParamSpec_ParamHolder*) parampair.second;
                o << ParamHolderNames[param.index];
                constraints = param.constraints;
                break; }
            case SubFunction:
              { const ParamSpec_SubFunction& param = *(const ParamSpec_SubFunction*) parampair.second;
                constraints = param.constraints;
                if(param.data.match_type == GroupFunction)
                {
                    if(param.data.subfunc_opcode == cNeg)
                        { o << "-"; DumpParams(param.data.param_list, param.data.param_count, o); }
                    else if(param.data.subfunc_opcode == cInv)
                        { o << "/"; DumpParams(param.data.param_list, param.data.param_count, o); }
                    else
                    {
                        std::string opcode = FP_GetOpcodeName(param.data.subfunc_opcode).substr(1);
                        for(size_t a=0; a<opcode.size(); ++a) opcode[a] = (char) std::toupper(opcode[a]);
                        o << opcode << "( ";
                        DumpParams(param.data.param_list, param.data.param_count, o);
                        o << " )";
                    }
                }
                else
                {
                    o << '(' << FP_GetOpcodeName(param.data.subfunc_opcode) << ' ';
                    if(param.data.match_type == PositionalParams) o << '[';
                    if(param.data.match_type == SelectedParams) o << '{';
                    DumpParams(param.data.param_list, param.data.param_count, o);
                    if(param.data.restholder_index != 0)
                        o << " <" << param.data.restholder_index << '>';
                    if(param.data.match_type == PositionalParams) o << "]";
                    if(param.data.match_type == SelectedParams) o << "}";
                    o << ')';
                }
                break; }
        }
        switch( ImmedConstraint_Value(constraints & ValueMask) )
        {
            case ValueMask: break;
            case Value_AnyNum: break;
            case Value_EvenInt:   o << "@E"; break;
            case Value_OddInt:    o << "@O"; break;
            case Value_IsInteger: o << "@I"; break;
            case Value_NonInteger:o << "@F"; break;
            case Value_Logical:   o << "@L"; break;
        }
        switch( ImmedConstraint_Sign(constraints & SignMask) )
        {
            case SignMask: break;
            case Sign_AnySign: break;
            case Sign_Positive:   o << "@P"; break;
            case Sign_Negative:   o << "@N"; break;
        }
        switch( ImmedConstraint_Oneness(constraints & OnenessMask) )
        {
            case OnenessMask: break;
            case Oneness_Any: break;
            case Oneness_One:     o << "@1"; break;
            case Oneness_NotOne:  o << "@M"; break;
        }
    }

    void DumpParams(unsigned paramlist, unsigned count, std::ostream& o)
    {
        for(unsigned a=0; a<count; ++a)
        {
            if(a > 0) o << ' ';
            const ParamSpec& param = ParamSpec_Extract(paramlist,a);
            DumpParam(param, o);
            unsigned depcode = ParamSpec_GetDepCode(param);
            if(depcode != 0)
                o << "@D" << depcode;
        }
    }
}

#line 1 "fpoptimizer/fpoptimizer_grammar_data.c"
/* This file is automatically generated. Do not edit... */
// line removed for fpoptimizer.c: #include "fpoptimizer_consts.h"
#include "fpconfig.h"
#include "fptypes.h"
#include <algorithm>

#define P1(a) a
#define P2(a,b) (P1(a) | (b << PARAM_INDEX_BITS))
#define P3(a,b,c) (P2(a,b) | (c << (PARAM_INDEX_BITS*2)))

#define grammar_optimize_abslogical grammar_optimize_abslogical_tweak
#define grammar_optimize_ignore_if_sideeffects grammar_optimize_ignore_if_sideeffects_tweak
#define grammar_optimize_nonshortcut_logical_evaluation grammar_optimize_nonshortcut_logical_evaluation_tweak
#define grammar_optimize_round1 grammar_optimize_round1_tweak
#define grammar_optimize_round2 grammar_optimize_round2_tweak
#define grammar_optimize_round3 grammar_optimize_round3_tweak
#define grammar_optimize_round4 grammar_optimize_round4_tweak
#define grammar_optimize_shortcut_logical_evaluation grammar_optimize_shortcut_logical_evaluation_tweak
// line removed for fpoptimizer.c: #include "fpoptimizer_grammar.h"
#undef grammar_optimize_abslogical
#undef grammar_optimize_ignore_if_sideeffects
#undef grammar_optimize_nonshortcut_logical_evaluation
#undef grammar_optimize_round1
#undef grammar_optimize_round2
#undef grammar_optimize_round3
#undef grammar_optimize_round4
#undef grammar_optimize_shortcut_logical_evaluation

using namespace FPoptimizer_Grammar;
using namespace FUNCTIONPARSERTYPES;

namespace
{
    const struct ParamSpec_List
    {
        ParamSpec_ParamHolder plist_p[33];
#define P(n) (n)
        ParamSpec_NumConstant plist_n[17];
#define N(n) (n+33)
        ParamSpec_SubFunction plist_s[422];
#define S(n) (n+33+17)
    } /*PACKED_GRAMMAR_ATTRIBUTE*/ plist =
    {
        { /* plist_p - ParamSpec_ParamHolder[33] */
        /* 0	*/ {0, Sign_Negative | Constness_Const, 0x0}, /* %@N */
        /* 1	*/ {0, Constness_Const, 0x0}, /* % */
        /* 2	*/ {0, Sign_Positive | Constness_Const, 0x0}, /* %@P */
        /* 3	*/ {0, Constness_Const, 0x1}, /* % */
        /* 4	*/ {0, Value_IsInteger | Sign_Positive | Constness_Const, 0x0}, /* %@I@P */
        /* 5	*/ {0, Oneness_NotOne | Constness_Const, 0x1}, /* %@M */
        /* 6	*/ {0, Oneness_NotOne | Constness_Const, 0x0}, /* %@M */
        /* 7	*/ {0, Oneness_One | Constness_Const, 0x0}, /* %@1 */
        /* 8	*/ {0, Value_Logical | Constness_Const, 0x0}, /* %@L */
        /* 9	*/ {1, Constness_Const, 0x0}, /* & */
        /* 10	*/ {1, Value_EvenInt | Constness_Const, 0x0}, /* &@E */
        /* 11	*/ {1, Oneness_NotOne | Constness_Const, 0x0}, /* &@M */
        /* 12	*/ {1, Value_IsInteger | Constness_Const, 0x0}, /* &@I */
        /* 13	*/ {1, Sign_Positive | Constness_Const, 0x0}, /* &@P */
        /* 14	*/ {2, 0, 0x0}, /* x */
        /* 15	*/ {2, 0, 0x4}, /* x */
        /* 16	*/ {2, Value_IsInteger, 0x0}, /* x@I */
        /* 17	*/ {2, Sign_Positive, 0x0}, /* x@P */
        /* 18	*/ {2, Sign_NoIdea, 0x0}, /* x */
        /* 19	*/ {2, Value_Logical, 0x0}, /* x@L */
        /* 20	*/ {3, Sign_NoIdea, 0x0}, /* y */
        /* 21	*/ {3, 0, 0x0}, /* y */
        /* 22	*/ {3, Value_Logical, 0x0}, /* y@L */
        /* 23	*/ {3, 0, 0x8}, /* y */
        /* 24	*/ {3, Value_OddInt, 0x0}, /* y@O */
        /* 25	*/ {3, Value_NonInteger, 0x0}, /* y@F */
        /* 26	*/ {3, Value_EvenInt, 0x0}, /* y@E */
        /* 27	*/ {3, Sign_Positive, 0x0}, /* y@P */
        /* 28	*/ {4, 0, 0x0}, /* z */
        /* 29	*/ {4, 0, 0x16}, /* z */
        /* 30	*/ {4, Value_IsInteger, 0x0}, /* z@I */
        /* 31	*/ {5, 0, 0x0}, /* a */
        /* 32	*/ {6, 0, 0x0}, /* b */
        },

        { /* plist_n - ParamSpec_NumConstant[17] */
        /* 0	*/ {-2}, /* -2 */
        /* 1	*/ {-CONSTANT_PIHALF}, /* -1.57079632679 */
        /* 2	*/ {-1}, /* -1 */
        /* 3	*/ {-0.5}, /* -0.5 */
        /* 4	*/ {0}, /* 0 */
        /* 5	*/ {CONSTANT_RD}, /* 0.0174532925199 */
        /* 6	*/ {CONSTANT_EI}, /* 0.367879441171 */
        /* 7	*/ {CONSTANT_L10I}, /* 0.434294481903 */
        /* 8	*/ {0.5}, /* 0.5 */
        /* 9	*/ {CONSTANT_L2}, /* 0.69314718056 */
        /* 10	*/ {1}, /* 1 */
        /* 11	*/ {CONSTANT_L2I}, /* 1.44269504089 */
        /* 12	*/ {CONSTANT_PIHALF}, /* 1.57079632679 */
        /* 13	*/ {2}, /* 2 */
        /* 14	*/ {CONSTANT_L10}, /* 2.30258509299 */
        /* 15	*/ {CONSTANT_E}, /* 2.71828182846 */
        /* 16	*/ {CONSTANT_DR}, /* 57.2957795131 */
        },

        { /* plist_s - ParamSpec_SubFunction[422] */
        /* 0	*/ {{1,P1(P(14))               , cAbs        ,PositionalParams,0}, 0, 0x0}, /* (cAbs [x]) */
        /* 1	*/ {{1,P1(P(21))               , cAbs        ,PositionalParams,0}, 0, 0x0}, /* (cAbs [y]) */
        /* 2	*/ {{1,P1(S(283))              , cAbs        ,PositionalParams,0}, 0, 0x0}, /* (cAbs [(cMul {x y})]) */
        /* 3	*/ {{1,P1(P(14))               , cAcos       ,PositionalParams,0}, 0, 0x0}, /* (cAcos [x]) */
        /* 4	*/ {{1,P1(P(14))               , cAcosh      ,PositionalParams,0}, 0, 0x0}, /* (cAcosh [x]) */
        /* 5	*/ {{1,P1(P(14))               , cAsin       ,PositionalParams,0}, 0, 0x0}, /* (cAsin [x]) */
        /* 6	*/ {{1,P1(P(14))               , cAsinh      ,PositionalParams,0}, 0, 0x0}, /* (cAsinh [x]) */
        /* 7	*/ {{1,P1(S(242))              , cAsinh      ,PositionalParams,0}, 0, 0x0}, /* (cAsinh [(cAdd  <1>)]) */
        /* 8	*/ {{1,P1(P(14))               , cAtan       ,PositionalParams,0}, 0, 0x0}, /* (cAtan [x]) */
        /* 9	*/ {{2,P2(P(14),P(21))         , cAtan2      ,PositionalParams,0}, 0, 0x0}, /* (cAtan2 [x y]) */
        /* 10	*/ {{2,P2(P(14),S(109))        , cAtan2      ,PositionalParams,0}, 0, 0x0}, /* (cAtan2 [x (cPow [y -%])]) */
        /* 11	*/ {{1,P1(P(14))               , cAtanh      ,PositionalParams,0}, 0, 0x0}, /* (cAtanh [x]) */
        /* 12	*/ {{1,P1(P(14))               , cCeil       ,PositionalParams,0}, 0, 0x4}, /* (cCeil [x]) */
        /* 13	*/ {{1,P1(S(331))              , cCeil       ,PositionalParams,0}, 0, 0x0}, /* (cCeil [(cMul  <1>)]) */
        /* 14	*/ {{1,P1(P(14))               , cCos        ,PositionalParams,0}, 0, 0x0}, /* (cCos [x]) */
        /* 15	*/ {{1,P1(P(14))               , cCos        ,PositionalParams,0}, 0, 0x4}, /* (cCos [x]) */
        /* 16	*/ {{1,P1(P(21))               , cCos        ,PositionalParams,0}, 0, 0x0}, /* (cCos [y]) */
        /* 17	*/ {{1,P1(S(215))              , cCos        ,PositionalParams,0}, 0, 0x0}, /* (cCos [(cAdd {x y})]) */
        /* 18	*/ {{1,P1(S(221))              , cCos        ,PositionalParams,0}, 0, 0x0}, /* (cCos [(cAdd {x (cMul {-1 y})})]) */
        /* 19	*/ {{1,P1(S(242))              , cCos        ,PositionalParams,0}, 0, 0x0}, /* (cCos [(cAdd  <1>)]) */
        /* 20	*/ {{1,P1(S(301))              , cCos        ,PositionalParams,0}, 0, 0x0}, /* (cCos [(cMul {-% x})]) */
        /* 21	*/ {{1,P1(S(358))              , cCos        ,PositionalParams,0}, 0, 0x0}, /* (cCos [(cMul -% <1>)]) */
        /* 22	*/ {{1,P1(P(14))               , cCosh       ,PositionalParams,0}, 0, 0x4}, /* (cCosh [x]) */
        /* 23	*/ {{1,P1(P(14))               , cCosh       ,PositionalParams,0}, 0, 0x0}, /* (cCosh [x]) */
        /* 24	*/ {{1,P1(S(50))               , cCosh       ,PositionalParams,0}, 0, 0x0}, /* (cCosh [(cLog [(cPow [& x])])]) */
        /* 25	*/ {{1,P1(S(295))              , cCosh       ,PositionalParams,0}, 0, 0x0}, /* (cCosh [(cMul {x LOG( & )})]) */
        /* 26	*/ {{1,P1(S(301))              , cCosh       ,PositionalParams,0}, 0, 0x0}, /* (cCosh [(cMul {-% x})]) */
        /* 27	*/ {{1,P1(P(14))               , cCot        ,PositionalParams,0}, 0, 0x0}, /* (cCot [x]) */
        /* 28	*/ {{1,P1(P(14))               , cCsc        ,PositionalParams,0}, 0, 0x0}, /* (cCsc [x]) */
        /* 29	*/ {{1,P1(P(14))               , cFloor      ,PositionalParams,0}, 0, 0x4}, /* (cFloor [x]) */
        /* 30	*/ {{1,P1(S(331))              , cFloor      ,PositionalParams,0}, 0, 0x0}, /* (cFloor [(cMul  <1>)]) */
        /* 31	*/ {{3,P3(P(14),N(10),S(400))  , cIf         ,PositionalParams,0}, 0, 0x0}, /* (cIf [x 1 (cNotNot [y])]) */
        /* 32	*/ {{3,P3(P(14),P(21),P(28))   , cIf         ,PositionalParams,0}, 0, 0x4}, /* (cIf [x y z]) */
        /* 33	*/ {{3,P3(P(14),P(21),P(28))   , cIf         ,PositionalParams,0}, 0, 0x0}, /* (cIf [x y z]) */
        /* 34	*/ {{3,P3(P(14),P(8),S(376))   , cIf         ,PositionalParams,0}, 0, 0x0}, /* (cIf [x %@L (cNot [y])]) */
        /* 35	*/ {{3,P3(P(14),P(21),S(374))  , cIf         ,PositionalParams,0}, 0, 0x0}, /* (cIf [x y (cNot [%])]) */
        /* 36	*/ {{3,P3(P(14),P(31),P(32))   , cIf         ,PositionalParams,0}, 0, 0x4}, /* (cIf [x a b]) */
        /* 37	*/ {{3,P3(P(14),S(374),P(21))  , cIf         ,PositionalParams,0}, 0, 0x0}, /* (cIf [x (cNot [%]) y]) */
        /* 38	*/ {{3,P3(P(14),S(376),P(8))   , cIf         ,PositionalParams,0}, 0, 0x0}, /* (cIf [x (cNot [y]) %@L]) */
        /* 39	*/ {{3,P3(P(14),S(58),S(59))   , cIf         ,PositionalParams,0}, 0, 0x0}, /* (cIf [x (cMax [y a]) (cMax [z b])]) */
        /* 40	*/ {{3,P3(P(14),S(61),S(62))   , cIf         ,PositionalParams,0}, 0, 0x0}, /* (cIf [x (cMin [y a]) (cMin [z b])]) */
        /* 41	*/ {{3,P3(P(14),S(222),S(224)) , cIf         ,PositionalParams,0}, 0, 0x0}, /* (cIf [x (cAdd {y a}) (cAdd {z b})]) */
        /* 42	*/ {{3,P3(P(14),S(288),S(289)) , cIf         ,PositionalParams,0}, 0, 0x0}, /* (cIf [x (cMul {y a}) (cMul {z b})]) */
        /* 43	*/ {{3,P3(P(14),S(386),S(387)) , cIf         ,PositionalParams,0}, 0, 0x0}, /* (cIf [x (cAnd {y a}) (cAnd {z b})]) */
        /* 44	*/ {{3,P3(P(14),S(393),S(394)) , cIf         ,PositionalParams,0}, 0, 0x0}, /* (cIf [x (cOr {y a}) (cOr {z b})]) */
        /* 45	*/ {{3,P3(P(14),S(400),N(4))   , cIf         ,PositionalParams,0}, 0, 0x0}, /* (cIf [x (cNotNot [y]) 0]) */
        /* 46	*/ {{1,P1(S(242))              , cInt        ,PositionalParams,0}, 0, 0x0}, /* (cInt [(cAdd  <1>)]) */
        /* 47	*/ {{1,P1(P(14))               , cLog        ,PositionalParams,0}, 0, 0x0}, /* (cLog [x]) */
        /* 48	*/ {{1,P1(P(21))               , cLog        ,PositionalParams,0}, 0, 0x0}, /* (cLog [y]) */
        /* 49	*/ {{1,P1(P(28))               , cLog        ,PositionalParams,0}, 0, 0x0}, /* (cLog [z]) */
        /* 50	*/ {{1,P1(S(83))               , cLog        ,PositionalParams,0}, 0, 0x0}, /* (cLog [(cPow [& x])]) */
        /* 51	*/ {{1,P1(S(283))              , cLog        ,PositionalParams,0}, 0, 0x0}, /* (cLog [(cMul {x y})]) */
        /* 52	*/ {{1,P1(S(331))              , cLog        ,PositionalParams,0}, 0, 0x0}, /* (cLog [(cMul  <1>)]) */
        /* 53	*/ {{1,P1(P(1))                , cLog        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* LOG( % ) */
        /* 54	*/ {{1,P1(P(9))                , cLog        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* LOG( & ) */
        /* 55	*/ {{1,P1(P(14))               , cLog10      ,PositionalParams,0}, 0, 0x0}, /* (cLog10 [x]) */
        /* 56	*/ {{1,P1(P(14))               , cLog2       ,PositionalParams,0}, 0, 0x0}, /* (cLog2 [x]) */
        /* 57	*/ {{2,P2(P(14),P(21))         , cMax        ,PositionalParams,0}, 0, 0x0}, /* (cMax [x y]) */
        /* 58	*/ {{2,P2(P(21),P(31))         , cMax        ,PositionalParams,0}, 0, 0x0}, /* (cMax [y a]) */
        /* 59	*/ {{2,P2(P(28),P(32))         , cMax        ,PositionalParams,0}, 0, 0x0}, /* (cMax [z b]) */
        /* 60	*/ {{2,P2(P(14),P(21))         , cMin        ,PositionalParams,0}, 0, 0x0}, /* (cMin [x y]) */
        /* 61	*/ {{2,P2(P(21),P(31))         , cMin        ,PositionalParams,0}, 0, 0x0}, /* (cMin [y a]) */
        /* 62	*/ {{2,P2(P(28),P(32))         , cMin        ,PositionalParams,0}, 0, 0x0}, /* (cMin [z b]) */
        /* 63	*/ {{2,P2(P(1),N(10))          , cMin        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* MIN( % 1 ) */
        /* 64	*/ {{2,P2(P(1),P(9))           , cMin        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* MIN( % & ) */
        /* 65	*/ {{2,P2(N(6),P(14))          , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [0.367879441171 x]) */
        /* 66	*/ {{2,P2(N(6),P(14))          , cPow        ,PositionalParams,0}, 0, 0x4}, /* (cPow [0.367879441171 x]) */
        /* 67	*/ {{2,P2(N(15),P(14))         , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [2.71828182846 x]) */
        /* 68	*/ {{2,P2(N(15),P(14))         , cPow        ,PositionalParams,0}, 0, 0x4}, /* (cPow [2.71828182846 x]) */
        /* 69	*/ {{2,P2(N(6),S(6))           , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [0.367879441171 (cAsinh [x])]) */
        /* 70	*/ {{2,P2(N(15),S(6))          , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [2.71828182846 (cAsinh [x])]) */
        /* 71	*/ {{2,P2(N(6),S(7))           , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [0.367879441171 (cAsinh [(cAdd  <1>)])]) */
        /* 72	*/ {{2,P2(N(15),S(7))          , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [2.71828182846 (cAsinh [(cAdd  <1>)])]) */
        /* 73	*/ {{2,P2(N(15),S(332))        , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [2.71828182846 (cMul  <2>)]) */
        /* 74	*/ {{2,P2(P(1),P(14))          , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [% x]) */
        /* 75	*/ {{2,P2(P(1),P(14))          , cPow        ,PositionalParams,0}, 0, 0x1}, /* (cPow [% x]) */
        /* 76	*/ {{2,P2(P(1),P(21))          , cPow        ,PositionalParams,0}, 0, 0x1}, /* (cPow [% y]) */
        /* 77	*/ {{2,P2(P(14),P(1))          , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [x %]) */
        /* 78	*/ {{2,P2(P(14),P(2))          , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [x %@P]) */
        /* 79	*/ {{2,P2(P(21),P(13))         , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [y &@P]) */
        /* 80	*/ {{2,P2(P(21),P(9))          , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [y &]) */
        /* 81	*/ {{2,P2(P(2),P(14))          , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [%@P x]) */
        /* 82	*/ {{2,P2(P(9),P(14))          , cPow        ,PositionalParams,0}, 0, 0x6}, /* (cPow [& x]) */
        /* 83	*/ {{2,P2(P(9),P(14))          , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [& x]) */
        /* 84	*/ {{2,P2(P(9),P(21))          , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [& y]) */
        /* 85	*/ {{2,P2(P(9),S(223))         , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [& (cAdd {y (cMul {z (cLog [x]) /LOG( & )})})]) */
        /* 86	*/ {{2,P2(P(14),P(4))          , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [x %@I@P]) */
        /* 87	*/ {{2,P2(P(14),P(12))         , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [x &@I]) */
        /* 88	*/ {{2,P2(P(1),P(21))          , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [% y]) */
        /* 89	*/ {{2,P2(P(14),N(13))         , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [x 2]) */
        /* 90	*/ {{2,P2(P(14),P(24))         , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [x y@O]) */
        /* 91	*/ {{2,P2(P(14),P(25))         , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [x y@F]) */
        /* 92	*/ {{2,P2(P(17),P(21))         , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [x@P y]) */
        /* 93	*/ {{2,P2(P(18),P(21))         , cPow        ,PositionalParams,0}, Sign_Positive, 0x0}, /* (cPow [x y])@P */
        /* 94	*/ {{2,P2(P(14),P(21))         , cPow        ,PositionalParams,0}, 0, 0x4}, /* (cPow [x y]) */
        /* 95	*/ {{2,P2(P(14),P(28))         , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [x z]) */
        /* 96	*/ {{2,P2(P(17),P(28))         , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [x@P z]) */
        /* 97	*/ {{2,P2(P(14),S(63))         , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [x MIN( % 1 )]) */
        /* 98	*/ {{2,P2(P(14),S(64))         , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [x MIN( % & )]) */
        /* 99	*/ {{2,P2(P(14),S(208))        , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [x (cAdd {1 -MIN( % 1 )})]) */
        /* 100	*/ {{2,P2(P(14),S(216))        , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [x (cAdd {y z})]) */
        /* 101	*/ {{2,P2(P(14),S(218))        , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [x (cAdd {% -MIN( % 1 )})]) */
        /* 102	*/ {{2,P2(P(14),S(219))        , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [x (cAdd {% -MIN( % & )})]) */
        /* 103	*/ {{2,P2(P(14),S(220))        , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [x (cAdd {& -MIN( % & )})]) */
        /* 104	*/ {{2,P2(P(21),N(2))          , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [y -1]) */
        /* 105	*/ {{2,P2(P(14),P(21))         , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [x y]) */
        /* 106	*/ {{2,P2(P(1),S(242))         , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [% (cAdd  <1>)]) */
        /* 107	*/ {{2,P2(P(14),S(53))         , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [x LOG( % )]) */
        /* 108	*/ {{2,P2(P(20),P(0))          , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [y %@N]) */
        /* 109	*/ {{2,P2(P(21),S(194))        , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [y -%]) */
        /* 110	*/ {{2,P2(P(28),S(242))        , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [z (cAdd  <1>)]) */
        /* 111	*/ {{2,P2(S(14),N(2))          , cPow        ,PositionalParams,0}, 0, 0x4}, /* (cPow [(cCos [x]) -1]) */
        /* 112	*/ {{2,P2(S(14),N(2))          , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cCos [x]) -1]) */
        /* 113	*/ {{2,P2(S(20),N(2))          , cPow        ,PositionalParams,0}, 0, 0x5}, /* (cPow [(cCos [(cMul {-% x})]) -1]) */
        /* 114	*/ {{2,P2(S(21),N(2))          , cPow        ,PositionalParams,0}, 0, 0x1}, /* (cPow [(cCos [(cMul -% <1>)]) -1]) */
        /* 115	*/ {{2,P2(S(23),N(2))          , cPow        ,PositionalParams,0}, 0, 0x4}, /* (cPow [(cCosh [x]) -1]) */
        /* 116	*/ {{2,P2(S(23),N(2))          , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cCosh [x]) -1]) */
        /* 117	*/ {{2,P2(S(14),N(13))         , cPow        ,PositionalParams,0}, 0, 0x4}, /* (cPow [(cCos [x]) 2]) */
        /* 118	*/ {{2,P2(S(14),N(13))         , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cCos [x]) 2]) */
        /* 119	*/ {{2,P2(S(23),P(1))          , cPow        ,PositionalParams,0}, 0, 0x4}, /* (cPow [(cCosh [x]) %]) */
        /* 120	*/ {{2,P2(S(23),S(259))        , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cCosh [x]) ADD( % 1 )]) */
        /* 121	*/ {{2,P2(S(26),N(2))          , cPow        ,PositionalParams,0}, 0, 0x5}, /* (cPow [(cCosh [(cMul {-% x})]) -1]) */
        /* 122	*/ {{2,P2(S(47),N(2))          , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cLog [x]) -1]) */
        /* 123	*/ {{2,P2(S(49),N(2))          , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cLog [z]) -1]) */
        /* 124	*/ {{2,P2(S(55),N(2))          , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cLog10 [x]) -1]) */
        /* 125	*/ {{2,P2(S(56),N(2))          , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cLog2 [x]) -1]) */
        /* 126	*/ {{2,P2(S(77),S(417))        , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cPow [x %]) /%]) */
        /* 127	*/ {{2,P2(S(77),S(421))        , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cPow [x %]) /MIN( % & )]) */
        /* 128	*/ {{2,P2(S(80),S(421))        , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cPow [y &]) /MIN( % & )]) */
        /* 129	*/ {{2,P2(S(158),N(2))         , cPow        ,PositionalParams,0}, 0, 0x4}, /* (cPow [(cSin [x]) -1]) */
        /* 130	*/ {{2,P2(S(170),N(2))         , cPow        ,PositionalParams,0}, 0, 0x4}, /* (cPow [(cSinh [x]) -1]) */
        /* 131	*/ {{2,P2(S(177),N(2))         , cPow        ,PositionalParams,0}, 0, 0x4}, /* (cPow [(cTan [x]) -1]) */
        /* 132	*/ {{2,P2(S(177),N(2))         , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cTan [x]) -1]) */
        /* 133	*/ {{2,P2(S(186),N(2))         , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cTanh [x]) -1]) */
        /* 134	*/ {{2,P2(S(186),N(2))         , cPow        ,PositionalParams,0}, 0, 0x4}, /* (cPow [(cTanh [x]) -1]) */
        /* 135	*/ {{2,P2(S(189),N(2))         , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cTanh [(cMul {x LOG( % ) 0.5})]) -1]) */
        /* 136	*/ {{2,P2(S(199),N(2))         , cPow        ,PositionalParams,0}, 0, 0x5}, /* (cPow [(cAdd {-1 (cPow [% x])}) -1]) */
        /* 137	*/ {{2,P2(S(201),N(2))         , cPow        ,PositionalParams,0}, 0, 0x5}, /* (cPow [(cAdd {1 (cPow [% x])}) -1]) */
        /* 138	*/ {{2,P2(S(332),N(2))         , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cMul  <2>) -1]) */
        /* 139	*/ {{2,P2(S(346),N(2))         , cPow        ,PositionalParams,0}, 0, 0x4}, /* (cPow [(cMul x <2>) -1]) */
        /* 140	*/ {{2,P2(S(158),N(13))        , cPow        ,PositionalParams,0}, 0, 0x4}, /* (cPow [(cSin [x]) 2]) */
        /* 141	*/ {{2,P2(S(158),N(13))        , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cSin [x]) 2]) */
        /* 142	*/ {{2,P2(S(205),N(3))         , cPow        ,PositionalParams,0}, 0, 0x4}, /* (cPow [(cAdd {1 (cPow [x 2])}) -0.5]) */
        /* 143	*/ {{2,P2(S(205),N(8))         , cPow        ,PositionalParams,0}, 0, 0x4}, /* (cPow [(cAdd {1 (cPow [x 2])}) 0.5]) */
        /* 144	*/ {{2,P2(S(207),N(3))         , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cAdd {1 (cPow [(cAdd  <1>) 2])}) -0.5]) */
        /* 145	*/ {{2,P2(S(207),N(8))         , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cAdd {1 (cPow [(cAdd  <1>) 2])}) 0.5]) */
        /* 146	*/ {{2,P2(S(209),N(2))         , cPow        ,PositionalParams,0}, 0, 0x4}, /* (cPow [(cAdd {1 (cMul {-1 x})}) -1]) */
        /* 147	*/ {{2,P2(S(203),N(8))         , cPow        ,PositionalParams,0}, 0, 0x4}, /* (cPow [(cAdd {-1 (cPow [x 2])}) 0.5]) */
        /* 148	*/ {{2,P2(S(226),N(8))         , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cAdd {(cPow [x 2]) -1}) 0.5]) */
        /* 149	*/ {{2,P2(S(236),N(8))         , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cAdd {(cPow [x 2]) 1}) 0.5]) */
        /* 150	*/ {{2,P2(S(237),N(8))         , cPow        ,PositionalParams,0}, 0, 0x4}, /* (cPow [(cAdd {(cMul {(cPow [x 2]) -1}) 1}) 0.5]) */
        /* 151	*/ {{2,P2(S(242),N(13))        , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cAdd  <1>) 2]) */
        /* 152	*/ {{2,P2(S(331),P(9))         , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [(cMul  <1>) &]) */
        /* 153	*/ {{2,P2(S(418),P(14))        , cPow        ,PositionalParams,0}, 0, 0x0}, /* (cPow [/& x]) */
        /* 154	*/ {{2,P2(S(418),P(14))        , cPow        ,PositionalParams,0}, 0, 0x6}, /* (cPow [/& x]) */
        /* 155	*/ {{2,P2(P(1),P(9))           , cPow        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* POW( % & ) */
        /* 156	*/ {{2,P2(P(9),S(417))         , cPow        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* POW( & /% ) */
        /* 157	*/ {{1,P1(P(14))               , cSec        ,PositionalParams,0}, 0, 0x0}, /* (cSec [x]) */
        /* 158	*/ {{1,P1(P(14))               , cSin        ,PositionalParams,0}, 0, 0x0}, /* (cSin [x]) */
        /* 159	*/ {{1,P1(P(14))               , cSin        ,PositionalParams,0}, 0, 0x4}, /* (cSin [x]) */
        /* 160	*/ {{1,P1(P(21))               , cSin        ,PositionalParams,0}, 0, 0x0}, /* (cSin [y]) */
        /* 161	*/ {{1,P1(S(215))              , cSin        ,PositionalParams,0}, 0, 0x0}, /* (cSin [(cAdd {x y})]) */
        /* 162	*/ {{1,P1(S(221))              , cSin        ,PositionalParams,0}, 0, 0x0}, /* (cSin [(cAdd {x (cMul {-1 y})})]) */
        /* 163	*/ {{1,P1(S(242))              , cSin        ,PositionalParams,0}, 0, 0x0}, /* (cSin [(cAdd  <1>)]) */
        /* 164	*/ {{1,P1(S(279))              , cSin        ,PositionalParams,0}, 0, 0x5}, /* (cSin [(cMul {% x})]) */
        /* 165	*/ {{1,P1(S(331))              , cSin        ,PositionalParams,0}, 0, 0x0}, /* (cSin [(cMul  <1>)]) */
        /* 166	*/ {{1,P1(S(335))              , cSin        ,PositionalParams,0}, 0, 0x0}, /* (cSin [(cMul %@N <1>)]) */
        /* 167	*/ {{1,P1(S(339))              , cSin        ,PositionalParams,0}, 0, 0x1}, /* (cSin [(cMul % <1>)]) */
        /* 168	*/ {{1,P1(S(358))              , cSin        ,PositionalParams,0}, 0, 0x0}, /* (cSin [(cMul -% <1>)]) */
        /* 169	*/ {{1,P1(P(14))               , cSinh       ,PositionalParams,0}, 0, 0x4}, /* (cSinh [x]) */
        /* 170	*/ {{1,P1(P(14))               , cSinh       ,PositionalParams,0}, 0, 0x0}, /* (cSinh [x]) */
        /* 171	*/ {{1,P1(S(50))               , cSinh       ,PositionalParams,0}, 0, 0x0}, /* (cSinh [(cLog [(cPow [& x])])]) */
        /* 172	*/ {{1,P1(S(279))              , cSinh       ,PositionalParams,0}, 0, 0x5}, /* (cSinh [(cMul {% x})]) */
        /* 173	*/ {{1,P1(S(295))              , cSinh       ,PositionalParams,0}, 0, 0x0}, /* (cSinh [(cMul {x LOG( & )})]) */
        /* 174	*/ {{1,P1(S(331))              , cSinh       ,PositionalParams,0}, 0, 0x0}, /* (cSinh [(cMul  <1>)]) */
        /* 175	*/ {{1,P1(S(335))              , cSinh       ,PositionalParams,0}, 0, 0x0}, /* (cSinh [(cMul %@N <1>)]) */
        /* 176	*/ {{1,P1(S(358))              , cSinh       ,PositionalParams,0}, 0, 0x0}, /* (cSinh [(cMul -% <1>)]) */
        /* 177	*/ {{1,P1(P(14))               , cTan        ,PositionalParams,0}, 0, 0x0}, /* (cTan [x]) */
        /* 178	*/ {{1,P1(P(14))               , cTan        ,PositionalParams,0}, 0, 0x4}, /* (cTan [x]) */
        /* 179	*/ {{1,P1(S(210))              , cTan        ,PositionalParams,0}, 0, 0x4}, /* (cTan [(cAdd {1.57079632679 (cMul {-1 x})})]) */
        /* 180	*/ {{1,P1(S(213))              , cTan        ,PositionalParams,0}, 0, 0x0}, /* (cTan [(cAdd {1.57079632679 (cMul -1 <1>)})]) */
        /* 181	*/ {{1,P1(S(279))              , cTan        ,PositionalParams,0}, 0, 0x0}, /* (cTan [(cMul {% x})]) */
        /* 182	*/ {{1,P1(S(331))              , cTan        ,PositionalParams,0}, 0, 0x0}, /* (cTan [(cMul  <1>)]) */
        /* 183	*/ {{1,P1(S(339))              , cTan        ,PositionalParams,0}, 0, 0x0}, /* (cTan [(cMul % <1>)]) */
        /* 184	*/ {{1,P1(S(335))              , cTan        ,PositionalParams,0}, 0, 0x0}, /* (cTan [(cMul %@N <1>)]) */
        /* 185	*/ {{1,P1(S(358))              , cTan        ,PositionalParams,0}, 0, 0x0}, /* (cTan [(cMul -% <1>)]) */
        /* 186	*/ {{1,P1(P(14))               , cTanh       ,PositionalParams,0}, 0, 0x0}, /* (cTanh [x]) */
        /* 187	*/ {{1,P1(P(14))               , cTanh       ,PositionalParams,0}, 0, 0x4}, /* (cTanh [x]) */
        /* 188	*/ {{1,P1(S(279))              , cTanh       ,PositionalParams,0}, 0, 0x0}, /* (cTanh [(cMul {% x})]) */
        /* 189	*/ {{1,P1(S(286))              , cTanh       ,PositionalParams,0}, 0, 0x0}, /* (cTanh [(cMul {x LOG( % ) 0.5})]) */
        /* 190	*/ {{1,P1(S(331))              , cTanh       ,PositionalParams,0}, 0, 0x0}, /* (cTanh [(cMul  <1>)]) */
        /* 191	*/ {{1,P1(S(335))              , cTanh       ,PositionalParams,0}, 0, 0x0}, /* (cTanh [(cMul %@N <1>)]) */
        /* 192	*/ {{1,P1(S(358))              , cTanh       ,PositionalParams,0}, 0, 0x0}, /* (cTanh [(cMul -% <1>)]) */
        /* 193	*/ {{1,P1(P(14))               , cTrunc      ,PositionalParams,0}, 0, 0x0}, /* (cTrunc [x]) */
        /* 194	*/ {{1,P1(P(1))                , cNeg        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* -% */
        /* 195	*/ {{1,P1(P(1))                , cNeg        ,GroupFunction   ,0}, Constness_Const, 0x1}, /* -% */
        /* 196	*/ {{1,P1(S(63))               , cNeg        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* -MIN( % 1 ) */
        /* 197	*/ {{1,P1(S(64))               , cNeg        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* -MIN( % & ) */
        /* 198	*/ {{2,P2(N(2),S(74))          , cAdd        ,SelectedParams  ,0}, 0, 0x5}, /* (cAdd {-1 (cPow [% x])}) */
        /* 199	*/ {{2,P2(N(2),S(74))          , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {-1 (cPow [% x])}) */
        /* 200	*/ {{2,P2(N(10),P(14))         , cAdd        ,SelectedParams  ,0}, 0, 0x4}, /* (cAdd {1 x}) */
        /* 201	*/ {{2,P2(N(10),S(74))         , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {1 (cPow [% x])}) */
        /* 202	*/ {{2,P2(N(10),S(74))         , cAdd        ,SelectedParams  ,0}, 0, 0x5}, /* (cAdd {1 (cPow [% x])}) */
        /* 203	*/ {{2,P2(N(2),S(89))          , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {-1 (cPow [x 2])}) */
        /* 204	*/ {{2,P2(N(2),S(96))          , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {-1 (cPow [x@P z])}) */
        /* 205	*/ {{2,P2(N(10),S(89))         , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {1 (cPow [x 2])}) */
        /* 206	*/ {{2,P2(N(10),S(96))         , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {1 (cPow [x@P z])}) */
        /* 207	*/ {{2,P2(N(10),S(151))        , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {1 (cPow [(cAdd  <1>) 2])}) */
        /* 208	*/ {{2,P2(N(10),S(196))        , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {1 -MIN( % 1 )}) */
        /* 209	*/ {{2,P2(N(10),S(261))        , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {1 (cMul {-1 x})}) */
        /* 210	*/ {{2,P2(N(12),S(261))        , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {1.57079632679 (cMul {-1 x})}) */
        /* 211	*/ {{2,P2(N(10),S(303))        , cAdd        ,SelectedParams  ,0}, 0, 0x1}, /* (cAdd {1 (cMul {(cLog [x]) /%})}) */
        /* 212	*/ {{2,P2(N(10),S(334))        , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {1 (cMul -1 <2>)}) */
        /* 213	*/ {{2,P2(N(12),S(333))        , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {1.57079632679 (cMul -1 <1>)}) */
        /* 214	*/ {{2,P2(N(12),S(335))        , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {1.57079632679 (cMul %@N <1>)}) */
        /* 215	*/ {{2,P2(P(14),P(21))         , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {x y}) */
        /* 216	*/ {{2,P2(P(21),P(28))         , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {y z}) */
        /* 217	*/ {{2,P2(P(7),S(95))          , cAdd        ,SelectedParams  ,0}, 0, 0x4}, /* (cAdd {%@1 (cPow [x z])}) */
        /* 218	*/ {{2,P2(P(1),S(196))         , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {% -MIN( % 1 )}) */
        /* 219	*/ {{2,P2(P(1),S(197))         , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {% -MIN( % & )}) */
        /* 220	*/ {{2,P2(P(9),S(197))         , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {& -MIN( % & )}) */
        /* 221	*/ {{2,P2(P(14),S(262))        , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {x (cMul {-1 y})}) */
        /* 222	*/ {{2,P2(P(21),P(31))         , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {y a}) */
        /* 223	*/ {{2,P2(P(21),S(290))        , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {y (cMul {z (cLog [x]) /LOG( & )})}) */
        /* 224	*/ {{2,P2(P(28),P(32))         , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {z b}) */
        /* 225	*/ {{2,P2(S(47),P(1))          , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {(cLog [x]) %}) */
        /* 226	*/ {{2,P2(S(89),N(2))          , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {(cPow [x 2]) -1}) */
        /* 227	*/ {{2,P2(S(280),P(3))         , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {(cMul {% (cPow [x@P z])})@D1 %@D1}) */
        /* 228	*/ {{2,P2(S(47),P(9))          , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {(cLog [x]) &}) */
        /* 229	*/ {{2,P2(S(84),S(85))         , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {(cPow [& y]) (cPow [& (cAdd {y (cMul {z (cLog [x]) /LOG( & )})})])}) */
        /* 230	*/ {{2,P2(S(147),P(15))        , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {(cPow [(cAdd {-1 (cPow [x 2])}) 0.5])@D4 x@D4}) */
        /* 231	*/ {{2,P2(S(280),S(195))       , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {(cMul {% (cPow [x@P z])})@D1 -%@D1}) */
        /* 232	*/ {{2,P2(S(298),S(85))        , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {(cMul {(cPow [& y]) -1}) (cPow [& (cAdd {y (cMul {z (cLog [x]) /LOG( & )})})])}) */
        /* 233	*/ {{2,P2(S(300),S(100))       , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {(cMul {(cPow [x y]) %}) (cPow [x (cAdd {y z})])}) */
        /* 234	*/ {{2,P2(S(316),S(52))        , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {(cMul {LOG( % ) y}) (cLog [(cMul  <1>)])}) */
        /* 235	*/ {{2,P2(S(52),S(53))         , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {(cLog [(cMul  <1>)]) LOG( % )}) */
        /* 236	*/ {{2,P2(S(89),N(10))         , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {(cPow [x 2]) 1}) */
        /* 237	*/ {{2,P2(S(317),N(10))        , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {(cMul {(cPow [x 2]) -1}) 1}) */
        /* 238	*/ {{2,P2(S(341),S(293))       , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {(cMul % & <1>) (cMul {& (cAdd  <2>)})}) */
        /* 239	*/ {{2,P2(S(362),S(294))       , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {MUL( % & ) (cMul {& (cAdd  <1>)})}) */
        /* 240	*/ {{2,P2(S(354),S(353))       , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {(cMul (cPow [x (cAdd {% -MIN( % 1 )})]) <1>) (cMul (cPow [x (cAdd {1 -MIN( % 1 )})]) <2>)}) */
        /* 241	*/ {{2,P2(S(355),S(356))       , cAdd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAdd {(cMul (cPow [x (cAdd {% -MIN( % & )})]) <1>) (cMul (cPow [x (cAdd {& -MIN( % & )})]) <2>)}) */
        /* 242	*/ {{0,0                       , cAdd        ,AnyParams       ,1}, 0, 0x0}, /* (cAdd  <1>) */
        /* 243	*/ {{0,0                       , cAdd        ,AnyParams       ,2}, 0, 0x0}, /* (cAdd  <2>) */
        /* 244	*/ {{0,0                       , cAdd        ,AnyParams       ,1}, Sign_Positive, 0x0}, /* (cAdd  <1>)@P */
        /* 245	*/ {{1,P1(N(1))                , cAdd        ,AnyParams       ,1}, 0, 0x0}, /* (cAdd -1.57079632679 <1>) */
        /* 246	*/ {{1,P1(N(8))                , cAdd        ,AnyParams       ,1}, 0, 0x0}, /* (cAdd 0.5 <1>) */
        /* 247	*/ {{1,P1(P(6))                , cAdd        ,AnyParams       ,1}, 0, 0x0}, /* (cAdd %@M <1>) */
        /* 248	*/ {{1,P1(P(1))                , cAdd        ,AnyParams       ,1}, 0, 0x0}, /* (cAdd % <1>) */
        /* 249	*/ {{1,P1(P(11))               , cAdd        ,AnyParams       ,1}, 0, 0x0}, /* (cAdd &@M <1>) */
        /* 250	*/ {{1,P1(P(9))                , cAdd        ,AnyParams       ,2}, 0, 0x0}, /* (cAdd & <2>) */
        /* 251	*/ {{1,P1(P(14))               , cAdd        ,AnyParams       ,1}, 0, 0x4}, /* (cAdd x <1>) */
        /* 252	*/ {{1,P1(P(14))               , cAdd        ,AnyParams       ,2}, 0, 0x4}, /* (cAdd x <2>) */
        /* 253	*/ {{1,P1(P(14))               , cAdd        ,AnyParams       ,1}, 0, 0x0}, /* (cAdd x <1>) */
        /* 254	*/ {{2,P2(P(9),S(194))         , cAdd        ,AnyParams       ,2}, 0, 0x0}, /* (cAdd & -% <2>) */
        /* 255	*/ {{1,P1(S(123))              , cAdd        ,AnyParams       ,1}, 0, 0x16}, /* (cAdd (cPow [(cLog [z]) -1]) <1>) */
        /* 256	*/ {{1,P1(S(334))              , cAdd        ,AnyParams       ,1}, 0, 0x0}, /* (cAdd (cMul -1 <2>) <1>) */
        /* 257	*/ {{1,P1(S(338))              , cAdd        ,AnyParams       ,2}, 0, 0x0}, /* (cAdd (cMul %@M <1>) <2>) */
        /* 258	*/ {{1,P1(S(357))              , cAdd        ,AnyParams       ,1}, 0, 0x16}, /* (cAdd (cMul (cPow [(cLog [z]) -1]) <2>) <1>) */
        /* 259	*/ {{2,P2(P(1),N(10))          , cAdd        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* ADD( % 1 ) */
        /* 260	*/ {{2,P2(P(9),S(194))         , cAdd        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* ADD( & -% ) */
        /* 261	*/ {{2,P2(N(2),P(14))          , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {-1 x}) */
        /* 262	*/ {{2,P2(N(2),P(21))          , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {-1 y}) */
        /* 263	*/ {{2,P2(N(2),S(13))          , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {-1 (cCeil [(cMul  <1>)])}) */
        /* 264	*/ {{2,P2(N(2),S(18))          , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {-1 (cCos [(cAdd {x (cMul {-1 y})})])}) */
        /* 265	*/ {{2,P2(N(2),S(19))          , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {-1 (cCos [(cAdd  <1>)])}) */
        /* 266	*/ {{2,P2(N(2),S(23))          , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {-1 (cCosh [x])}) */
        /* 267	*/ {{2,P2(N(2),S(30))          , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {-1 (cFloor [(cMul  <1>)])}) */
        /* 268	*/ {{2,P2(N(2),S(83))          , cMul        ,SelectedParams  ,0}, 0, 0x6}, /* (cMul {-1 (cPow [& x])}) */
        /* 269	*/ {{2,P2(N(2),S(118))         , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {-1 (cPow [(cCos [x]) 2])}) */
        /* 270	*/ {{2,P2(N(2),S(141))         , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {-1 (cPow [(cSin [x]) 2])}) */
        /* 271	*/ {{2,P2(N(2),S(153))         , cMul        ,SelectedParams  ,0}, 0, 0x6}, /* (cMul {-1 (cPow [/& x])}) */
        /* 272	*/ {{2,P2(N(2),S(165))         , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {-1 (cSin [(cMul  <1>)])}) */
        /* 273	*/ {{2,P2(N(2),S(170))         , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {-1 (cSinh [x])}) */
        /* 274	*/ {{2,P2(N(2),S(174))         , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {-1 (cSinh [(cMul  <1>)])}) */
        /* 275	*/ {{2,P2(N(2),S(182))         , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {-1 (cTan [(cMul  <1>)])}) */
        /* 276	*/ {{2,P2(N(2),S(190))         , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {-1 (cTanh [(cMul  <1>)])}) */
        /* 277	*/ {{2,P2(N(15),S(110))        , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {2.71828182846 (cPow [z (cAdd  <1>)])}) */
        /* 278	*/ {{3,P3(P(14),N(8),S(417))   , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {x 0.5 /%}) */
        /* 279	*/ {{2,P2(P(1),P(14))          , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {% x}) */
        /* 280	*/ {{2,P2(P(1),S(96))          , cMul        ,SelectedParams  ,0}, 0, 0x1}, /* (cMul {% (cPow [x@P z])}) */
        /* 281	*/ {{2,P2(P(1),S(212))         , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {% (cAdd {1 (cMul -1 <2>)})}) */
        /* 282	*/ {{2,P2(P(1),S(256))         , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {% (cAdd (cMul -1 <2>) <1>)}) */
        /* 283	*/ {{2,P2(P(14),P(21))         , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {x y}) */
        /* 284	*/ {{2,P2(P(14),S(104))        , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {x (cPow [y -1])}) */
        /* 285	*/ {{2,P2(P(14),S(391))        , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {x (cAnd  <1>)}) */
        /* 286	*/ {{3,P3(P(14),S(53),N(8))    , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {x LOG( % ) 0.5}) */
        /* 287	*/ {{2,P2(P(21),P(28))         , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {y z}) */
        /* 288	*/ {{2,P2(P(21),P(31))         , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {y a}) */
        /* 289	*/ {{2,P2(P(28),P(32))         , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {z b}) */
        /* 290	*/ {{3,P3(P(28),S(47),S(420))  , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {z (cLog [x]) /LOG( & )}) */
        /* 291	*/ {{2,P2(P(1),S(83))          , cMul        ,SelectedParams  ,0}, 0, 0x7}, /* (cMul {% (cPow [& x])}) */
        /* 292	*/ {{2,P2(P(1),S(153))         , cMul        ,SelectedParams  ,0}, 0, 0x7}, /* (cMul {% (cPow [/& x])}) */
        /* 293	*/ {{2,P2(P(9),S(243))         , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {& (cAdd  <2>)}) */
        /* 294	*/ {{2,P2(P(9),S(242))         , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {& (cAdd  <1>)}) */
        /* 295	*/ {{2,P2(P(14),S(54))         , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {x LOG( & )}) */
        /* 296	*/ {{2,P2(P(14),S(108))        , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {x (cPow [y %@N])}) */
        /* 297	*/ {{2,P2(S(23),N(2))          , cMul        ,SelectedParams  ,0}, 0, 0x4}, /* (cMul {(cCosh [x]) -1}) */
        /* 298	*/ {{2,P2(S(84),N(2))          , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {(cPow [& y]) -1}) */
        /* 299	*/ {{2,P2(S(11),N(13))         , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {(cAtanh [x]) 2}) */
        /* 300	*/ {{2,P2(S(105),P(1))         , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {(cPow [x y]) %}) */
        /* 301	*/ {{2,P2(S(194),P(14))        , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {-% x}) */
        /* 302	*/ {{2,P2(S(14),S(16))         , cMul        ,SelectedParams  ,0}, 0, 0x12}, /* (cMul {(cCos [x]) (cCos [y])}) */
        /* 303	*/ {{2,P2(S(47),S(417))        , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {(cLog [x]) /%}) */
        /* 304	*/ {{2,P2(S(65),N(2))          , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {(cPow [0.367879441171 x]) -1}) */
        /* 305	*/ {{2,P2(S(65),N(2))          , cMul        ,SelectedParams  ,0}, 0, 0x4}, /* (cMul {(cPow [0.367879441171 x]) -1}) */
        /* 306	*/ {{2,P2(S(67),N(2))          , cMul        ,SelectedParams  ,0}, 0, 0x4}, /* (cMul {(cPow [2.71828182846 x]) -1}) */
        /* 307	*/ {{2,P2(S(170),N(2))         , cMul        ,SelectedParams  ,0}, 0, 0x4}, /* (cMul {(cSinh [x]) -1}) */
        /* 308	*/ {{3,P3(S(25),N(13),P(1))    , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {(cCosh [(cMul {x LOG( & )})]) 2 %}) */
        /* 309	*/ {{2,P2(S(173),N(0))         , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {(cSinh [(cMul {x LOG( & )})]) -2}) */
        /* 310	*/ {{2,P2(S(24),N(13))         , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {(cCosh [(cLog [(cPow [& x])])]) 2}) */
        /* 311	*/ {{2,P2(S(171),N(13))        , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {(cSinh [(cLog [(cPow [& x])])]) 2}) */
        /* 312	*/ {{3,P3(S(173),N(13),P(1))   , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {(cSinh [(cMul {x LOG( & )})]) 2 %}) */
        /* 313	*/ {{2,P2(S(194),S(153))       , cMul        ,SelectedParams  ,0}, 0, 0x7}, /* (cMul {-% (cPow [/& x])}) */
        /* 314	*/ {{3,P3(S(14),S(16),N(2))    , cMul        ,SelectedParams  ,0}, 0, 0x12}, /* (cMul {(cCos [x]) (cCos [y]) -1}) */
        /* 315	*/ {{2,P2(S(14),S(160))        , cMul        ,SelectedParams  ,0}, 0, 0x12}, /* (cMul {(cCos [x]) (cSin [y])}) */
        /* 316	*/ {{2,P2(S(53),P(21))         , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {LOG( % ) y}) */
        /* 317	*/ {{2,P2(S(89),N(2))          , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {(cPow [x 2]) -1}) */
        /* 318	*/ {{2,P2(S(158),S(16))        , cMul        ,SelectedParams  ,0}, 0, 0x12}, /* (cMul {(cSin [x]) (cCos [y])}) */
        /* 319	*/ {{2,P2(S(110),S(73))        , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {(cPow [z (cAdd  <1>)]) (cPow [2.71828182846 (cMul  <2>)])}) */
        /* 320	*/ {{2,P2(S(155),S(106))       , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {POW( % & ) (cPow [% (cAdd  <1>)])}) */
        /* 321	*/ {{2,P2(S(155),S(107))       , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {POW( % & ) (cPow [x LOG( % )])}) */
        /* 322	*/ {{2,P2(S(155),S(152))       , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {POW( % & ) (cPow [(cMul  <1>) &])}) */
        /* 323	*/ {{2,P2(S(158),S(160))       , cMul        ,SelectedParams  ,0}, 0, 0x12}, /* (cMul {(cSin [x]) (cSin [y])}) */
        /* 324	*/ {{3,P3(S(14),S(160),N(2))   , cMul        ,SelectedParams  ,0}, 0, 0x12}, /* (cMul {(cCos [x]) (cSin [y]) -1}) */
        /* 325	*/ {{3,P3(S(158),S(160),N(2))  , cMul        ,SelectedParams  ,0}, 0, 0x12}, /* (cMul {(cSin [x]) (cSin [y]) -1}) */
        /* 326	*/ {{2,P2(S(97),S(240))        , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {(cPow [x MIN( % 1 )]) (cAdd {(cMul (cPow [x (cAdd {% -MIN( % 1 )})]) <1>) (cMul (cPow [x (cAdd {1 -MIN( % 1 )})]) <2>)})}) */
        /* 327	*/ {{2,P2(S(98),S(241))        , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {(cPow [x MIN( % & )]) (cAdd {(cMul (cPow [x (cAdd {% -MIN( % & )})]) <1>) (cMul (cPow [x (cAdd {& -MIN( % & )})]) <2>)})}) */
        /* 328	*/ {{2,P2(S(200),S(146))       , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {(cAdd {1 x})@D4 (cPow [(cAdd {1 (cMul {-1 x})}) -1])@D4}) */
        /* 329	*/ {{2,P2(S(375),P(21))        , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {(cNot [x]) y}) */
        /* 330	*/ {{2,P2(S(399),P(21))        , cMul        ,SelectedParams  ,0}, 0, 0x0}, /* (cMul {(cNotNot [x]) y}) */
        /* 331	*/ {{0,0                       , cMul        ,AnyParams       ,1}, 0, 0x0}, /* (cMul  <1>) */
        /* 332	*/ {{0,0                       , cMul        ,AnyParams       ,2}, 0, 0x0}, /* (cMul  <2>) */
        /* 333	*/ {{1,P1(N(2))                , cMul        ,AnyParams       ,1}, 0, 0x0}, /* (cMul -1 <1>) */
        /* 334	*/ {{1,P1(N(2))                , cMul        ,AnyParams       ,2}, 0, 0x0}, /* (cMul -1 <2>) */
        /* 335	*/ {{1,P1(P(0))                , cMul        ,AnyParams       ,1}, 0, 0x0}, /* (cMul %@N <1>) */
        /* 336	*/ {{1,P1(P(2))                , cMul        ,AnyParams       ,1}, 0, 0x0}, /* (cMul %@P <1>) */
        /* 337	*/ {{1,P1(P(2))                , cMul        ,AnyParams       ,1}, 0, 0x1}, /* (cMul %@P <1>) */
        /* 338	*/ {{1,P1(P(6))                , cMul        ,AnyParams       ,1}, 0, 0x0}, /* (cMul %@M <1>) */
        /* 339	*/ {{1,P1(P(1))                , cMul        ,AnyParams       ,1}, 0, 0x0}, /* (cMul % <1>) */
        /* 340	*/ {{1,P1(P(9))                , cMul        ,AnyParams       ,2}, 0, 0x0}, /* (cMul & <2>) */
        /* 341	*/ {{2,P2(P(1),P(9))           , cMul        ,AnyParams       ,1}, 0, 0x0}, /* (cMul % & <1>) */
        /* 342	*/ {{2,P2(P(9),S(417))         , cMul        ,AnyParams       ,2}, 0, 0x0}, /* (cMul & /% <2>) */
        /* 343	*/ {{1,P1(P(14))               , cMul        ,AnyParams       ,1}, 0, 0x4}, /* (cMul x <1>) */
        /* 344	*/ {{1,P1(P(14))               , cMul        ,AnyParams       ,2}, 0, 0x4}, /* (cMul x <2>) */
        /* 345	*/ {{1,P1(P(14))               , cMul        ,AnyParams       ,1}, 0, 0x0}, /* (cMul x <1>) */
        /* 346	*/ {{1,P1(P(14))               , cMul        ,AnyParams       ,2}, 0, 0x0}, /* (cMul x <2>) */
        /* 347	*/ {{1,P1(S(0))                , cMul        ,AnyParams       ,1}, 0, 0x0}, /* (cMul (cAbs [x]) <1>) */
        /* 348	*/ {{1,P1(S(47))               , cMul        ,AnyParams       ,1}, 0, 0x0}, /* (cMul (cLog [x]) <1>) */
        /* 349	*/ {{1,P1(S(53))               , cMul        ,AnyParams       ,1}, 0, 0x0}, /* (cMul LOG( % ) <1>) */
        /* 350	*/ {{1,P1(S(86))               , cMul        ,AnyParams       ,1}, 0, 0x4}, /* (cMul (cPow [x %@I@P]) <1>) */
        /* 351	*/ {{1,P1(S(87))               , cMul        ,AnyParams       ,2}, 0, 0x4}, /* (cMul (cPow [x &@I]) <2>) */
        /* 352	*/ {{1,P1(S(88))               , cMul        ,AnyParams       ,1}, 0, 0x0}, /* (cMul (cPow [% y]) <1>) */
        /* 353	*/ {{1,P1(S(99))               , cMul        ,AnyParams       ,2}, 0, 0x0}, /* (cMul (cPow [x (cAdd {1 -MIN( % 1 )})]) <2>) */
        /* 354	*/ {{1,P1(S(101))              , cMul        ,AnyParams       ,1}, 0, 0x0}, /* (cMul (cPow [x (cAdd {% -MIN( % 1 )})]) <1>) */
        /* 355	*/ {{1,P1(S(102))              , cMul        ,AnyParams       ,1}, 0, 0x0}, /* (cMul (cPow [x (cAdd {% -MIN( % & )})]) <1>) */
        /* 356	*/ {{1,P1(S(103))              , cMul        ,AnyParams       ,2}, 0, 0x0}, /* (cMul (cPow [x (cAdd {& -MIN( % & )})]) <2>) */
        /* 357	*/ {{1,P1(S(123))              , cMul        ,AnyParams       ,2}, 0, 0x0}, /* (cMul (cPow [(cLog [z]) -1]) <2>) */
        /* 358	*/ {{1,P1(S(194))              , cMul        ,AnyParams       ,1}, 0, 0x0}, /* (cMul -% <1>) */
        /* 359	*/ {{1,P1(S(194))              , cMul        ,AnyParams       ,2}, 0, 0x1}, /* (cMul -% <2>) */
        /* 360	*/ {{2,P2(S(123),S(47))        , cMul        ,AnyParams       ,1}, 0, 0x16}, /* (cMul (cPow [(cLog [z]) -1]) (cLog [x]) <1>) */
        /* 361	*/ {{2,P2(S(419),S(47))        , cMul        ,AnyParams       ,1}, 0, 0x1}, /* (cMul /LOG( % ) (cLog [x]) <1>) */
        /* 362	*/ {{2,P2(P(1),P(9))           , cMul        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* MUL( % & ) */
        /* 363	*/ {{2,P2(P(9),S(417))         , cMul        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* MUL( & /% ) */
        /* 364	*/ {{2,P2(S(54),S(419))        , cMul        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* MUL( LOG( & ) /LOG( % ) ) */
        /* 365	*/ {{2,P2(P(14),P(21))         , cEqual      ,PositionalParams,0}, 0, 0x12}, /* (cEqual [x y]) */
        /* 366	*/ {{2,P2(P(14),P(21))         , cEqual      ,PositionalParams,0}, 0, 0x0}, /* (cEqual [x y]) */
        /* 367	*/ {{2,P2(P(14),P(28))         , cEqual      ,PositionalParams,0}, 0, 0x20}, /* (cEqual [x z]) */
        /* 368	*/ {{2,P2(P(21),P(28))         , cEqual      ,PositionalParams,0}, 0, 0x24}, /* (cEqual [y z]) */
        /* 369	*/ {{2,P2(P(21),P(28))         , cEqual      ,PositionalParams,0}, 0, 0x0}, /* (cEqual [y z]) */
        /* 370	*/ {{2,P2(N(4),P(14))          , cLess       ,PositionalParams,0}, 0, 0x4}, /* (cLess [0 x]) */
        /* 371	*/ {{2,P2(P(21),P(14))         , cLess       ,PositionalParams,0}, 0, 0x0}, /* (cLess [y x]) */
        /* 372	*/ {{2,P2(P(14),P(21))         , cLess       ,PositionalParams,0}, 0, 0x12}, /* (cLess [x y]) */
        /* 373	*/ {{2,P2(P(14),P(21))         , cLessOrEq   ,PositionalParams,0}, 0, 0x0}, /* (cLessOrEq [x y]) */
        /* 374	*/ {{1,P1(P(1))                , cNot        ,PositionalParams,0}, 0, 0x0}, /* (cNot [%]) */
        /* 375	*/ {{1,P1(P(14))               , cNot        ,PositionalParams,0}, 0, 0x0}, /* (cNot [x]) */
        /* 376	*/ {{1,P1(P(21))               , cNot        ,PositionalParams,0}, 0, 0x0}, /* (cNot [y]) */
        /* 377	*/ {{1,P1(P(28))               , cNot        ,PositionalParams,0}, 0, 0x0}, /* (cNot [z]) */
        /* 378	*/ {{1,P1(S(278))              , cNot        ,PositionalParams,0}, 0, 0x0}, /* (cNot [(cMul {x 0.5 /%})]) */
        /* 379	*/ {{1,P1(S(385))              , cNot        ,PositionalParams,0}, 0, 0x0}, /* (cNot [(cAnd {x y})]) */
        /* 380	*/ {{1,P1(S(388))              , cNot        ,PositionalParams,0}, 0, 0x0}, /* (cNot [(cAnd {z (cIf [x y (cNot [%])])})]) */
        /* 381	*/ {{1,P1(S(389))              , cNot        ,PositionalParams,0}, 0, 0x0}, /* (cNot [(cAnd {z (cIf [x (cNot [%]) y])})]) */
        /* 382	*/ {{1,P1(S(392))              , cNot        ,PositionalParams,0}, 0, 0x0}, /* (cNot [(cOr {x y})]) */
        /* 383	*/ {{1,P1(S(395))              , cNot        ,PositionalParams,0}, 0, 0x0}, /* (cNot [(cOr {z (cIf [x y (cNot [%])])})]) */
        /* 384	*/ {{1,P1(S(396))              , cNot        ,PositionalParams,0}, 0, 0x0}, /* (cNot [(cOr {z (cIf [x (cNot [%]) y])})]) */
        /* 385	*/ {{2,P2(P(14),P(21))         , cAnd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAnd {x y}) */
        /* 386	*/ {{2,P2(P(21),P(31))         , cAnd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAnd {y a}) */
        /* 387	*/ {{2,P2(P(28),P(32))         , cAnd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAnd {z b}) */
        /* 388	*/ {{2,P2(P(28),S(35))         , cAnd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAnd {z (cIf [x y (cNot [%])])}) */
        /* 389	*/ {{2,P2(P(28),S(37))         , cAnd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAnd {z (cIf [x (cNot [%]) y])}) */
        /* 390	*/ {{2,P2(S(375),P(21))        , cAnd        ,SelectedParams  ,0}, 0, 0x0}, /* (cAnd {(cNot [x]) y}) */
        /* 391	*/ {{0,0                       , cAnd        ,AnyParams       ,1}, 0, 0x0}, /* (cAnd  <1>) */
        /* 392	*/ {{2,P2(P(14),P(21))         , cOr         ,SelectedParams  ,0}, 0, 0x0}, /* (cOr {x y}) */
        /* 393	*/ {{2,P2(P(21),P(31))         , cOr         ,SelectedParams  ,0}, 0, 0x0}, /* (cOr {y a}) */
        /* 394	*/ {{2,P2(P(28),P(32))         , cOr         ,SelectedParams  ,0}, 0, 0x0}, /* (cOr {z b}) */
        /* 395	*/ {{2,P2(P(28),S(35))         , cOr         ,SelectedParams  ,0}, 0, 0x0}, /* (cOr {z (cIf [x y (cNot [%])])}) */
        /* 396	*/ {{2,P2(P(28),S(37))         , cOr         ,SelectedParams  ,0}, 0, 0x0}, /* (cOr {z (cIf [x (cNot [%]) y])}) */
        /* 397	*/ {{2,P2(S(375),P(21))        , cOr         ,SelectedParams  ,0}, 0, 0x0}, /* (cOr {(cNot [x]) y}) */
        /* 398	*/ {{0,0                       , cOr         ,AnyParams       ,1}, 0, 0x0}, /* (cOr  <1>) */
        /* 399	*/ {{1,P1(P(14))               , cNotNot     ,PositionalParams,0}, 0, 0x0}, /* (cNotNot [x]) */
        /* 400	*/ {{1,P1(P(21))               , cNotNot     ,PositionalParams,0}, 0, 0x0}, /* (cNotNot [y]) */
        /* 401	*/ {{1,P1(S(215))              , cNotNot     ,PositionalParams,0}, 0, 0x0}, /* (cNotNot [(cAdd {x y})]) */
        /* 402	*/ {{1,P1(S(253))              , cNotNot     ,PositionalParams,0}, 0, 0x0}, /* (cNotNot [(cAdd x <1>)]) */
        /* 403	*/ {{1,P1(S(278))              , cNotNot     ,PositionalParams,0}, 0, 0x0}, /* (cNotNot [(cMul {x 0.5 /%})]) */
        /* 404	*/ {{1,P1(S(285))              , cNotNot     ,PositionalParams,0}, 0, 0x0}, /* (cNotNot [(cMul {x (cAnd  <1>)})]) */
        /* 405	*/ {{1,P1(S(331))              , cDeg        ,PositionalParams,0}, 0, 0x0}, /* (cDeg [(cMul  <1>)]) */
        /* 406	*/ {{1,P1(S(331))              , cRad        ,PositionalParams,0}, 0, 0x0}, /* (cRad [(cMul  <1>)]) */
        /* 407	*/ {{3,P3(P(14),P(21),S(391))  , cAbsAnd     ,SelectedParams  ,0}, 0, 0x0}, /* (cAbsAnd {x y (cAnd  <1>)}) */
        /* 408	*/ {{3,P3(P(14),P(21),S(398))  , cAbsOr      ,SelectedParams  ,0}, 0, 0x0}, /* (cAbsOr {x y (cOr  <1>)}) */
        /* 409	*/ {{1,P1(P(14))               , cAbsNot     ,PositionalParams,0}, 0, 0x0}, /* (cAbsNot [x]) */
        /* 410	*/ {{1,P1(P(14))               , cAbsNotNot  ,PositionalParams,0}, 0, 0x0}, /* (cAbsNotNot [x]) */
        /* 411	*/ {{1,P1(P(21))               , cAbsNotNot  ,PositionalParams,0}, 0, 0x0}, /* (cAbsNotNot [y]) */
        /* 412	*/ {{1,P1(P(28))               , cAbsNotNot  ,PositionalParams,0}, 0, 0x0}, /* (cAbsNotNot [z]) */
        /* 413	*/ {{3,P3(P(14),N(10),S(411))  , cAbsIf      ,PositionalParams,0}, 0, 0x0}, /* (cAbsIf [x 1 (cAbsNotNot [y])]) */
        /* 414	*/ {{3,P3(P(14),P(21),P(28))   , cAbsIf      ,PositionalParams,0}, 0, 0x0}, /* (cAbsIf [x y z]) */
        /* 415	*/ {{3,P3(P(14),S(331),N(4))   , cAbsIf      ,PositionalParams,0}, 0, 0x0}, /* (cAbsIf [x (cMul  <1>) 0]) */
        /* 416	*/ {{3,P3(P(14),S(411),N(4))   , cAbsIf      ,PositionalParams,0}, 0, 0x0}, /* (cAbsIf [x (cAbsNotNot [y]) 0]) */
        /* 417	*/ {{1,P1(P(1))                , cInv        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* /% */
        /* 418	*/ {{1,P1(P(9))                , cInv        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* /& */
        /* 419	*/ {{1,P1(S(53))               , cInv        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* /LOG( % ) */
        /* 420	*/ {{1,P1(S(54))               , cInv        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* /LOG( & ) */
        /* 421	*/ {{1,P1(S(64))               , cInv        ,GroupFunction   ,0}, Constness_Const, 0x0}, /* /MIN( % & ) */
        },

    };
}
namespace FPoptimizer_Grammar
{
    const Rule grammar_rules[241] =
    {
        /* 0:	@L (cAbs [x])
         *	->	x
         */		 {ProduceNewTree, true , 1,P1(P(14))               , {1,P1(P(14))               , cAbs        ,PositionalParams,0}},
        /* 1:	(cAtan [(cMul {x (cPow [y %@N])})])
         *	->	(cAtan2 [x (cPow [y -%])])
         */		 {ProduceNewTree, false, 1,P1(S(10))               , {1,P1(S(296))              , cAtan       ,PositionalParams,0}},
        /* 2:	(cAtan2 [(cPow [(cAdd {(cMul {(cPow [x 2]) -1}) 1}) 0.5])@D4 x@D4])
         *	->	(cAcos [x])
         */		 {ProduceNewTree, false, 1,P1(S(3))                , {2,P2(S(150),P(15))        , cAtan2      ,PositionalParams,0}},
        /* 3:	(cAtan2 [x@D4 (cPow [(cAdd {(cMul {(cPow [x 2]) -1}) 1}) 0.5])@D4])
         *	->	(cAsin [x])
         */		 {ProduceNewTree, false, 1,P1(S(5))                , {2,P2(P(15),S(150))        , cAtan2      ,PositionalParams,0}},
        /* 4:	(cAtan2 [(cMul x <1>)@D4 (cMul x <2>)@D4])
         *	:	(cMul  <1>) (cMul  <2>)
         */		 {ReplaceParams , false, 2,P2(S(331),S(332))       , {2,P2(S(343),S(344))       , cAtan2      ,PositionalParams,0}},
        /* 5:	(cCeil [x@I])
         *	->	x
         */		 {ProduceNewTree, false, 1,P1(P(14))               , {1,P1(P(16))               , cCeil       ,PositionalParams,0}},
        /* 6:	(cCeil [(cMul -1 <1>)])
         *	->	(cMul {-1 (cFloor [(cMul  <1>)])})
         */		 {ProduceNewTree, false, 1,P1(S(267))              , {1,P1(S(333))              , cCeil       ,PositionalParams,0}},
        /* 7:	(cCos [(cAdd {1.57079632679 (cMul %@N <1>)})])
         *	->	(cSin [(cMul -% <1>)])
         */		 {ProduceNewTree, false, 1,P1(S(168))              , {1,P1(S(214))              , cCos        ,PositionalParams,0}},
        /* 8:	(cCos [(cAdd -1.57079632679 <1>)])
         *	->	(cSin [(cAdd  <1>)])
         */		 {ProduceNewTree, false, 1,P1(S(163))              , {1,P1(S(245))              , cCos        ,PositionalParams,0}},
        /* 9:	(cCos [(cAcos [x])])
         *	->	x
         */		 {ProduceNewTree, false, 1,P1(P(14))               , {1,P1(S(3))                , cCos        ,PositionalParams,0}},
        /* 10:	(cCos [(cAbs [x])])
         *	:	x
         */		 {ReplaceParams , false, 1,P1(P(14))               , {1,P1(S(0))                , cCos        ,PositionalParams,0}},
        /* 11:	(cCos [(cMul -1 <1>)])
         *	:	(cMul  <1>)
         */		 {ReplaceParams , false, 1,P1(S(331))              , {1,P1(S(333))              , cCos        ,PositionalParams,0}},
        /* 12:	(cCosh [(cAsinh [x])])
         *	->	(cPow [(cAdd {(cPow [x 2]) 1}) 0.5])
         */		 {ProduceNewTree, false, 1,P1(S(149))              , {1,P1(S(6))                , cCosh       ,PositionalParams,0}},
        /* 13:	(cCosh [(cAbs [x])])
         *	:	x
         */		 {ReplaceParams , false, 1,P1(P(14))               , {1,P1(S(0))                , cCosh       ,PositionalParams,0}},
        /* 14:	(cCosh [(cMul -1 <1>)])
         *	:	(cMul  <1>)
         */		 {ReplaceParams , false, 1,P1(S(331))              , {1,P1(S(333))              , cCosh       ,PositionalParams,0}},
        /* 15:	(cFloor [x@I])
         *	->	x
         */		 {ProduceNewTree, false, 1,P1(P(14))               , {1,P1(P(16))               , cFloor      ,PositionalParams,0}},
        /* 16:	(cFloor [(cMul -1 <1>)])
         *	->	(cMul {-1 (cCeil [(cMul  <1>)])})
         */		 {ProduceNewTree, false, 1,P1(S(263))              , {1,P1(S(333))              , cFloor      ,PositionalParams,0}},
        /* 17:	(cFloor [(cAdd 0.5 <1>)])
         *	->	(cInt [(cAdd  <1>)])
         */		 {ProduceNewTree, false, 1,P1(S(46))               , {1,P1(S(246))              , cFloor      ,PositionalParams,0}},
        /* 18:	(cIf [x 0 y])
         *	->	(cMul {(cNot [x]) y})
         */		 {ProduceNewTree, false, 1,P1(S(329))              , {3,P3(P(14),N(4),P(21))    , cIf         ,PositionalParams,0}},
        /* 19:	(cIf [x@P y z])
         *	->	(cAbsIf [x y z])
         */		 {ProduceNewTree, false, 1,P1(S(414))              , {3,P3(P(17),P(21),P(28))   , cIf         ,PositionalParams,0}},
        /* 20:	(cIf [x 0 y@L])
         *	->	(cAnd {(cNot [x]) y})
         */		 {ProduceNewTree, false, 1,P1(S(390))              , {3,P3(P(14),N(4),P(22))    , cIf         ,PositionalParams,0}},
        /* 21:	(cIf [x 1 y@L])
         *	->	(cOr {x y})
         */		 {ProduceNewTree, false, 1,P1(S(392))              , {3,P3(P(14),N(10),P(22))   , cIf         ,PositionalParams,0}},
        /* 22:	(cIf [x y 0])
         *	->	(cMul {(cNotNot [x]) y})
         */		 {ProduceNewTree, false, 1,P1(S(330))              , {3,P3(P(14),P(21),N(4))    , cIf         ,PositionalParams,0}},
        /* 23:	(cIf [x y@L 0])
         *	->	(cAnd {x y})
         */		 {ProduceNewTree, false, 1,P1(S(385))              , {3,P3(P(14),P(22),N(4))    , cIf         ,PositionalParams,0}},
        /* 24:	(cIf [x y@L 1])
         *	->	(cOr {(cNot [x]) y})
         */		 {ProduceNewTree, false, 1,P1(S(397))              , {3,P3(P(14),P(22),N(10))   , cIf         ,PositionalParams,0}},
        /* 25:	(cIf [(cLess [x y])@D12 y@D8 x@D4])
         *	->	(cMax [x y])
         */		 {ProduceNewTree, false, 1,P1(S(57))               , {3,P3(S(372),P(23),P(15))  , cIf         ,PositionalParams,0}},
        /* 26:	(cIf [(cLess [x y])@D12 x@D4 y@D8])
         *	->	(cMin [x y])
         */		 {ProduceNewTree, false, 1,P1(S(60))               , {3,P3(S(372),P(15),P(23))  , cIf         ,PositionalParams,0}},
        /* 27:	(cIf [(cLess [0 x])@D4 (cFloor [x])@D4 (cCeil [x])@D4])
         *	->	(cTrunc [x])
         */		 {ProduceNewTree, false, 1,P1(S(193))              , {3,P3(S(370),S(29),S(12))  , cIf         ,PositionalParams,0}},
        /* 28:	(cIf [(cLessOrEq [x y]) z a])
         *	:	(cLess [y x]) a z
         */		 {ReplaceParams , false, 3,P3(S(371),P(31),P(28))  , {3,P3(S(373),P(28),P(31))  , cIf         ,PositionalParams,0}},
        /* 29:	@L (cIf [x (cAbsNotNot [y]) z])
         *	:	x y z
         */		 {ReplaceParams , true , 3,P3(P(14),P(21),P(28))   , {3,P3(P(14),S(411),P(28))  , cIf         ,PositionalParams,0}},
        /* 30:	@L (cIf [x y (cAbsNotNot [z])])
         *	:	x y z
         */		 {ReplaceParams , true , 3,P3(P(14),P(21),P(28))   , {3,P3(P(14),P(21),S(412))  , cIf         ,PositionalParams,0}},
        /* 31:	(cInt [x@I])
         *	->	x
         */		 {ProduceNewTree, false, 1,P1(P(14))               , {1,P1(P(16))               , cInt        ,PositionalParams,0}},
        /* 32:	(cLog [(cMul %@P <1>)])
         *	->	(cAdd {(cLog [(cMul  <1>)]) LOG( % )})
         */		 {ProduceNewTree, false, 1,P1(S(235))              , {1,P1(S(336))              , cLog        ,PositionalParams,0}},
        /* 33:	(cLog [(cMul (cPow [% y]) <1>)])
         *	->	(cAdd {(cMul {LOG( % ) y}) (cLog [(cMul  <1>)])})
         */		 {ProduceNewTree, false, 1,P1(S(234))              , {1,P1(S(352))              , cLog        ,PositionalParams,0}},
        /* 34:	(cLog [(cAdd {(cPow [(cAdd {-1 (cPow [x 2])}) 0.5])@D4 x@D4})])
         *	->	(cAcosh [x])
         */		 {ProduceNewTree, false, 1,P1(S(4))                , {1,P1(S(230))              , cLog        ,PositionalParams,0}},
        /* 35:	(cLog [(cMul {(cAdd {1 x})@D4 (cPow [(cAdd {1 (cMul {-1 x})}) -1])@D4})])
         *	->	(cMul {(cAtanh [x]) 2})
         */		 {ProduceNewTree, false, 1,P1(S(299))              , {1,P1(S(328))              , cLog        ,PositionalParams,0}},
        /* 36:	(cMax x@D4 x@D4)
         *	:	x
         */		 {ReplaceParams , false, 1,P1(P(14))               , {2,P2(P(15),P(15))         , cMax        ,AnyParams       ,0}},
        /* 37:	(cMax (cIf [x y z])@D4 (cIf [x a b])@D4)
         *	:	(cIf [x (cMax [y a]) (cMax [z b])])
         */		 {ReplaceParams , false, 1,P1(S(39))               , {2,P2(S(32),S(36))         , cMax        ,AnyParams       ,0}},
        /* 38:	(cMin x@D4 x@D4)
         *	:	x
         */		 {ReplaceParams , false, 1,P1(P(14))               , {2,P2(P(15),P(15))         , cMin        ,AnyParams       ,0}},
        /* 39:	(cMin (cIf [x y z])@D4 (cIf [x a b])@D4)
         *	:	(cIf [x (cMin [y a]) (cMin [z b])])
         */		 {ReplaceParams , false, 1,P1(S(40))               , {2,P2(S(32),S(36))         , cMin        ,AnyParams       ,0}},
        /* 40:	(cPow [(cMul %@P <1>) &])
         *	->	(cMul {POW( % & ) (cPow [(cMul  <1>) &])})
         */		 {ProduceNewTree, false, 1,P1(S(322))              , {2,P2(S(336),P(9))         , cPow        ,PositionalParams,0}},
        /* 41:	(cPow [% (cAdd {(cLog [x]) &})])
         *	->	(cMul {POW( % & ) (cPow [x LOG( % )])})
         */		 {ProduceNewTree, false, 1,P1(S(321))              , {2,P2(P(1),S(228))         , cPow        ,PositionalParams,0}},
        /* 42:	(cPow [(cMul %@N <1>) &@E])
         *	->	(cMul {POW( % & ) (cPow [(cMul  <1>) &])})
         */		 {ProduceNewTree, false, 1,P1(S(322))              , {2,P2(S(335),P(10))        , cPow        ,PositionalParams,0}},
        /* 43:	(cPow [z@D16 (cAdd (cMul (cPow [(cLog [z]) -1]) <2>) <1>)@D16])
         *	->	(cMul {(cPow [z (cAdd  <1>)]) (cPow [2.71828182846 (cMul  <2>)])})
         */		 {ProduceNewTree, false, 1,P1(S(319))              , {2,P2(P(29),S(258))        , cPow        ,PositionalParams,0}},
        /* 44:	(cPow [z@D16 (cAdd (cPow [(cLog [z]) -1]) <1>)@D16])
         *	->	(cMul {2.71828182846 (cPow [z (cAdd  <1>)])})
         */		 {ProduceNewTree, false, 1,P1(S(277))              , {2,P2(P(29),S(255))        , cPow        ,PositionalParams,0}},
        /* 45:	(cPow [% (cAdd &@M <1>)])
         *	->	(cMul {POW( % & ) (cPow [% (cAdd  <1>)])})
         */		 {ProduceNewTree, false, 1,P1(S(320))              , {2,P2(P(1),S(249))         , cPow        ,PositionalParams,0}},
        /* 46:	(cPow [(cPow [x y@O]) z])
         *	:	x (cMul {y z})
         */		 {ReplaceParams , false, 2,P2(P(14),S(287))        , {2,P2(S(90),P(28))         , cPow        ,PositionalParams,0}},
        /* 47:	(cPow [(cPow [x y@F]) z])
         *	:	x (cMul {y z})
         */		 {ReplaceParams , false, 2,P2(P(14),S(287))        , {2,P2(S(91),P(28))         , cPow        ,PositionalParams,0}},
        /* 48:	(cPow [(cPow [x@P y]) z])
         *	:	x (cMul {y z})
         */		 {ReplaceParams , false, 2,P2(P(14),S(287))        , {2,P2(S(92),P(28))         , cPow        ,PositionalParams,0}},
        /* 49:	(cPow [(cPow [x y])@P z])
         *	:	(cAbs [x]) (cMul {y z})
         */		 {ReplaceParams , false, 2,P2(S(0),S(287))         , {2,P2(S(93),P(28))         , cPow        ,PositionalParams,0}},
        /* 50:	(cPow [(cSin [x]) %@N])
         *	:	(cCsc [x]) -%
         */		 {ReplaceParams , false, 2,P2(S(28),S(194))        , {2,P2(S(158),P(0))         , cPow        ,PositionalParams,0}},
        /* 51:	(cPow [(cCos [x]) %@N])
         *	:	(cSec [x]) -%
         */		 {ReplaceParams , false, 2,P2(S(157),S(194))       , {2,P2(S(14),P(0))          , cPow        ,PositionalParams,0}},
        /* 52:	(cPow [(cTan [x]) %@N])
         *	:	(cCot [x]) -%
         */		 {ReplaceParams , false, 2,P2(S(27),S(194))        , {2,P2(S(177),P(0))         , cPow        ,PositionalParams,0}},
        /* 53:	(cPow [% (cLog [x])])
         *	:	x LOG( % )
         */		 {ReplaceParams , false, 2,P2(P(14),S(53))         , {2,P2(P(1),S(47))          , cPow        ,PositionalParams,0}},
        /* 54:	(cPow [% (cMul (cLog [x]) <1>)])
         *	:	x (cMul LOG( % ) <1>)
         */		 {ReplaceParams , false, 2,P2(P(14),S(349))        , {2,P2(P(1),S(348))         , cPow        ,PositionalParams,0}},
        /* 55:	(cPow [z@D16 (cMul (cPow [(cLog [z]) -1]) (cLog [x]) <1>)@D16])
         *	:	x (cMul  <1>)
         */		 {ReplaceParams , false, 2,P2(P(14),S(331))        , {2,P2(P(29),S(360))        , cPow        ,PositionalParams,0}},
        /* 56:	(cPow [%@D1 (cMul /LOG( % ) (cLog [x]) <1>)@D1])
         *	:	x (cMul  <1>)
         */		 {ReplaceParams , false, 2,P2(P(14),S(331))        , {2,P2(P(3),S(361))         , cPow        ,PositionalParams,0}},
        /* 57:	(cPow [(cPow [x y]) z@I])
         *	:	x (cMul {y z})
         */		 {ReplaceParams , false, 2,P2(P(14),S(287))        , {2,P2(S(105),P(30))        , cPow        ,PositionalParams,0}},
        /* 58:	(cPow [(cAbs [x]) y@E])
         *	:	x y
         */		 {ReplaceParams , false, 2,P2(P(14),P(21))         , {2,P2(S(0),P(26))          , cPow        ,PositionalParams,0}},
        /* 59:	(cPow [(cMul (cAbs [x]) <1>) y@E])
         *	:	(cMul x <1>) y
         */		 {ReplaceParams , false, 2,P2(S(345),P(21))        , {2,P2(S(347),P(26))        , cPow        ,PositionalParams,0}},
        /* 60:	(cSin [(cMul -1 <1>)])
         *	->	(cMul {-1 (cSin [(cMul  <1>)])})
         */		 {ProduceNewTree, false, 1,P1(S(272))              , {1,P1(S(333))              , cSin        ,PositionalParams,0}},
        /* 61:	(cSin [(cAdd {1.57079632679 (cMul %@N <1>)})])
         *	->	(cCos [(cMul -% <1>)])
         */		 {ProduceNewTree, false, 1,P1(S(21))               , {1,P1(S(214))              , cSin        ,PositionalParams,0}},
        /* 62:	(cSin [(cAdd -1.57079632679 <1>)])
         *	->	(cMul {-1 (cCos [(cAdd  <1>)])})
         */		 {ProduceNewTree, false, 1,P1(S(265))              , {1,P1(S(245))              , cSin        ,PositionalParams,0}},
        /* 63:	(cSin [(cAsin [x])])
         *	->	x
         */		 {ProduceNewTree, false, 1,P1(P(14))               , {1,P1(S(5))                , cSin        ,PositionalParams,0}},
        /* 64:	(cSinh [(cMul -1 <1>)])
         *	->	(cMul {-1 (cSinh [(cMul  <1>)])})
         */		 {ProduceNewTree, false, 1,P1(S(274))              , {1,P1(S(333))              , cSinh       ,PositionalParams,0}},
        /* 65:	(cSinh [(cAcosh [x])])
         *	->	(cPow [(cAdd {(cPow [x 2]) -1}) 0.5])
         */		 {ProduceNewTree, false, 1,P1(S(148))              , {1,P1(S(4))                , cSinh       ,PositionalParams,0}},
        /* 66:	(cTan [(cMul -1 <1>)])
         *	->	(cMul {-1 (cTan [(cMul  <1>)])})
         */		 {ProduceNewTree, false, 1,P1(S(275))              , {1,P1(S(333))              , cTan        ,PositionalParams,0}},
        /* 67:	(cTan [(cAtan [x])])
         *	->	x
         */		 {ProduceNewTree, false, 1,P1(P(14))               , {1,P1(S(8))                , cTan        ,PositionalParams,0}},
        /* 68:	(cTan [(cAtan2 [x y])])
         *	->	(cMul {x (cPow [y -1])})
         */		 {ProduceNewTree, false, 1,P1(S(284))              , {1,P1(S(9))                , cTan        ,PositionalParams,0}},
        /* 69:	(cTanh [(cMul -1 <1>)])
         *	->	(cMul {-1 (cTanh [(cMul  <1>)])})
         */		 {ProduceNewTree, false, 1,P1(S(276))              , {1,P1(S(333))              , cTanh       ,PositionalParams,0}},
        /* 70:	(cTrunc [x@I])
         *	->	x
         */		 {ProduceNewTree, false, 1,P1(P(14))               , {1,P1(P(16))               , cTrunc      ,PositionalParams,0}},
        /* 71:	(cAdd (cPow [(cAdd {1 (cPow [(cAdd  <1>) 2])}) 0.5]) <1>)
         *	->	(cPow [2.71828182846 (cAsinh [(cAdd  <1>)])])
         */		 {ProduceNewTree, false, 1,P1(S(72))               , {1,P1(S(145))              , cAdd        ,AnyParams       ,1}},
        /* 72:	(cAdd (cPow [(cAdd {1 (cPow [(cAdd  <1>) 2])}) -0.5]) <1>)
         *	->	(cPow [0.367879441171 (cAsinh [(cAdd  <1>)])])
         */		 {ProduceNewTree, false, 1,P1(S(71))               , {1,P1(S(144))              , cAdd        ,AnyParams       ,1}},
        /* 73:	(cAdd (cPow [(cAdd {1 (cPow [x 2])}) 0.5])@D4 x@D4)
         *	:	(cPow [2.71828182846 (cAsinh [x])])
         */		 {ReplaceParams , false, 1,P1(S(70))               , {2,P2(S(143),P(15))        , cAdd        ,AnyParams       ,0}},
        /* 74:	(cAdd (cPow [(cAdd {1 (cPow [x 2])}) -0.5])@D4 x@D4)
         *	:	(cPow [0.367879441171 (cAsinh [x])])
         */		 {ReplaceParams , false, 1,P1(S(69))               , {2,P2(S(142),P(15))        , cAdd        ,AnyParams       ,0}},
        /* 75:	(cAdd (cIf [x y z])@D4 (cIf [x a b])@D4)
         *	:	(cIf [x (cAdd {y a}) (cAdd {z b})])
         */		 {ReplaceParams , false, 1,P1(S(41))               , {2,P2(S(32),S(36))         , cAdd        ,AnyParams       ,0}},
        /* 76:	(cAdd (cLog [x]) (cLog [y]))
         *	:	(cLog [(cMul {x y})])
         */		 {ReplaceParams , false, 1,P1(S(51))               , {2,P2(S(47),S(48))         , cAdd        ,AnyParams       ,0}},
        /* 77:	(cAdd (cMul (cPow [x %@I@P]) <1>)@D4 (cMul (cPow [x &@I]) <2>)@D4)
         *	:	(cMul {(cPow [x MIN( % & )]) (cAdd {(cMul (cPow [x (cAdd {% -MIN( % & )})]) <1>) (cMul (cPow [x (cAdd {& -MIN( % & )})]) <2>)})})
         */		 {ReplaceParams , false, 1,P1(S(327))              , {2,P2(S(350),S(351))       , cAdd        ,AnyParams       ,0}},
        /* 78:	(cAdd (cMul (cPow [x %@I@P]) <1>)@D4 (cMul x <2>)@D4)
         *	:	(cMul {(cPow [x MIN( % 1 )]) (cAdd {(cMul (cPow [x (cAdd {% -MIN( % 1 )})]) <1>) (cMul (cPow [x (cAdd {1 -MIN( % 1 )})]) <2>)})})
         */		 {ReplaceParams , false, 1,P1(S(326))              , {2,P2(S(350),S(344))       , cAdd        ,AnyParams       ,0}},
        /* 79:	(cAdd (cMul %@P <1>)@D1 (cMul -% <2>)@D1)
         *	:	(cMul {% (cAdd (cMul -1 <2>) <1>)})
         */		 {ReplaceParams , false, 1,P1(S(282))              , {2,P2(S(337),S(359))       , cAdd        ,AnyParams       ,0}},
        /* 80:	(cAdd %@M@D1 (cMul -% <2>)@D1)
         *	:	(cMul {% (cAdd {1 (cMul -1 <2>)})})
         */		 {ReplaceParams , false, 1,P1(S(281))              , {2,P2(P(5),S(359))         , cAdd        ,AnyParams       ,0}},
        /* 81:	(cAdd (cPow [(cSin [x]) 2])@D4 (cPow [(cCos [x]) 2])@D4)
         *	:	1
         */		 {ReplaceParams , false, 1,P1(N(10))               , {2,P2(S(140),S(117))       , cAdd        ,AnyParams       ,0}},
        /* 82:	(cAdd 1 (cMul {-1 (cPow [(cSin [x]) 2])}))
         *	:	(cPow [(cCos [x]) 2])
         */		 {ReplaceParams , false, 1,P1(S(118))              , {2,P2(N(10),S(270))        , cAdd        ,AnyParams       ,0}},
        /* 83:	(cAdd 1 (cMul {-1 (cPow [(cCos [x]) 2])}))
         *	:	(cPow [(cSin [x]) 2])
         */		 {ReplaceParams , false, 1,P1(S(141))              , {2,P2(N(10),S(269))        , cAdd        ,AnyParams       ,0}},
        /* 84:	(cAdd (cMul {(cSin [x]) (cCos [y])})@D12 (cMul {(cCos [x]) (cSin [y])})@D12)
         *	:	(cSin [(cAdd {x y})])
         */		 {ReplaceParams , false, 1,P1(S(161))              , {2,P2(S(318),S(315))       , cAdd        ,AnyParams       ,0}},
        /* 85:	(cAdd (cMul {(cSin [x]) (cCos [y])})@D12 (cMul {(cCos [x]) (cSin [y]) -1})@D12)
         *	:	(cSin [(cAdd {x (cMul {-1 y})})])
         */		 {ReplaceParams , false, 1,P1(S(162))              , {2,P2(S(318),S(324))       , cAdd        ,AnyParams       ,0}},
        /* 86:	(cAdd (cMul {(cCos [x]) (cCos [y])})@D12 (cMul {(cSin [x]) (cSin [y])})@D12)
         *	:	(cCos [(cAdd {x y})])
         */		 {ReplaceParams , false, 1,P1(S(17))               , {2,P2(S(302),S(323))       , cAdd        ,AnyParams       ,0}},
        /* 87:	(cAdd (cMul {(cCos [x]) (cCos [y]) -1})@D12 (cMul {(cSin [x]) (cSin [y])})@D12)
         *	:	(cMul {-1 (cCos [(cAdd {x (cMul {-1 y})})])})
         */		 {ReplaceParams , false, 1,P1(S(264))              , {2,P2(S(314),S(323))       , cAdd        ,AnyParams       ,0}},
        /* 88:	(cAdd (cMul {(cCos [x]) (cCos [y])})@D12 (cMul {(cSin [x]) (cSin [y]) -1})@D12)
         *	:	(cCos [(cAdd {x (cMul {-1 y})})])
         */		 {ReplaceParams , false, 1,P1(S(18))               , {2,P2(S(302),S(325))       , cAdd        ,AnyParams       ,0}},
        /* 89:	(cAdd (cPow [& x])@D6 (cMul {-1 (cPow [/& x])})@D6)
         *	:	(cMul {(cSinh [(cLog [(cPow [& x])])]) 2})
         */		 {ReplaceParams , false, 1,P1(S(311))              , {2,P2(S(82),S(271))        , cAdd        ,AnyParams       ,0}},
        /* 90:	(cAdd (cPow [& x])@D6 (cPow [/& x])@D6)
         *	:	(cMul {(cCosh [(cLog [(cPow [& x])])]) 2})
         */		 {ReplaceParams , false, 1,P1(S(310))              , {2,P2(S(82),S(154))        , cAdd        ,AnyParams       ,0}},
        /* 91:	(cAdd (cMul {-1 (cPow [& x])})@D6 (cPow [/& x])@D6)
         *	:	(cMul {(cSinh [(cMul {x LOG( & )})]) -2})
         */		 {ReplaceParams , false, 1,P1(S(309))              , {2,P2(S(268),S(154))       , cAdd        ,AnyParams       ,0}},
        /* 92:	(cAdd (cMul {% (cPow [& x])})@D7 (cMul {-% (cPow [/& x])})@D7)
         *	:	(cMul {(cSinh [(cMul {x LOG( & )})]) 2 %})
         */		 {ReplaceParams , false, 1,P1(S(312))              , {2,P2(S(291),S(313))       , cAdd        ,AnyParams       ,0}},
        /* 93:	(cAdd (cMul {% (cPow [& x])})@D7 (cMul {% (cPow [/& x])})@D7)
         *	:	(cMul {(cCosh [(cMul {x LOG( & )})]) 2 %})
         */		 {ReplaceParams , false, 1,P1(S(308))              , {2,P2(S(291),S(292))       , cAdd        ,AnyParams       ,0}},
        /* 94:	(cAdd (cCosh [x])@D4 (cSinh [x])@D4)
         *	:	(cPow [2.71828182846 x])
         */		 {ReplaceParams , false, 1,P1(S(67))               , {2,P2(S(22),S(169))        , cAdd        ,AnyParams       ,0}},
        /* 95:	(cAdd (cMul {(cCosh [x]) -1})@D4 (cSinh [x])@D4)
         *	:	(cMul {(cPow [0.367879441171 x]) -1})
         */		 {ReplaceParams , false, 1,P1(S(304))              , {2,P2(S(297),S(169))       , cAdd        ,AnyParams       ,0}},
        /* 96:	(cAdd (cCosh [x])@D4 (cMul {(cPow [2.71828182846 x]) -1})@D4)
         *	:	(cMul {-1 (cSinh [x])})
         */		 {ReplaceParams , false, 1,P1(S(273))              , {2,P2(S(22),S(306))        , cAdd        ,AnyParams       ,0}},
        /* 97:	(cAdd (cSinh [x])@D4 (cMul {(cPow [2.71828182846 x]) -1})@D4)
         *	:	(cMul {-1 (cCosh [x])})
         */		 {ReplaceParams , false, 1,P1(S(266))              , {2,P2(S(169),S(306))       , cAdd        ,AnyParams       ,0}},
        /* 98:	(cAdd (cCosh [x])@D4 (cMul {(cSinh [x]) -1})@D4)
         *	:	(cPow [0.367879441171 x])
         */		 {ReplaceParams , false, 1,P1(S(65))               , {2,P2(S(22),S(307))        , cAdd        ,AnyParams       ,0}},
        /* 99:	(cAdd (cMul {(cSinh [x]) -1})@D4 (cPow [2.71828182846 x])@D4)
         *	:	(cCosh [x])
         */		 {ReplaceParams , false, 1,P1(S(23))               , {2,P2(S(307),S(68))        , cAdd        ,AnyParams       ,0}},
        /* 100:	(cAdd (cMul {(cCosh [x]) -1})@D4 (cPow [2.71828182846 x])@D4)
         *	:	(cSinh [x])
         */		 {ReplaceParams , false, 1,P1(S(170))              , {2,P2(S(297),S(68))        , cAdd        ,AnyParams       ,0}},
        /* 101:	(cAdd (cCosh [x])@D4 (cMul {(cPow [0.367879441171 x]) -1})@D4)
         *	:	(cSinh [x])
         */		 {ReplaceParams , false, 1,P1(S(170))              , {2,P2(S(22),S(305))        , cAdd        ,AnyParams       ,0}},
        /* 102:	(cAdd (cSinh [x])@D4 (cPow [0.367879441171 x])@D4)
         *	:	(cCosh [x])
         */		 {ReplaceParams , false, 1,P1(S(23))               , {2,P2(S(169),S(66))        , cAdd        ,AnyParams       ,0}},
        /* 103:	(cAdd (cMul {(cCosh [x]) -1})@D4 (cPow [0.367879441171 x])@D4)
         *	:	(cMul {-1 (cSinh [x])})
         */		 {ReplaceParams , false, 1,P1(S(273))              , {2,P2(S(297),S(66))        , cAdd        ,AnyParams       ,0}},
        /* 104:	(cMul {%@D1 (cAdd {1 (cMul {(cLog [x]) /%})})@D1})
         *	->	(cAdd {(cLog [x]) %})
         */		 {ProduceNewTree, false, 1,P1(S(225))              , {2,P2(P(3),S(211))         , cMul        ,SelectedParams  ,0}},
        /* 105:	(cMul x@L <1>)
         *	->	(cAbsIf [x (cMul  <1>) 0])
         */		 {ProduceNewTree, false, 1,P1(S(415))              , {1,P1(P(19))               , cMul        ,AnyParams       ,1}},
        /* 106:	(cMul 57.2957795131 <1>)
         *	->	(cDeg [(cMul  <1>)])
         */		 {ProduceNewTree, false, 1,P1(S(405))              , {1,P1(N(16))               , cMul        ,AnyParams       ,1}},
        /* 107:	(cMul 0.0174532925199 <1>)
         *	->	(cRad [(cMul  <1>)])
         */		 {ProduceNewTree, false, 1,P1(S(406))              , {1,P1(N(5))                , cMul        ,AnyParams       ,1}},
        /* 108:	(cMul (cPow [(cMul x <2>) -1])@D4 x@D4)
         *	:	(cPow [(cMul  <2>) -1])
         */		 {ReplaceParams , false, 1,P1(S(138))              , {2,P2(S(139),P(15))        , cMul        ,AnyParams       ,0}},
        /* 109:	(cMul (cIf [x y z])@D4 (cIf [x a b])@D4)
         *	:	(cIf [x (cMul {y a}) (cMul {z b})])
         */		 {ReplaceParams , false, 1,P1(S(42))               , {2,P2(S(32),S(36))         , cMul        ,AnyParams       ,0}},
        /* 110:	(cMul (cNot [x]) (cNot [y]))
         *	:	(cNot [(cOr {x y})])
         */		 {ReplaceParams , false, 1,P1(S(382))              , {2,P2(S(375),S(376))       , cMul        ,AnyParams       ,0}},
        /* 111:	(cMul (cNotNot [x]) (cNotNot [y]))
         *	:	(cAnd {x y})
         */		 {ReplaceParams , false, 1,P1(S(385))              , {2,P2(S(399),S(400))       , cMul        ,AnyParams       ,0}},
        /* 112:	(cMul (cNot [x]) (cNotNot [y]))
         *	:	(cAnd {(cNot [x]) y})
         */		 {ReplaceParams , false, 1,P1(S(390))              , {2,P2(S(375),S(400))       , cMul        ,AnyParams       ,0}},
        /* 113:	(cMul (cAbs [x]) (cAbs [y]))
         *	:	(cAbs [(cMul {x y})])
         */		 {ReplaceParams , false, 1,P1(S(2))                , {2,P2(S(0),S(1))           , cMul        ,AnyParams       ,0}},
        /* 114:	(cMul (cAdd (cMul %@M <1>) <2>) &)
         *	:	(cAdd {(cMul % & <1>) (cMul {& (cAdd  <2>)})})
         */		 {ReplaceParams , false, 1,P1(S(238))              , {2,P2(S(257),P(9))         , cMul        ,AnyParams       ,0}},
        /* 115:	(cMul (cAdd %@M <1>) &)
         *	:	(cAdd {MUL( % & ) (cMul {& (cAdd  <1>)})})
         */		 {ReplaceParams , false, 1,P1(S(239))              , {2,P2(S(247),P(9))         , cMul        ,AnyParams       ,0}},
        /* 116:	(cMul (cPow [x y])@D4 (cAdd {%@1 (cPow [x z])})@D4)
         *	:	(cAdd {(cMul {(cPow [x y]) %}) (cPow [x (cAdd {y z})])})
         */		 {ReplaceParams , false, 1,P1(S(233))              , {2,P2(S(94),S(217))        , cMul        ,AnyParams       ,0}},
        /* 117:	(cMul (cPow [& y]) (cAdd {1 (cPow [x@P z])}))
         *	:	(cAdd {(cPow [& y]) (cPow [& (cAdd {y (cMul {z (cLog [x]) /LOG( & )})})])})
         */		 {ReplaceParams , false, 1,P1(S(229))              , {2,P2(S(84),S(206))        , cMul        ,AnyParams       ,0}},
        /* 118:	(cMul (cPow [& y]) (cAdd {-1 (cPow [x@P z])}))
         *	:	(cAdd {(cMul {(cPow [& y]) -1}) (cPow [& (cAdd {y (cMul {z (cLog [x]) /LOG( & )})})])})
         */		 {ReplaceParams , false, 1,P1(S(232))              , {2,P2(S(84),S(204))        , cMul        ,AnyParams       ,0}},
        /* 119:	(cMul -1 (cSin [(cMul %@N <1>)]))
         *	:	(cSin [(cMul -% <1>)])
         */		 {ReplaceParams , false, 1,P1(S(168))              , {2,P2(N(2),S(166))         , cMul        ,AnyParams       ,0}},
        /* 120:	(cMul -1 (cSinh [(cMul %@N <1>)]))
         *	:	(cSinh [(cMul -% <1>)])
         */		 {ReplaceParams , false, 1,P1(S(176))              , {2,P2(N(2),S(175))         , cMul        ,AnyParams       ,0}},
        /* 121:	(cMul (cPow [(cSinh [x]) -1])@D4 (cCosh [x])@D4)
         *	:	(cPow [(cTanh [x]) -1])
         */		 {ReplaceParams , false, 1,P1(S(133))              , {2,P2(S(130),S(22))        , cMul        ,AnyParams       ,0}},
        /* 122:	(cMul (cTanh [x])@D4 (cCosh [x])@D4)
         *	:	(cSinh [x])
         */		 {ReplaceParams , false, 1,P1(S(170))              , {2,P2(S(187),S(22))        , cMul        ,AnyParams       ,0}},
        /* 123:	(cMul (cPow [(cTanh [x]) -1])@D4 (cSinh [x])@D4)
         *	:	(cCosh [x])
         */		 {ReplaceParams , false, 1,P1(S(23))               , {2,P2(S(134),S(169))       , cMul        ,AnyParams       ,0}},
        /* 124:	(cMul (cPow [(cTan [x]) -1])@D4 (cSin [x])@D4)
         *	:	(cCos [x])
         */		 {ReplaceParams , false, 1,P1(S(14))               , {2,P2(S(131),S(159))       , cMul        ,AnyParams       ,0}},
        /* 125:	(cMul (cSin [x])@D4 (cPow [(cCos [x]) -1])@D4)
         *	:	(cTan [x])
         */		 {ReplaceParams , false, 1,P1(S(177))              , {2,P2(S(159),S(111))       , cMul        ,AnyParams       ,0}},
        /* 126:	(cMul (cTan [x])@D4 (cPow [(cSin [x]) -1])@D4)
         *	:	(cPow [(cCos [x]) -1])
         */		 {ReplaceParams , false, 1,P1(S(112))              , {2,P2(S(178),S(129))       , cMul        ,AnyParams       ,0}},
        /* 127:	(cMul (cPow [(cSin [x]) -1])@D4 (cCos [x])@D4)
         *	:	(cPow [(cTan [x]) -1])
         */		 {ReplaceParams , false, 1,P1(S(132))              , {2,P2(S(129),S(15))        , cMul        ,AnyParams       ,0}},
        /* 128:	(cMul (cTan [x])@D4 (cCos [x])@D4)
         *	:	(cSin [x])
         */		 {ReplaceParams , false, 1,P1(S(158))              , {2,P2(S(178),S(15))        , cMul        ,AnyParams       ,0}},
        /* 129:	(cMul (cTan [(cAdd {1.57079632679 (cMul {-1 x})})])@D4 (cTan [x])@D4)
         *	:	1
         */		 {ReplaceParams , false, 1,P1(N(10))               , {2,P2(S(179),S(178))       , cMul        ,AnyParams       ,0}},
        /* 130:	(cMul (cSin [(cMul % <1>)])@D1 (cPow [(cCos [(cMul -% <1>)]) -1])@D1)
         *	:	(cTan [(cMul % <1>)])
         */		 {ReplaceParams , false, 1,P1(S(183))              , {2,P2(S(167),S(114))       , cMul        ,AnyParams       ,0}},
        /* 131:	(cMul (cTan [(cAdd {1.57079632679 (cMul -1 <1>)})]) (cTan [(cMul  <1>)]))
         *	:	1
         */		 {ReplaceParams , false, 1,P1(N(10))               , {2,P2(S(180),S(182))       , cMul        ,AnyParams       ,0}},
        /* 132:	(cMul -1 (cTan [(cMul %@N <1>)]))
         *	:	(cTan [(cMul -% <1>)])
         */		 {ReplaceParams , false, 1,P1(S(185))              , {2,P2(N(2),S(184))         , cMul        ,AnyParams       ,0}},
        /* 133:	(cMul (cSinh [x])@D4 (cPow [(cCosh [x]) -1])@D4)
         *	:	(cTanh [x])
         */		 {ReplaceParams , false, 1,P1(S(186))              , {2,P2(S(169),S(115))       , cMul        ,AnyParams       ,0}},
        /* 134:	(cMul (cTanh [x])@D4 (cPow [(cSinh [x]) -1])@D4)
         *	:	(cPow [(cCosh [x]) -1])
         */		 {ReplaceParams , false, 1,P1(S(116))              , {2,P2(S(187),S(130))       , cMul        ,AnyParams       ,0}},
        /* 135:	(cMul (cSinh [(cMul {% x})])@D5 (cPow [(cCosh [(cMul {-% x})]) -1])@D5)
         *	:	(cTanh [(cMul {% x})])
         */		 {ReplaceParams , false, 1,P1(S(188))              , {2,P2(S(172),S(121))       , cMul        ,AnyParams       ,0}},
        /* 136:	(cMul (cSin [(cMul {% x})])@D5 (cPow [(cCos [(cMul {-% x})]) -1])@D5)
         *	:	(cTan [(cMul {% x})])
         */		 {ReplaceParams , false, 1,P1(S(181))              , {2,P2(S(164),S(113))       , cMul        ,AnyParams       ,0}},
        /* 137:	(cMul -1 (cTanh [(cMul %@N <1>)]))
         *	:	(cTanh [(cMul -% <1>)])
         */		 {ReplaceParams , false, 1,P1(S(192))              , {2,P2(N(2),S(191))         , cMul        ,AnyParams       ,0}},
        /* 138:	(cMul (cAdd {-1 (cPow [% x])})@D5 (cPow [(cAdd {1 (cPow [% x])}) -1])@D5)
         *	:	(cTanh [(cMul {x LOG( % ) 0.5})])
         */		 {ReplaceParams , false, 1,P1(S(189))              , {2,P2(S(198),S(137))       , cMul        ,AnyParams       ,0}},
        /* 139:	(cMul (cAdd {1 (cPow [% x])})@D5 (cPow [(cAdd {-1 (cPow [% x])}) -1])@D5)
         *	:	(cPow [(cTanh [(cMul {x LOG( % ) 0.5})]) -1])
         */		 {ReplaceParams , false, 1,P1(S(135))              , {2,P2(S(202),S(136))       , cMul        ,AnyParams       ,0}},
        /* 140:	(cMul (cLog [x]) 0.434294481903)
         *	:	(cLog10 [x])
         */		 {ReplaceParams , false, 1,P1(S(55))               , {2,P2(S(47),N(7))          , cMul        ,AnyParams       ,0}},
        /* 141:	(cMul (cPow [(cLog [x]) -1]) 2.30258509299)
         *	:	(cPow [(cLog10 [x]) -1])
         */		 {ReplaceParams , false, 1,P1(S(124))              , {2,P2(S(122),N(14))        , cMul        ,AnyParams       ,0}},
        /* 142:	(cMul (cLog [x]) 1.44269504089)
         *	:	(cLog2 [x])
         */		 {ReplaceParams , false, 1,P1(S(56))               , {2,P2(S(47),N(11))         , cMul        ,AnyParams       ,0}},
        /* 143:	(cMul (cPow [(cLog [x]) -1]) 0.69314718056)
         *	:	(cPow [(cLog2 [x]) -1])
         */		 {ReplaceParams , false, 1,P1(S(125))              , {2,P2(S(122),N(9))         , cMul        ,AnyParams       ,0}},
        /* 144:	(cMul (cPow [& y]) (cAdd {(cMul {% (cPow [x@P z])})@D1 %@D1}))
         *	:	% (cAdd {(cPow [& y]) (cPow [& (cAdd {y (cMul {z (cLog [x]) /LOG( & )})})])})
         */		 {ReplaceParams , false, 2,P2(P(1),S(229))         , {2,P2(S(84),S(227))        , cMul        ,AnyParams       ,0}},
        /* 145:	(cMul (cPow [& y]) (cAdd {(cMul {% (cPow [x@P z])})@D1 -%@D1}))
         *	:	% (cAdd {(cMul {(cPow [& y]) -1}) (cPow [& (cAdd {y (cMul {z (cLog [x]) /LOG( & )})})])})
         */		 {ReplaceParams , false, 2,P2(P(1),S(232))         , {2,P2(S(84),S(231))        , cMul        ,AnyParams       ,0}},
        /* 146:	(cMul (cSinh [x])@D4 (cPow [(cCosh [x]) %])@D4)
         *	:	(cTanh [x]) (cPow [(cCosh [x]) ADD( % 1 )])
         */		 {ReplaceParams , false, 2,P2(S(186),S(120))       , {2,P2(S(169),S(119))       , cMul        ,AnyParams       ,0}},
        /* 147:	@L (cMul (cAbs [x]))
         *	:	x
         */		 {ReplaceParams , true , 1,P1(P(14))               , {1,P1(S(0))                , cMul        ,AnyParams       ,0}},
        /* 148:	@L (cMul %@N)
         *	:	-%
         */		 {ReplaceParams , true , 1,P1(S(194))              , {1,P1(P(0))                , cMul        ,AnyParams       ,0}},
        /* 149:	(cEqual [(cAbs [x]) 0])
         *	:	x 0
         */		 {ReplaceParams , false, 2,P2(P(14),N(4))          , {2,P2(S(0),N(4))           , cEqual      ,PositionalParams,0}},
        /* 150:	(cEqual [(cAdd % <1>) &])
         *	:	(cAdd  <1>) ADD( & -% )
         */		 {ReplaceParams , false, 2,P2(S(242),S(260))       , {2,P2(S(248),P(9))         , cEqual      ,PositionalParams,0}},
        /* 151:	(cEqual [(cMul % <1>) &])
         *	:	(cMul  <1>) MUL( & /% )
         */		 {ReplaceParams , false, 2,P2(S(331),S(363))       , {2,P2(S(339),P(9))         , cEqual      ,PositionalParams,0}},
        /* 152:	(cEqual [(cPow [x %]) &])
         *	:	(cPow [(cPow [x %]) /%]) POW( & /% )
         */		 {ReplaceParams , false, 2,P2(S(126),S(156))       , {2,P2(S(77),P(9))          , cEqual      ,PositionalParams,0}},
        /* 153:	(cEqual [(cAdd % <1>) (cAdd & <2>)])
         *	:	(cAdd  <1>) (cAdd & -% <2>)
         */		 {ReplaceParams , false, 2,P2(S(242),S(254))       , {2,P2(S(248),S(250))       , cEqual      ,PositionalParams,0}},
        /* 154:	(cEqual [(cAdd x <1>)@D4 (cAdd x <2>)@D4])
         *	:	(cAdd  <1>) (cAdd  <2>)
         */		 {ReplaceParams , false, 2,P2(S(242),S(243))       , {2,P2(S(251),S(252))       , cEqual      ,PositionalParams,0}},
        /* 155:	(cEqual [(cMul % <1>) (cMul & <2>)])
         *	:	(cMul  <1>) (cMul & /% <2>)
         */		 {ReplaceParams , false, 2,P2(S(331),S(342))       , {2,P2(S(339),S(340))       , cEqual      ,PositionalParams,0}},
        /* 156:	(cEqual [(cPow [%@P x]) &@P])
         *	:	x MUL( LOG( & ) /LOG( % ) )
         */		 {ReplaceParams , false, 2,P2(P(14),S(364))        , {2,P2(S(81),P(13))         , cEqual      ,PositionalParams,0}},
        /* 157:	(cEqual [(cPow [x %@P]) (cPow [y &@P])])
         *	:	(cPow [(cPow [x %]) /MIN( % & )]) (cPow [(cPow [y &]) /MIN( % & )])
         */		 {ReplaceParams , false, 2,P2(S(127),S(128))       , {2,P2(S(78),S(79))         , cEqual      ,PositionalParams,0}},
        /* 158:	(cEqual [(cPow [% x])@D1 (cPow [% y])@D1])
         *	:	x y
         */		 {ReplaceParams , false, 2,P2(P(14),P(21))         , {2,P2(S(75),S(76))         , cEqual      ,PositionalParams,0}},
        /* 159:	(cNEqual [(cAbs [x]) 0])
         *	:	x 0
         */		 {ReplaceParams , false, 2,P2(P(14),N(4))          , {2,P2(S(0),N(4))           , cNEqual     ,PositionalParams,0}},
        /* 160:	(cNEqual [(cAdd % <1>) &])
         *	:	(cAdd  <1>) ADD( & -% )
         */		 {ReplaceParams , false, 2,P2(S(242),S(260))       , {2,P2(S(248),P(9))         , cNEqual     ,PositionalParams,0}},
        /* 161:	(cNEqual [(cMul % <1>) &])
         *	:	(cMul  <1>) MUL( & /% )
         */		 {ReplaceParams , false, 2,P2(S(331),S(363))       , {2,P2(S(339),P(9))         , cNEqual     ,PositionalParams,0}},
        /* 162:	(cNEqual [(cPow [x %]) &])
         *	:	(cPow [(cPow [x %]) /%]) POW( & /% )
         */		 {ReplaceParams , false, 2,P2(S(126),S(156))       , {2,P2(S(77),P(9))          , cNEqual     ,PositionalParams,0}},
        /* 163:	(cNEqual [(cAdd % <1>) (cAdd & <2>)])
         *	:	(cAdd  <1>) (cAdd & -% <2>)
         */		 {ReplaceParams , false, 2,P2(S(242),S(254))       , {2,P2(S(248),S(250))       , cNEqual     ,PositionalParams,0}},
        /* 164:	(cNEqual [(cAdd x <1>)@D4 (cAdd x <2>)@D4])
         *	:	(cAdd  <1>) (cAdd  <2>)
         */		 {ReplaceParams , false, 2,P2(S(242),S(243))       , {2,P2(S(251),S(252))       , cNEqual     ,PositionalParams,0}},
        /* 165:	(cNEqual [(cMul % <1>) (cMul & <2>)])
         *	:	(cMul  <1>) (cMul & /% <2>)
         */		 {ReplaceParams , false, 2,P2(S(331),S(342))       , {2,P2(S(339),S(340))       , cNEqual     ,PositionalParams,0}},
        /* 166:	(cNEqual [(cPow [%@P x]) &@P])
         *	:	x MUL( LOG( & ) /LOG( % ) )
         */		 {ReplaceParams , false, 2,P2(P(14),S(364))        , {2,P2(S(81),P(13))         , cNEqual     ,PositionalParams,0}},
        /* 167:	(cNEqual [(cPow [x %@P]) (cPow [y &@P])])
         *	:	(cPow [(cPow [x %]) /MIN( % & )]) (cPow [(cPow [y &]) /MIN( % & )])
         */		 {ReplaceParams , false, 2,P2(S(127),S(128))       , {2,P2(S(78),S(79))         , cNEqual     ,PositionalParams,0}},
        /* 168:	(cNEqual [(cPow [% x])@D1 (cPow [% y])@D1])
         *	:	x y
         */		 {ReplaceParams , false, 2,P2(P(14),P(21))         , {2,P2(S(75),S(76))         , cNEqual     ,PositionalParams,0}},
        /* 169:	(cLess [x 0.5])
         *	->	(cAbsNot [x])
         */		 {ProduceNewTree, false, 1,P1(S(409))              , {2,P2(P(14),N(8))          , cLess       ,PositionalParams,0}},
        /* 170:	(cLess [(cAbs [x]) %])
         *	->	(cNot [(cMul {x 0.5 /%})])
         */		 {ProduceNewTree, false, 1,P1(S(378))              , {2,P2(S(0),P(1))           , cLess       ,PositionalParams,0}},
        /* 171:	(cLess [(cMul %@P <1>) &])
         *	:	(cMul  <1>) MUL( & /% )
         */		 {ReplaceParams , false, 2,P2(S(331),S(363))       , {2,P2(S(336),P(9))         , cLess       ,PositionalParams,0}},
        /* 172:	(cLess [(cMul %@N <1>) &])
         *	:	MUL( & /% ) (cMul  <1>)
         */		 {ReplaceParams , false, 2,P2(S(363),S(331))       , {2,P2(S(335),P(9))         , cLess       ,PositionalParams,0}},
        /* 173:	(cLess [(cAdd % <1>) &])
         *	:	(cAdd  <1>) ADD( & -% )
         */		 {ReplaceParams , false, 2,P2(S(242),S(260))       , {2,P2(S(248),P(9))         , cLess       ,PositionalParams,0}},
        /* 174:	(cLess [(cPow [x %]) &])
         *	:	(cPow [(cPow [x %]) /%]) POW( & /% )
         */		 {ReplaceParams , false, 2,P2(S(126),S(156))       , {2,P2(S(77),P(9))          , cLess       ,PositionalParams,0}},
        /* 175:	(cLess [(cAdd % <1>) (cAdd & <2>)])
         *	:	(cAdd  <1>) (cAdd & -% <2>)
         */		 {ReplaceParams , false, 2,P2(S(242),S(254))       , {2,P2(S(248),S(250))       , cLess       ,PositionalParams,0}},
        /* 176:	(cLess [(cAdd x <1>)@D4 (cAdd x <2>)@D4])
         *	:	(cAdd  <1>) (cAdd  <2>)
         */		 {ReplaceParams , false, 2,P2(S(242),S(243))       , {2,P2(S(251),S(252))       , cLess       ,PositionalParams,0}},
        /* 177:	(cLess [(cMul %@P <1>) (cMul & <2>)])
         *	:	(cMul  <1>) (cMul & /% <2>)
         */		 {ReplaceParams , false, 2,P2(S(331),S(342))       , {2,P2(S(336),S(340))       , cLess       ,PositionalParams,0}},
        /* 178:	(cLess [(cMul %@N <1>) (cMul & <2>)])
         *	:	(cMul & /% <2>) (cMul  <1>)
         */		 {ReplaceParams , false, 2,P2(S(342),S(331))       , {2,P2(S(335),S(340))       , cLess       ,PositionalParams,0}},
        /* 179:	(cLess [(cPow [%@P x]) &@P])
         *	:	x MUL( LOG( & ) /LOG( % ) )
         */		 {ReplaceParams , false, 2,P2(P(14),S(364))        , {2,P2(S(81),P(13))         , cLess       ,PositionalParams,0}},
        /* 180:	(cLess [(cPow [x %@P]) (cPow [y &@P])])
         *	:	(cPow [(cPow [x %]) /MIN( % & )]) (cPow [(cPow [y &]) /MIN( % & )])
         */		 {ReplaceParams , false, 2,P2(S(127),S(128))       , {2,P2(S(78),S(79))         , cLess       ,PositionalParams,0}},
        /* 181:	(cLess [(cPow [% x])@D1 (cPow [% y])@D1])
         *	:	x y
         */		 {ReplaceParams , false, 2,P2(P(14),P(21))         , {2,P2(S(75),S(76))         , cLess       ,PositionalParams,0}},
        /* 182:	(cLessOrEq [% (cAbs [x])])
         *	->	(cNotNot [(cMul {x 0.5 /%})])
         */		 {ProduceNewTree, false, 1,P1(S(403))              , {2,P2(P(1),S(0))           , cLessOrEq   ,PositionalParams,0}},
        /* 183:	(cLessOrEq [(cMul %@P <1>) &])
         *	:	(cMul  <1>) MUL( & /% )
         */		 {ReplaceParams , false, 2,P2(S(331),S(363))       , {2,P2(S(336),P(9))         , cLessOrEq   ,PositionalParams,0}},
        /* 184:	(cLessOrEq [(cMul %@N <1>) &])
         *	:	MUL( & /% ) (cMul  <1>)
         */		 {ReplaceParams , false, 2,P2(S(363),S(331))       , {2,P2(S(335),P(9))         , cLessOrEq   ,PositionalParams,0}},
        /* 185:	(cLessOrEq [(cAdd % <1>) &])
         *	:	(cAdd  <1>) ADD( & -% )
         */		 {ReplaceParams , false, 2,P2(S(242),S(260))       , {2,P2(S(248),P(9))         , cLessOrEq   ,PositionalParams,0}},
        /* 186:	(cLessOrEq [(cPow [x %]) &])
         *	:	(cPow [(cPow [x %]) /%]) POW( & /% )
         */		 {ReplaceParams , false, 2,P2(S(126),S(156))       , {2,P2(S(77),P(9))          , cLessOrEq   ,PositionalParams,0}},
        /* 187:	(cLessOrEq [(cAdd % <1>) (cAdd & <2>)])
         *	:	(cAdd  <1>) (cAdd & -% <2>)
         */		 {ReplaceParams , false, 2,P2(S(242),S(254))       , {2,P2(S(248),S(250))       , cLessOrEq   ,PositionalParams,0}},
        /* 188:	(cLessOrEq [(cAdd x <1>)@D4 (cAdd x <2>)@D4])
         *	:	(cAdd  <1>) (cAdd  <2>)
         */		 {ReplaceParams , false, 2,P2(S(242),S(243))       , {2,P2(S(251),S(252))       , cLessOrEq   ,PositionalParams,0}},
        /* 189:	(cLessOrEq [(cMul %@P <1>) (cMul & <2>)])
         *	:	(cMul  <1>) (cMul & /% <2>)
         */		 {ReplaceParams , false, 2,P2(S(331),S(342))       , {2,P2(S(336),S(340))       , cLessOrEq   ,PositionalParams,0}},
        /* 190:	(cLessOrEq [(cMul %@N <1>) (cMul & <2>)])
         *	:	(cMul & /% <2>) (cMul  <1>)
         */		 {ReplaceParams , false, 2,P2(S(342),S(331))       , {2,P2(S(335),S(340))       , cLessOrEq   ,PositionalParams,0}},
        /* 191:	(cLessOrEq [(cPow [%@P x]) &@P])
         *	:	x MUL( LOG( & ) /LOG( % ) )
         */		 {ReplaceParams , false, 2,P2(P(14),S(364))        , {2,P2(S(81),P(13))         , cLessOrEq   ,PositionalParams,0}},
        /* 192:	(cLessOrEq [(cPow [x %@P]) (cPow [y &@P])])
         *	:	(cPow [(cPow [x %]) /MIN( % & )]) (cPow [(cPow [y &]) /MIN( % & )])
         */		 {ReplaceParams , false, 2,P2(S(127),S(128))       , {2,P2(S(78),S(79))         , cLessOrEq   ,PositionalParams,0}},
        /* 193:	(cLessOrEq [(cPow [% x])@D1 (cPow [% y])@D1])
         *	:	x y
         */		 {ReplaceParams , false, 2,P2(P(14),P(21))         , {2,P2(S(75),S(76))         , cLessOrEq   ,PositionalParams,0}},
        /* 194:	(cGreater [(cMul %@P <1>) &])
         *	:	(cMul  <1>) MUL( & /% )
         */		 {ReplaceParams , false, 2,P2(S(331),S(363))       , {2,P2(S(336),P(9))         , cGreater    ,PositionalParams,0}},
        /* 195:	(cGreater [(cMul %@N <1>) &])
         *	:	MUL( & /% ) (cMul  <1>)
         */		 {ReplaceParams , false, 2,P2(S(363),S(331))       , {2,P2(S(335),P(9))         , cGreater    ,PositionalParams,0}},
        /* 196:	(cGreater [(cAdd % <1>) &])
         *	:	(cAdd  <1>) ADD( & -% )
         */		 {ReplaceParams , false, 2,P2(S(242),S(260))       , {2,P2(S(248),P(9))         , cGreater    ,PositionalParams,0}},
        /* 197:	(cGreater [(cPow [x %]) &])
         *	:	(cPow [(cPow [x %]) /%]) POW( & /% )
         */		 {ReplaceParams , false, 2,P2(S(126),S(156))       , {2,P2(S(77),P(9))          , cGreater    ,PositionalParams,0}},
        /* 198:	(cGreater [(cAdd % <1>) (cAdd & <2>)])
         *	:	(cAdd  <1>) (cAdd & -% <2>)
         */		 {ReplaceParams , false, 2,P2(S(242),S(254))       , {2,P2(S(248),S(250))       , cGreater    ,PositionalParams,0}},
        /* 199:	(cGreater [(cAdd x <1>)@D4 (cAdd x <2>)@D4])
         *	:	(cAdd  <1>) (cAdd  <2>)
         */		 {ReplaceParams , false, 2,P2(S(242),S(243))       , {2,P2(S(251),S(252))       , cGreater    ,PositionalParams,0}},
        /* 200:	(cGreater [(cMul %@P <1>) (cMul & <2>)])
         *	:	(cMul  <1>) (cMul & /% <2>)
         */		 {ReplaceParams , false, 2,P2(S(331),S(342))       , {2,P2(S(336),S(340))       , cGreater    ,PositionalParams,0}},
        /* 201:	(cGreater [(cMul %@N <1>) (cMul & <2>)])
         *	:	(cMul & /% <2>) (cMul  <1>)
         */		 {ReplaceParams , false, 2,P2(S(342),S(331))       , {2,P2(S(335),S(340))       , cGreater    ,PositionalParams,0}},
        /* 202:	(cGreater [(cPow [%@P x]) &@P])
         *	:	x MUL( LOG( & ) /LOG( % ) )
         */		 {ReplaceParams , false, 2,P2(P(14),S(364))        , {2,P2(S(81),P(13))         , cGreater    ,PositionalParams,0}},
        /* 203:	(cGreater [(cPow [x %@P]) (cPow [y &@P])])
         *	:	(cPow [(cPow [x %]) /MIN( % & )]) (cPow [(cPow [y &]) /MIN( % & )])
         */		 {ReplaceParams , false, 2,P2(S(127),S(128))       , {2,P2(S(78),S(79))         , cGreater    ,PositionalParams,0}},
        /* 204:	(cGreater [(cPow [% x])@D1 (cPow [% y])@D1])
         *	:	x y
         */		 {ReplaceParams , false, 2,P2(P(14),P(21))         , {2,P2(S(75),S(76))         , cGreater    ,PositionalParams,0}},
        /* 205:	(cGreaterOrEq [x 0.5])
         *	->	(cAbsNotNot [x])
         */		 {ProduceNewTree, false, 1,P1(S(410))              , {2,P2(P(14),N(8))          , cGreaterOrEq,PositionalParams,0}},
        /* 206:	(cGreaterOrEq [(cMul %@P <1>) &])
         *	:	(cMul  <1>) MUL( & /% )
         */		 {ReplaceParams , false, 2,P2(S(331),S(363))       , {2,P2(S(336),P(9))         , cGreaterOrEq,PositionalParams,0}},
        /* 207:	(cGreaterOrEq [(cMul %@N <1>) &])
         *	:	MUL( & /% ) (cMul  <1>)
         */		 {ReplaceParams , false, 2,P2(S(363),S(331))       , {2,P2(S(335),P(9))         , cGreaterOrEq,PositionalParams,0}},
        /* 208:	(cGreaterOrEq [(cAdd % <1>) &])
         *	:	(cAdd  <1>) ADD( & -% )
         */		 {ReplaceParams , false, 2,P2(S(242),S(260))       , {2,P2(S(248),P(9))         , cGreaterOrEq,PositionalParams,0}},
        /* 209:	(cGreaterOrEq [(cPow [x %]) &])
         *	:	(cPow [(cPow [x %]) /%]) POW( & /% )
         */		 {ReplaceParams , false, 2,P2(S(126),S(156))       , {2,P2(S(77),P(9))          , cGreaterOrEq,PositionalParams,0}},
        /* 210:	(cGreaterOrEq [(cAdd % <1>) (cAdd & <2>)])
         *	:	(cAdd  <1>) (cAdd & -% <2>)
         */		 {ReplaceParams , false, 2,P2(S(242),S(254))       , {2,P2(S(248),S(250))       , cGreaterOrEq,PositionalParams,0}},
        /* 211:	(cGreaterOrEq [(cAdd x <1>)@D4 (cAdd x <2>)@D4])
         *	:	(cAdd  <1>) (cAdd  <2>)
         */		 {ReplaceParams , false, 2,P2(S(242),S(243))       , {2,P2(S(251),S(252))       , cGreaterOrEq,PositionalParams,0}},
        /* 212:	(cGreaterOrEq [(cMul %@P <1>) (cMul & <2>)])
         *	:	(cMul  <1>) (cMul & /% <2>)
         */		 {ReplaceParams , false, 2,P2(S(331),S(342))       , {2,P2(S(336),S(340))       , cGreaterOrEq,PositionalParams,0}},
        /* 213:	(cGreaterOrEq [(cMul %@N <1>) (cMul & <2>)])
         *	:	(cMul & /% <2>) (cMul  <1>)
         */		 {ReplaceParams , false, 2,P2(S(342),S(331))       , {2,P2(S(335),S(340))       , cGreaterOrEq,PositionalParams,0}},
        /* 214:	(cGreaterOrEq [(cPow [%@P x]) &@P])
         *	:	x MUL( LOG( & ) /LOG( % ) )
         */		 {ReplaceParams , false, 2,P2(P(14),S(364))        , {2,P2(S(81),P(13))         , cGreaterOrEq,PositionalParams,0}},
        /* 215:	(cGreaterOrEq [(cPow [x %@P]) (cPow [y &@P])])
         *	:	(cPow [(cPow [x %]) /MIN( % & )]) (cPow [(cPow [y &]) /MIN( % & )])
         */		 {ReplaceParams , false, 2,P2(S(127),S(128))       , {2,P2(S(78),S(79))         , cGreaterOrEq,PositionalParams,0}},
        /* 216:	(cGreaterOrEq [(cPow [% x])@D1 (cPow [% y])@D1])
         *	:	x y
         */		 {ReplaceParams , false, 2,P2(P(14),P(21))         , {2,P2(S(75),S(76))         , cGreaterOrEq,PositionalParams,0}},
        /* 217:	(cNot [x@P])
         *	->	(cAbsNot [x])
         */		 {ProduceNewTree, false, 1,P1(S(409))              , {1,P1(P(17))               , cNot        ,PositionalParams,0}},
        /* 218:	(cAnd x@L <1>)
         *	->	(cNotNot [(cMul {x (cAnd  <1>)})])
         */		 {ProduceNewTree, false, 1,P1(S(404))              , {1,P1(P(19))               , cAnd        ,AnyParams       ,1}},
        /* 219:	(cAnd x@P y@P <1>)
         *	->	(cAbsAnd {x y (cAnd  <1>)})
         */		 {ProduceNewTree, false, 1,P1(S(407))              , {2,P2(P(17),P(27))         , cAnd        ,AnyParams       ,1}},
        /* 220:	(cAnd x y)
         *	->	(cIf [x (cNotNot [y]) 0])
         */		 {ProduceNewTree, false, 1,P1(S(45))               , {2,P2(P(14),P(21))         , cAnd        ,AnyParams       ,0}},
        /* 221:	(cAnd (cIf [x y z])@D4 (cIf [x a b])@D4)
         *	:	(cIf [x (cAnd {y a}) (cAnd {z b})])
         */		 {ReplaceParams , false, 1,P1(S(43))               , {2,P2(S(32),S(36))         , cAnd        ,AnyParams       ,0}},
        /* 222:	(cAnd (cNot [x]) (cNot [y]))
         *	:	(cNot [(cOr {x y})])
         */		 {ReplaceParams , false, 1,P1(S(382))              , {2,P2(S(375),S(376))       , cAnd        ,AnyParams       ,0}},
        /* 223:	(cAnd (cNot [z]) (cIf [x (cNot [y]) %@L]))
         *	:	(cNot [(cOr {z (cIf [x y (cNot [%])])})])
         */		 {ReplaceParams , false, 1,P1(S(383))              , {2,P2(S(377),S(38))        , cAnd        ,AnyParams       ,0}},
        /* 224:	(cAnd (cNot [z]) (cIf [x %@L (cNot [y])]))
         *	:	(cNot [(cOr {z (cIf [x (cNot [%]) y])})])
         */		 {ReplaceParams , false, 1,P1(S(384))              , {2,P2(S(377),S(34))        , cAnd        ,AnyParams       ,0}},
        /* 225:	(cAnd (cEqual [x y])@D12 (cEqual [y z])@D24 (cEqual [x z])@D20)
         *	:	(cEqual [x y]) (cEqual [y z])
         */		 {ReplaceParams , false, 2,P2(S(366),S(369))       , {3,P3(S(365),S(368),S(367)), cAnd        ,AnyParams       ,0}},
        /* 226:	(cOr x@P y@P <1>)
         *	->	(cAbsOr {x y (cOr  <1>)})
         */		 {ProduceNewTree, false, 1,P1(S(408))              , {2,P2(P(17),P(27))         , cOr         ,AnyParams       ,1}},
        /* 227:	(cOr x y)
         *	->	(cIf [x 1 (cNotNot [y])])
         */		 {ProduceNewTree, false, 1,P1(S(31))               , {2,P2(P(14),P(21))         , cOr         ,AnyParams       ,0}},
        /* 228:	(cOr x@L y@L)
         *	:	(cNotNot [(cAdd {x y})])
         */		 {ReplaceParams , false, 1,P1(S(401))              , {2,P2(P(19),P(22))         , cOr         ,AnyParams       ,0}},
        /* 229:	(cOr (cIf [x y z])@D4 (cIf [x a b])@D4)
         *	:	(cIf [x (cOr {y a}) (cOr {z b})])
         */		 {ReplaceParams , false, 1,P1(S(44))               , {2,P2(S(32),S(36))         , cOr         ,AnyParams       ,0}},
        /* 230:	(cOr (cNot [x]) (cNot [y]))
         *	:	(cNot [(cAnd {x y})])
         */		 {ReplaceParams , false, 1,P1(S(379))              , {2,P2(S(375),S(376))       , cOr         ,AnyParams       ,0}},
        /* 231:	(cOr (cNot [z]) (cIf [x (cNot [y]) %@L]))
         *	:	(cNot [(cAnd {z (cIf [x y (cNot [%])])})])
         */		 {ReplaceParams , false, 1,P1(S(380))              , {2,P2(S(377),S(38))        , cOr         ,AnyParams       ,0}},
        /* 232:	(cOr (cNot [z]) (cIf [x %@L (cNot [y])]))
         *	:	(cNot [(cAnd {z (cIf [x (cNot [%]) y])})])
         */		 {ReplaceParams , false, 1,P1(S(381))              , {2,P2(S(377),S(34))        , cOr         ,AnyParams       ,0}},
        /* 233:	(cOr x@L (cAdd  <1>)@P)
         *	:	(cNotNot [(cAdd x <1>)])
         */		 {ReplaceParams , false, 1,P1(S(402))              , {2,P2(P(19),S(244))        , cOr         ,AnyParams       ,0}},
        /* 234:	(cNotNot [x@P])
         *	->	(cAbsNotNot [x])
         */		 {ProduceNewTree, false, 1,P1(S(410))              , {1,P1(P(17))               , cNotNot     ,PositionalParams,0}},
        /* 235:	@L (cNotNot [x])
         *	->	x
         */		 {ProduceNewTree, true , 1,P1(P(14))               , {1,P1(P(14))               , cNotNot     ,PositionalParams,0}},
        /* 236:	(cAbsAnd x y)
         *	->	(cAbsIf [x (cAbsNotNot [y]) 0])
         */		 {ProduceNewTree, false, 1,P1(S(416))              , {2,P2(P(14),P(21))         , cAbsAnd     ,AnyParams       ,0}},
        /* 237:	(cAbsOr x y)
         *	->	(cAbsIf [x 1 (cAbsNotNot [y])])
         */		 {ProduceNewTree, false, 1,P1(S(413))              , {2,P2(P(14),P(21))         , cAbsOr      ,AnyParams       ,0}},
        /* 238:	@L (cAbsNotNot [x])
         *	->	x
         */		 {ProduceNewTree, true , 1,P1(P(14))               , {1,P1(P(14))               , cAbsNotNot  ,PositionalParams,0}},
        /* 239:	(cAbsIf [(cNotNot [x]) y z])
         *	->	(cIf [x y z])
         */		 {ProduceNewTree, false, 1,P1(S(33))               , {3,P3(S(399),P(21),P(28))  , cAbsIf      ,PositionalParams,0}},
        /* 240:	(cAbsIf [(cLessOrEq [x y]) z a])
         *	:	(cLess [y x]) a z
         */		 {ReplaceParams , false, 3,P3(S(371),P(31),P(28))  , {3,P3(S(373),P(28),P(31))  , cAbsIf      ,PositionalParams,0}},
    };
#undef P
#undef N
#undef S

    struct grammar_optimize_abslogical_type
    {
        unsigned c;
        unsigned char l[12];
    };
    extern "C"
    {
        grammar_optimize_abslogical_type grammar_optimize_abslogical =
        {
            12,
            { 19,29,30,169,205,217,219,226,234,238,
              239,240
    }   };  }
    struct grammar_optimize_ignore_if_sideeffects_type
    {
        unsigned c;
        unsigned char l[34];
    };
    extern "C"
    {
        grammar_optimize_ignore_if_sideeffects_type grammar_optimize_ignore_if_sideeffects =
        {
            34,
            { 0,18,20,21,22,23,24,25,26,27,
              28,37,39,75,109,110,111,112,147,148,
              149,159,170,182,221,222,223,224,225,229,
              230,231,232,235
    }   };  }
    struct grammar_optimize_nonshortcut_logical_evaluation_type
    {
        unsigned c;
        unsigned char l[31];
    };
    extern "C"
    {
        grammar_optimize_nonshortcut_logical_evaluation_type grammar_optimize_nonshortcut_logical_evaluation =
        {
            31,
            { 0,25,26,27,28,37,39,75,109,110,
              111,112,147,148,149,159,170,182,218,221,
              222,223,224,225,228,229,230,231,232,233,
              235
    }   };  }
    struct grammar_optimize_round1_type
    {
        unsigned c;
        unsigned char l[92];
    };
    extern "C"
    {
        grammar_optimize_round1_type grammar_optimize_round1 =
        {
            92,
            { 0,1,2,3,4,5,6,7,8,9,
              10,11,15,16,25,26,27,28,31,32,
              33,36,37,38,39,40,41,42,43,44,
              45,46,47,48,49,53,54,55,56,57,
              58,59,60,61,62,63,70,75,76,77,
              78,79,80,81,82,83,84,85,86,87,
              88,104,108,109,110,111,112,113,114,115,
              116,117,118,119,144,145,147,148,149,159,
              170,182,221,222,223,224,225,229,230,231,
              232,235
    }   };  }
    struct grammar_optimize_round2_type
    {
        unsigned c;
        unsigned char l[81];
    };
    extern "C"
    {
        grammar_optimize_round2_type grammar_optimize_round2 =
        {
            81,
            { 0,5,12,13,14,15,25,26,27,28,
              31,34,35,36,37,38,39,40,42,43,
              44,45,46,47,48,49,57,58,59,64,
              65,70,71,72,73,74,75,89,90,91,
              92,93,94,95,96,97,98,99,100,101,
              102,103,108,109,110,111,112,113,116,117,
              118,120,144,145,146,147,148,149,159,170,
              182,221,222,223,224,225,229,230,231,232,
              235
    }   };  }
    struct grammar_optimize_round3_type
    {
        unsigned c;
        unsigned char l[85];
    };
    extern "C"
    {
        grammar_optimize_round3_type grammar_optimize_round3 =
        {
            85,
            { 66,67,68,69,121,122,123,124,125,126,
              127,128,129,130,131,132,133,134,135,136,
              137,138,139,150,151,152,153,154,155,156,
              157,158,160,161,162,163,164,165,166,167,
              168,171,172,173,174,175,176,177,178,179,
              180,181,183,184,185,186,187,188,189,190,
              191,192,193,194,195,196,197,198,199,200,
              201,202,203,204,206,207,208,209,210,211,
              212,213,214,215,216
    }   };  }
    struct grammar_optimize_round4_type
    {
        unsigned c;
        unsigned char l[10];
    };
    extern "C"
    {
        grammar_optimize_round4_type grammar_optimize_round4 =
        {
            10,
            { 17,50,51,52,106,107,140,141,142,143
    }   };  }
    struct grammar_optimize_shortcut_logical_evaluation_type
    {
        unsigned c;
        unsigned char l[33];
    };
    extern "C"
    {
        grammar_optimize_shortcut_logical_evaluation_type grammar_optimize_shortcut_logical_evaluation =
        {
            33,
            { 0,25,26,27,28,37,39,75,105,109,
              110,111,112,147,148,149,159,170,182,220,
              221,222,223,224,225,227,229,230,231,232,
              235,236,237
    }   };  }
}
#undef P1
#undef P2
#undef P3

namespace FPoptimizer_Grammar
{
    ParamSpec ParamSpec_Extract(unsigned paramlist, unsigned index)
    {
        index = (paramlist >> (index * PARAM_INDEX_BITS)) % (1 << PARAM_INDEX_BITS);
        const unsigned p_begin = 0;
        const unsigned n_begin = p_begin + sizeof(plist.plist_p)/sizeof(*plist.plist_p);
        const unsigned s_begin = n_begin + sizeof(plist.plist_n)/sizeof(*plist.plist_n);
      /*const unsigned     end = s_begin + sizeof(plist.plist_s)/sizeof(*plist.plist_s);*/
        if(index < s_begin)
        {
            if(index < n_begin)
                return ParamSpec(ParamHolder,(const void*)&plist.plist_p[index-p_begin]);
            else
                return ParamSpec(NumConstant,(const void*)&plist.plist_n[index-n_begin]);
        }
        else
            return ParamSpec(SubFunction,(const void*)&plist.plist_s[index-s_begin]);
    }
}

#line 1 "fpoptimizer/fpoptimizer_optimize.c"
#include "fpconfig.h"
#include "fparser.h"
#include "fptypes.h"

#ifdef FP_SUPPORT_OPTIMIZER

// line removed for fpoptimizer.c: #include "fpoptimizer_grammar.h"
// line removed for fpoptimizer.c: #include "fpoptimizer_consts.h"
// line removed for fpoptimizer.c: #include "fpoptimizer_opcodename.h"
// line removed for fpoptimizer.c: #include "fpoptimizer_optimize.h"

#include <stdio.h>

#include <algorithm>
#include <map>
#include <sstream>

using namespace FUNCTIONPARSERTYPES;
using namespace FPoptimizer_Grammar;
using namespace FPoptimizer_CodeTree;
using namespace FPoptimizer_Optimize;

namespace
{
    /* I have heard that std::equal_range() is practically worthless
     * due to the insane limitation that the two parameters for Comp() must
     * be of the same type. Hence we must reinvent the wheel and implement
     * our own here. This is practically identical to the one from
     * GNU libstdc++, except rewritten. -Bisqwit
     */
    template<typename It, typename T, typename Comp>
    std::pair<It, It>
    MyEqualRange(It first, It last, const T& val, Comp comp)
    {
        size_t len = last-first;
        while(len > 0)
        {
            size_t half = len/2;
            It middle(first); middle += half;
            if(comp(*middle, val))
            {
                first = middle;
                ++first;
                len = len - half - 1;
            }
            else if(comp(val, *middle))
            {
                len = half;
            }
            else
            {
                // The following implements this:
                // // left = lower_bound(first, middle, val, comp);
                It left(first);
              {///
                It& first2 = left;
                It last2(middle);
                size_t len2 = last2-first2;
                while(len2 > 0)
                {
                    size_t half2 = len2 / 2;
                    It middle2(first2); middle2 += half2;
                    if(comp(*middle2, val))
                    {
                        first2 = middle2;
                        ++first2;
                        len2 = len2 - half2 - 1;
                    }
                    else
                        len2 = half2;
                }
                // left = first2;  - not needed, already happens due to reference
              }///
                first += len;
                // The following implements this:
                // // right = upper_bound(++middle, first, val, comp);
                It right(++middle);
              {///
                It& first2 = right;
                It& last2 = first;
                size_t len2 = last2-first2;
                while(len2 > 0)
                {
                    size_t half2 = len2 / 2;
                    It middle2(first2); middle2 += half2;
                    if(comp(val, *middle2))
                        len2 = half2;
                    else
                    {
                        first2 = middle2;
                        ++first2;
                        len2 = len2 - half2 - 1;
                    }
                }
                // right = first2;  - not needed, already happens due to reference
              }///
                return std::pair<It,It> (left,right);
            }
        }
        return std::pair<It,It> (first,first);
    }

    /* A helper for std::equal_range */
    struct OpcodeRuleCompare
    {
        bool operator() (const CodeTree& tree,
                         unsigned rulenumber) const
        {
            /* If this function returns true, len=half.
             */
            const Rule& rule = grammar_rules[rulenumber];
            return tree.GetOpcode() < rule.match_tree.subfunc_opcode;
        }
        bool operator() (unsigned rulenumber,
                         const CodeTree& tree) const
        {
            /* If this function returns true, rule will be excluded from the equal_range
             */
            const Rule& rule = grammar_rules[rulenumber];
            return rule.match_tree.subfunc_opcode < tree.GetOpcode();
        }
    };

    /* Test and apply a rule to a given CodeTree */
    bool TestRuleAndApplyIfMatch(
        const Rule& rule,
        CodeTree& tree,
        bool from_logical_context)
    {
        MatchInfo info;

        MatchResultType found(false, MatchPositionSpecBaseP());

        if(rule.logical_context && !from_logical_context)
        {
            /* If the rule only applies in logical contexts,
             * but we do not have a logical context, fail the rule
             */
            goto fail;
        }

        /*std::cout << "TESTING: ";
        DumpMatch(rule, *tree, info, false);*/

        for(;;)
        {
        #ifdef DEBUG_SUBSTITUTIONS
            //DumpMatch(rule, tree, info, "Testing");
        #endif
            found = TestParams(rule.match_tree, tree, found.specs, info, true);
            if(found.found) break;
            if(!&*found.specs)
            {
            fail:;
                // Did not match
        #ifdef DEBUG_SUBSTITUTIONS
                DumpMatch(rule, tree, info, false);
        #endif
                return false;
            }
        }
        // Matched
    #ifdef DEBUG_SUBSTITUTIONS
        DumpMatch(rule, tree, info, true);
    #endif
        SynthesizeRule(rule, tree, info);
        return true;
    }
}

namespace FPoptimizer_Optimize
{
    /* Apply the grammar to a given CodeTree */
    bool ApplyGrammar(
        const Grammar& grammar,
        CodeTree& tree,
        bool from_logical_context)
    {
        if(tree.GetOptimizedUsing() == &grammar)
        {
#ifdef DEBUG_SUBSTITUTIONS
            std::cout << "Already optimized:  ";
            DumpTree(tree);
            std::cout << "\n" << std::flush;
#endif
            return false;
        }

        /* First optimize all children */
        if(true)
        {
            bool changed = false;

            switch(tree.GetOpcode())
            {
                case cNot:
                case cNotNot:
                case cAnd:
                case cOr:
                    for(size_t a=0; a<tree.GetParamCount(); ++a)
                        if(ApplyGrammar( grammar, tree.GetParam(a), true))
                            changed = true;
                    break;
                case cIf:
                case cAbsIf:
                    if(ApplyGrammar( grammar, tree.GetParam(0), tree.GetOpcode() == cIf))
                        changed = true;
                    for(size_t a=1; a<tree.GetParamCount(); ++a)
                        if(ApplyGrammar( grammar, tree.GetParam(a), from_logical_context))
                            changed = true;
                    break;
                default:
                    for(size_t a=0; a<tree.GetParamCount(); ++a)
                        if(ApplyGrammar( grammar, tree.GetParam(a), false))
                            changed = true;
            }

            if(changed)
            {
                // Give the parent node a rerun at optimization
                tree.Mark_Incompletely_Hashed();
                return true;
            }
        }

        /* Figure out which rules _may_ match this tree */
        typedef const unsigned char* rulenumit;

        std::pair<rulenumit, rulenumit> range =
            MyEqualRange(grammar.rule_list,
                         grammar.rule_list + grammar.rule_count,
                         tree,
                         OpcodeRuleCompare());

        if(range.first != range.second)
        {
#ifdef DEBUG_SUBSTITUTIONS
            std::vector<unsigned char> rules;
            rules.reserve(range.second - range.first);
            for(rulenumit r = range.first; r != range.second; ++r)
            {
                //if(grammar_rules[*r].match_tree.subfunc_opcode != tree.GetOpcode()) continue;
                if(IsLogisticallyPlausibleParamsMatch(grammar_rules[*r].match_tree, tree))
                    rules.push_back(*r);
            }
            range.first = &rules[0];
            range.second = &rules[rules.size()-1]+1;

            if(range.first != range.second)
            {
                std::cout << "Input (" << FP_GetOpcodeName(tree.GetOpcode())
                          << ")[" << tree.GetParamCount()
                          << "]";
                if(from_logical_context)
                    std::cout << "(Logical)";

                unsigned first=~unsigned(0), prev=~unsigned(0);
                const char* sep = ", rules ";
                for(rulenumit r = range.first; r != range.second; ++r)
                {
                    if(first==~unsigned(0)) first=prev=*r;
                    else if(*r == prev+1) prev=*r;
                    else
                    {
                        std::cout << sep << first; sep=",";
                        if(prev != first) std::cout << '-' << prev;
                        first = prev = *r;
                    }
                }
                if(first != ~unsigned(0))
                {
                    std::cout << sep << first;
                    if(prev != first) std::cout << '-' << prev;
                }
                std::cout << ": ";
                DumpTree(tree);
                std::cout << "\n" << std::flush;
            }
#endif

            bool changed = false;

            for(rulenumit r = range.first; r != range.second; ++r)
            {
            #ifndef DEBUG_SUBSTITUTIONS
                if(!IsLogisticallyPlausibleParamsMatch(grammar_rules[*r].match_tree, tree))
                    continue;
            #endif
                if(TestRuleAndApplyIfMatch(grammar_rules[*r], tree, from_logical_context))
                {
                    changed = true;
                    break;
                }
            }

            if(changed)
            {
    #ifdef DEBUG_SUBSTITUTIONS
                std::cout << "Changed." << std::endl;
                std::cout << "Output: ";
                DumpTree(tree);
                std::cout << "\n" << std::flush;
    #endif
                // Give the parent node a rerun at optimization
                tree.Mark_Incompletely_Hashed();
                return true;
            }
        }

        // No changes, consider the tree properly optimized.
        tree.SetOptimizedUsing(&grammar);
        return false;
    }

    void ApplyGrammars(FPoptimizer_CodeTree::CodeTree& tree)
    {
    #ifdef FPOPTIMIZER_MERGED_FILE
        #define C *(const Grammar*)&
    #else
        #define C
    #endif
        #ifdef DEBUG_SUBSTITUTIONS
        std::cout << "Applying grammar_optimize_round1\n";
        #endif
        while(ApplyGrammar(C grammar_optimize_round1, tree))
            { //std::cout << "Rerunning 1\n";
                tree.FixIncompleteHashes();
            }

        #ifdef DEBUG_SUBSTITUTIONS
        std::cout << "Applying grammar_optimize_round2\n";
        #endif
        while(ApplyGrammar(C grammar_optimize_round2, tree))
            { //std::cout << "Rerunning 2\n";
                tree.FixIncompleteHashes();
            }

        #ifdef DEBUG_SUBSTITUTIONS
        std::cout << "Applying grammar_optimize_round3\n";
        #endif
        while(ApplyGrammar(C grammar_optimize_round3, tree))
            { //std::cout << "Rerunning 3\n";
                tree.FixIncompleteHashes();
            }

        #ifndef FP_ENABLE_SHORTCUT_LOGICAL_EVALUATION
        #ifdef DEBUG_SUBSTITUTIONS
        std::cout << "Applying grammar_optimize_nonshortcut_logical_evaluation\n";
        #endif
        while(ApplyGrammar(C grammar_optimize_nonshortcut_logical_evaluation, tree))
            { //std::cout << "Rerunning 3\n";
                tree.FixIncompleteHashes();
            }
        #endif

        #ifdef DEBUG_SUBSTITUTIONS
        std::cout << "Applying grammar_optimize_round4\n";
        #endif
        while(ApplyGrammar(C grammar_optimize_round4, tree))
            { //std::cout << "Rerunning 4\n";
                tree.FixIncompleteHashes();
            }

        #ifdef FP_ENABLE_SHORTCUT_LOGICAL_EVALUATION
        #ifdef DEBUG_SUBSTITUTIONS
        std::cout << "Applying grammar_optimize_shortcut_logical_evaluation\n";
        #endif
        while(ApplyGrammar(C grammar_optimize_shortcut_logical_evaluation, tree))
            { //std::cout << "Rerunning 3\n";
                tree.FixIncompleteHashes();
            }
        #endif

        #ifdef FP_ENABLE_IGNORE_IF_SIDEEFFECTS
        #ifdef DEBUG_SUBSTITUTIONS
        std::cout << "Applying grammar_optimize_ignore_if_sideeffects\n";
        #endif
        while(ApplyGrammar(C grammar_optimize_ignore_if_sideeffects, tree))
            { //std::cout << "Rerunning 3\n";
                tree.FixIncompleteHashes();
            }
        #endif

        #ifdef DEBUG_SUBSTITUTIONS
        std::cout << "Applying grammar_optimize_abslogical\n";
        #endif
        while(ApplyGrammar(C grammar_optimize_abslogical, tree))
            { //std::cout << "Rerunning 3\n";
                tree.FixIncompleteHashes();
            }

        #undef C
    }
}

#endif

#line 1 "fpoptimizer/fpoptimizer_optimize_match.c"
#include "fpconfig.h"
#include "fparser.h"
#include "fptypes.h"

#ifdef FP_SUPPORT_OPTIMIZER

#include <algorithm>
#include <assert.h>
#include <cstring>
#include <cmath>

#include <memory> /* for auto_ptr */

// line removed for fpoptimizer.c: #include "fpoptimizer_grammar.h"
// line removed for fpoptimizer.c: #include "fpoptimizer_optimize.h"

using namespace FUNCTIONPARSERTYPES;
using namespace FPoptimizer_Grammar;
using namespace FPoptimizer_CodeTree;
using namespace FPoptimizer_Optimize;

namespace
{
    /* Test the given constraints to a given CodeTree */
    bool TestImmedConstraints(unsigned bitmask, const CodeTree& tree)
    {
        switch(bitmask & ValueMask)
        {
            case Value_AnyNum: case ValueMask: break;
            case Value_EvenInt:
                if(tree.GetEvennessInfo() != CodeTree::IsAlways)
                    return false;
                break;
            case Value_OddInt:
                if(tree.GetEvennessInfo() != CodeTree::IsNever)
                    return false;
                break;
            case Value_IsInteger:
                if(!tree.IsAlwaysInteger(true)) return false;
                break;
            case Value_NonInteger:
                if(!tree.IsAlwaysInteger(false)) return false;
                break;
            case Value_Logical:
                if(!tree.IsLogicalValue()) return false;
                break;
        }
        switch(bitmask & SignMask)
        {
            case Sign_AnySign: /*case SignMask:*/ break;
            case Sign_Positive:
                if(!tree.IsAlwaysSigned(true)) return false;
                break;
            case Sign_Negative:
                if(!tree.IsAlwaysSigned(false)) return false;
                break;
            case Sign_NoIdea:
                if(tree.IsAlwaysSigned(true)) return false;
                if(tree.IsAlwaysSigned(false)) return false;
                break;
        }
        switch(bitmask & OnenessMask)
        {
            case Oneness_Any: case OnenessMask: break;
            case Oneness_One:
                if(!tree.IsImmed()) return false;
                if(!FloatEqual(fabs(tree.GetImmed()), 1.0)) return false;
                break;
            case Oneness_NotOne:
                if(!tree.IsImmed()) return false;
                if(FloatEqual(fabs(tree.GetImmed()), 1.0)) return false;
                break;
        }
        switch(bitmask & ConstnessMask)
        {
            case Constness_Any: /*case ConstnessMask:*/ break;
            case Constness_Const:
                if(!tree.IsImmed()) return false;
                break;
        }
        return true;
    }

    template<unsigned extent, unsigned nbits, typename item_type=unsigned int>
    struct nbitmap
    {
    private:
        static const unsigned bits_in_char = 8;
        static const unsigned per_item = (sizeof(item_type)*bits_in_char)/nbits;
        item_type data[(extent+per_item-1) / per_item];
    public:
        void inc(unsigned index, int by=1)
        {
            data[pos(index)] += by * item_type(1 << shift(index));
        }
        inline void dec(unsigned index) { inc(index, -1); }
        int get(unsigned index) const { return (data[pos(index)] >> shift(index)) & mask(); }

        static inline unsigned pos(unsigned index) { return index/per_item; }
        static inline unsigned shift(unsigned index) { return nbits * (index%per_item); }
        static inline unsigned mask() { return (1 << nbits)-1; }
        static inline unsigned mask(unsigned index) { return mask() << shift(index); }
    };

    struct Needs
    {
        int SubTrees     : 8; // This many subtrees
        int Others       : 8; // This many others (namedholder)
        int minimum_need : 8; // At least this many leaves (restholder may require more)
        int Immeds       : 8; // This many immeds

        nbitmap<VarBegin,2> SubTreesDetail; // This many subtrees of each opcode type

        Needs()
        {
            std::memset(this, 0, sizeof(*this));
        }
        Needs(const Needs& b)
        {
            std::memcpy(this, &b, sizeof(b));
        }
        Needs& operator= (const Needs& b)
        {
            std::memcpy(this, &b, sizeof(b));
            return *this;
        }
    };

    Needs CreateNeedList_uncached(const ParamSpec_SubFunctionData& params)
    {
        Needs NeedList;

        // Figure out what we need
        for(unsigned a = 0; a < params.param_count; ++a)
        {
            const ParamSpec& parampair = ParamSpec_Extract(params.param_list, a);
            switch(parampair.first)
            {
                case SubFunction:
                {
                    const ParamSpec_SubFunction& param = *(const ParamSpec_SubFunction*) parampair.second;
                    if(param.data.match_type == GroupFunction)
                        ++NeedList.Immeds;
                    else
                    {
                        ++NeedList.SubTrees;
                        assert( param.data.subfunc_opcode < VarBegin );
                        NeedList.SubTreesDetail.inc(param.data.subfunc_opcode);
                    }
                    ++NeedList.minimum_need;
                    break;
                }
                case NumConstant:
                case ParamHolder:
                    ++NeedList.Others;
                    ++NeedList.minimum_need;
                    break;
            }
        }

        return NeedList;
    }

    Needs& CreateNeedList(const ParamSpec_SubFunctionData& params)
    {
        typedef std::map<const ParamSpec_SubFunctionData*, Needs> needlist_cached_t;
        static needlist_cached_t needlist_cached;

        needlist_cached_t::iterator i = needlist_cached.lower_bound(&params);
        if(i != needlist_cached.end() && i->first == &params)
            return i->second;

        return
            needlist_cached.insert(i,
                 std::make_pair(&params, CreateNeedList_uncached(params))
            )->second;
    }
    /* Construct CodeTree from a GroupFunction, hopefully evaluating to a constant value */
    CodeTree CalculateGroupFunction(
        const ParamSpec& parampair,
        const MatchInfo& info)
    {
        switch( parampair.first )
        {
            case NumConstant:
            {
                const ParamSpec_NumConstant& param = *(const ParamSpec_NumConstant*) parampair.second;
                return CodeTree( param.constvalue ); // Note: calculates hash too.
            }
            case ParamHolder:
            {
                const ParamSpec_ParamHolder& param = *(const ParamSpec_ParamHolder*) parampair.second;
                return info.GetParamHolderValueIfFound( param.index );
                // If the ParamHolder is not defined, it will simply
                // return an Undefined tree. This is ok.
            }
            case SubFunction:
            {
                const ParamSpec_SubFunction& param = *(const ParamSpec_SubFunction*) parampair.second;
                /* Synthesize a CodeTree which will take care of
                 * constant-folding our expression. It will also
                 * indicate whether the result is, in fact,
                 * a constant at all. */
                CodeTree result;
                result.SetOpcode( param.data.subfunc_opcode );
                result.GetParams().reserve(param.data.param_count);
                for(unsigned a=0; a<param.data.param_count; ++a)
                {
                    CodeTree tmp(
                        CalculateGroupFunction
                        (ParamSpec_Extract(param.data.param_list, a), info)
                                );
                    result.AddParamMove(tmp);
                }
                result.Rehash(); // This will also call ConstantFolding().
                return result;
            }
        }
        // Issue an un-calculatable tree. (This should be unreachable)
        return CodeTree(); // cNop
    }
}

namespace FPoptimizer_Optimize
{
    /* Test the list of parameters to a given CodeTree */
    /* A helper function which simply checks whether the
     * basic shape of the tree matches what we are expecting
     * i.e. given number of numeric constants, etc.
     */
    bool IsLogisticallyPlausibleParamsMatch(
        const ParamSpec_SubFunctionData& params,
        const CodeTree& tree)
    {
        /* First, check if the tree has any chances of matching... */
        /* Figure out what we need. */
        Needs NeedList ( CreateNeedList(params) );

        size_t nparams = tree.GetParamCount();

        if(nparams < size_t(NeedList.minimum_need))
        {
            // Impossible to satisfy
            return false;
        }

        // Figure out what we have (note: we already assume that the opcode of the tree matches!)
        for(size_t a=0; a<nparams; ++a)
        {
            unsigned opcode = tree.GetParam(a).GetOpcode();
            switch(opcode)
            {
                case cImmed:
                    if(NeedList.Immeds > 0) --NeedList.Immeds;
                    else --NeedList.Others;
                    break;
                case VarBegin:
                case cFCall:
                case cPCall:
                    --NeedList.Others;
                    break;
                default:
                    assert( opcode < VarBegin );
                    if(NeedList.SubTrees > 0
                    && NeedList.SubTreesDetail.get(opcode) > 0)
                    {
                        --NeedList.SubTrees;
                        NeedList.SubTreesDetail.dec(opcode);
                    }
                    else --NeedList.Others;
            }
        }

        // Check whether all needs were satisfied
        if(NeedList.Immeds > 0
        || NeedList.SubTrees > 0
        || NeedList.Others > 0)
        {
            // Something came short, impossible to satisfy.
            return false;
        }

        if(params.match_type != AnyParams)
        {
            if(0
            //|| NeedList.Immeds < 0 - already checked
            || NeedList.SubTrees < 0
            || NeedList.Others < 0
            //|| params.count != nparams - already checked
              )
            {
                // Something was too much.
                return false;
            }
        }
        return true;
    }

    /* Test the given parameter to a given CodeTree */
    MatchResultType TestParam(
        const ParamSpec& parampair,
        const CodeTree& tree,
        const MatchPositionSpecBaseP& start_at,
        MatchInfo& info)
    {
        /*std::cout << "TestParam(";
        DumpParam(parampair);
        std::cout << ", ";
        DumpTree(tree);
        std::cout << ")\n";*/

        /* What kind of param are we expecting */
        switch( parampair.first )
        {
            case NumConstant: /* A particular numeric value */
            {
                const ParamSpec_NumConstant& param = *(const ParamSpec_NumConstant*) parampair.second;
                if(!tree.IsImmed()) return false;
                return FloatEqual(tree.GetImmed(), param.constvalue);
            }
            case ParamHolder: /* Any arbitrary node */
            {
                const ParamSpec_ParamHolder& param = *(const ParamSpec_ParamHolder*) parampair.second;
                if(!TestImmedConstraints(param.constraints, tree)) return false;
                return info.SaveOrTestParamHolder(param.index, tree);
            }
            case SubFunction:
            {
                const ParamSpec_SubFunction& param = *(const ParamSpec_SubFunction*) parampair.second;
                if(param.data.match_type == GroupFunction)
                { /* A constant value acquired from this formula */
                    if(!TestImmedConstraints(param.constraints, tree)) return false;
                    /* Construct the formula */
                    CodeTree  grammar_func = CalculateGroupFunction(parampair, info);
        #ifdef DEBUG_SUBSTITUTIONS
                    DumpHashes(grammar_func);
                    std::cout << *(const void**)&grammar_func.GetImmed();
                    std::cout << "\n";
                    std::cout << *(const void**)&tree.GetImmed();
                    std::cout << "\n";
                    DumpHashes(tree);
                    std::cout << "Comparing ";
                    DumpTree(grammar_func);
                    std::cout << " and ";
                    DumpTree(tree);
                    std::cout << ": ";
                    std::cout << (grammar_func.IsIdenticalTo(tree) ? "true" : "false");
                    std::cout << "\n";
        #endif
                    /* Evaluate it and compare */
                    return grammar_func.IsIdenticalTo(tree);
                }
                else /* A subtree conforming these specs */
                {
                    if(!&*start_at)
                    {
                        if(!TestImmedConstraints(param.constraints, tree)) return false;
                        if(tree.GetOpcode() != param.data.subfunc_opcode) return false;
                    }
                    return TestParams(param.data, tree, start_at, info, false);
                }
            }
        }
        return false;
    }

    struct PositionalParams_Rec
    {
        MatchPositionSpecBaseP start_at; /* child's start_at */
        MatchInfo              info;     /* backup of "info" at start */

        PositionalParams_Rec(): start_at(), info() { }
    };
    class MatchPositionSpec_PositionalParams
        : public MatchPositionSpecBase,
          public std::vector<PositionalParams_Rec>
    {
    public:
        explicit MatchPositionSpec_PositionalParams(size_t n)
            : MatchPositionSpecBase(),
              std::vector<PositionalParams_Rec> (n)
              { }
    };

    struct AnyWhere_Rec
    {
        MatchPositionSpecBaseP start_at; /* child's start_at */
        AnyWhere_Rec() : start_at() { }
    };
    class MatchPositionSpec_AnyWhere
        : public MatchPositionSpecBase,
          public std::vector<AnyWhere_Rec>
    {
    public:
        unsigned trypos;   /* which param index to try next */

        explicit MatchPositionSpec_AnyWhere(size_t n)
            : MatchPositionSpecBase(),
              std::vector<AnyWhere_Rec> (n),
              trypos(0)
              { }
    };

    MatchResultType TestParam_AnyWhere(
        const ParamSpec& parampair,
        const CodeTree& tree,
        const MatchPositionSpecBaseP& start_at,
        MatchInfo&         info,
        std::vector<bool>& used,
        bool TopLevel)
    {
        FPOPT_autoptr<MatchPositionSpec_AnyWhere> position;
        unsigned a;
        if(&*start_at)
        {
            position = (MatchPositionSpec_AnyWhere*) &*start_at;
            a = position->trypos;
            goto retry_anywhere_2;
        }
        else
        {
            position = new MatchPositionSpec_AnyWhere(tree.GetParamCount());
            a = 0;
        }
        for(; a < tree.GetParamCount(); ++a)
        {
            if(used[a]) continue;

        retry_anywhere:
          { MatchResultType r = TestParam(
                parampair,
                tree.GetParam(a),
                (*position)[a].start_at,
                info);

            (*position)[a].start_at = r.specs;
            if(r.found)
            {
                used[a]               = true; // matched
                if(TopLevel) info.SaveMatchedParamIndex(a);

                position->trypos = a; // in case of backtrack, try a again
                return MatchResultType(true, &*position);
            } }
        retry_anywhere_2:
            if(&*(*position)[a].start_at) // is there another try?
            {
                goto retry_anywhere;
            }
            // no, move on
        }
        return false;
    }

    struct AnyParams_Rec
    {
        MatchPositionSpecBaseP start_at; /* child's start_at */
        MatchInfo              info;     /* backup of "info" at start */
        std::vector<bool>      used;     /* which params are remaining */

        explicit AnyParams_Rec(size_t nparams)
            : start_at(), info(), used(nparams) { }
    };
    class MatchPositionSpec_AnyParams
        : public MatchPositionSpecBase,
          public std::vector<AnyParams_Rec>
    {
    public:
        explicit MatchPositionSpec_AnyParams(size_t n, size_t m)
            : MatchPositionSpecBase(),
              std::vector<AnyParams_Rec> (n, AnyParams_Rec(m))
              { }
    };

    /* Test the list of parameters to a given CodeTree */
    MatchResultType TestParams(
        const ParamSpec_SubFunctionData& model_tree,
        const CodeTree& tree,
        const MatchPositionSpecBaseP& start_at,
        MatchInfo& info,
        bool TopLevel)
    {
        /* When PositionalParams or SelectedParams, verify that
         * the number of parameters is exactly as expected.
         */
        if(model_tree.match_type != AnyParams)
        {
            if(model_tree.param_count != tree.GetParamCount())
                return false;
        }

        /* Verify that the tree basically conforms the shape we are expecting */
        /* This test is not necessary; it may just save us some work. */
        if(!IsLogisticallyPlausibleParamsMatch(model_tree, tree))
        {
            return false;
        }

        /* Verify each parameter that they are found in the tree as expected. */
        switch(model_tree.match_type)
        {
            case PositionalParams:
            {
                /* Simple: Test all given parameters in succession. */
                FPOPT_autoptr<MatchPositionSpec_PositionalParams> position;
                unsigned a;
                if(&*start_at)
                {
                    position = (MatchPositionSpec_PositionalParams*) &*start_at;
                    a = model_tree.param_count - 1;
                    goto retry_positionalparams_2;
                }
                else
                {
                    position = new MatchPositionSpec_PositionalParams(model_tree.param_count);
                    a = 0;
                }

                for(; a < model_tree.param_count; ++a)
                {
                    (*position)[a].info = info;
                retry_positionalparams:
                  { MatchResultType r = TestParam(
                        ParamSpec_Extract(model_tree.param_list, a),
                        tree.GetParam(a),
                        (*position)[a].start_at,
                        info);

                    (*position)[a].start_at = r.specs;
                    if(r.found)
                    {
                        continue;
                  } }
                retry_positionalparams_2:
                    // doesn't match
                    if(&*(*position)[a].start_at) // is there another try?
                    {
                        info = (*position)[a].info;
                        goto retry_positionalparams;
                    }
                    // no, backtrack
                    if(a > 0)
                    {
                        --a;
                        goto retry_positionalparams_2;
                    }
                    // cannot backtrack
                    info = (*position)[0].info;
                    return false;
                }
                if(TopLevel)
                    for(unsigned a = 0; a < model_tree.param_count; ++a)
                        info.SaveMatchedParamIndex(a);
                return MatchResultType(true, &*position);
            }
            case SelectedParams:
                // same as AnyParams, except that model_tree.count==tree.GetParamCount()
                //                       and that there are no RestHolders
            case AnyParams:
            {
                /* Ensure that all given parameters are found somewhere, in any order */

                FPOPT_autoptr<MatchPositionSpec_AnyParams> position;
                std::vector<bool> used( tree.GetParamCount() );
                std::vector<unsigned> depcodes( model_tree.param_count );
                std::vector<unsigned> test_order( model_tree.param_count );
                for(unsigned a=0; a<model_tree.param_count; ++a)
                {
                    const ParamSpec parampair = ParamSpec_Extract(model_tree.param_list, a);
                    depcodes[a] = ParamSpec_GetDepCode(parampair);
                }
                { unsigned b=0;
                for(unsigned a=0; a<model_tree.param_count; ++a)
                    if(depcodes[a] != 0)
                        test_order[b++] = a;
                for(unsigned a=0; a<model_tree.param_count; ++a)
                    if(depcodes[a] == 0)
                        test_order[b++] = a;
                }

                unsigned a;
                if(&*start_at)
                {
                    position = (MatchPositionSpec_AnyParams*) &*start_at;
                    if(model_tree.param_count == 0)
                    {
                        a = 0;
                        goto retry_anyparams_4;
                    }
                    a = model_tree.param_count - 1;
                    goto retry_anyparams_2;
                }
                else
                {
                    position = new MatchPositionSpec_AnyParams(model_tree.param_count,
                                                               tree.GetParamCount());
                    a = 0;
                    if(model_tree.param_count != 0)
                    {
                        (*position)[0].info   = info;
                        (*position)[0].used   = used;
                    }
                }
                // Match all but restholders
                for(; a < model_tree.param_count; ++a)
                {
                    if(a > 0) // this test is not necessary, but it saves from doing
                    {         // duplicate work, because [0] was already saved above.
                        (*position)[a].info   = info;
                        (*position)[a].used   = used;
                    }
                retry_anyparams:
                  { MatchResultType r = TestParam_AnyWhere(
                        ParamSpec_Extract(model_tree.param_list, test_order[a]),
                        tree,
                        (*position)[a].start_at,
                        info,
                        used,
                        TopLevel);
                    (*position)[a].start_at = r.specs;
                    if(r.found)
                    {
                        continue;
                  } }
                retry_anyparams_2:
                    // doesn't match
                    if(&*(*position)[a].start_at) // is there another try?
                    {
                        info = (*position)[a].info;
                        used = (*position)[a].used;
                        goto retry_anyparams;
                    }
                    // no, backtrack
                retry_anyparams_3:
                    if(a > 0)
                    {
                        --a;
                        goto retry_anyparams_2;
                    }
                    // cannot backtrack
                    info = (*position)[0].info;
                    return false;
                }
            retry_anyparams_4:
                // Capture anything remaining in the restholder
                if(model_tree.restholder_index != 0)
                {
                    //std::vector<bool> used_backup(used);
                    //MatchInfo         info_backup(info);

                    if(!TopLevel
                    || !info.HasRestHolder(model_tree.restholder_index))
                    {
                        std::vector<CodeTree> matches;
                        matches.reserve(tree.GetParamCount());
                        for(unsigned b = 0; b < tree.GetParamCount(); ++b)
                        {
                            if(used[b]) continue; // Ignore subtrees that were already used
                            // Save this tree to this restholder

                            matches.push_back(tree.GetParam(b));
                            used[b] = true;
                            if(TopLevel) info.SaveMatchedParamIndex(b);
                        }
                        if(!info.SaveOrTestRestHolder(model_tree.restholder_index, matches))
                        {
                            // Failure at restholder matching. Backtrack if possible.
                            //used.swap(used_backup);
                            //info.swap(info_backup);
                            goto retry_anyparams_3;
                        }
                        //std::cout << "Saved restholder " << model_tree.restholder_index << "\n";
                    }
                    else
                    {
                        const std::vector<CodeTree>& matches
                            = info.GetRestHolderValues(model_tree.restholder_index);
                        //std::cout << "Testing restholder " << model_tree.restholder_index << std::flush;
                        for(size_t a=0; a<matches.size(); ++a)
                        {
                            bool found = false;
                            for(unsigned b = 0; b < tree.GetParamCount(); ++b)
                            {
                                if(used[b]) continue;
                                if(matches[a].IsIdenticalTo(tree.GetParam(b)))
                                {
                                    used[b] = true;
                                    if(TopLevel) info.SaveMatchedParamIndex(b);
                                    found = true;
                                    break;
                                }
                            }
                            if(!found)
                            {
                                //std::cout << " ... failed\n";
                                // Failure at restholder matching. Backtrack if possible.
                                //used.swap(used_backup);
                                //info.swap(info_backup);
                                goto retry_anyparams_3;
                            }
                        }
                        //std::cout << " ... ok\n";
                    }
                }
                return MatchResultType(true, model_tree.param_count ? &*position : 0);
            }
            case GroupFunction: // never occurs
                break;
        }
        return false; // doesn't match
    }
}

#endif

#line 1 "fpoptimizer/fpoptimizer_optimize_synth.c"
#include "fpconfig.h"
#include "fparser.h"
#include "fptypes.h"

#ifdef FP_SUPPORT_OPTIMIZER

#include <algorithm>
#include <assert.h>

// line removed for fpoptimizer.c: #include "fpoptimizer_optimize.h"

namespace FPoptimizer_Optimize
{
    /* Synthesize the given grammatic parameter into the codetree */
    void SynthesizeParam(
        const ParamSpec& parampair,
        CodeTree& tree,
        MatchInfo& info,
        bool inner = true)
    {
        switch( parampair.first )
        {
            case NumConstant:
              { const ParamSpec_NumConstant& param = *(const ParamSpec_NumConstant*) parampair.second;
                tree.SetImmed( param.constvalue );
                if(inner) tree.Rehash(false);
                break; }
            case ParamHolder:
              { const ParamSpec_ParamHolder& param = *(const ParamSpec_ParamHolder*) parampair.second;
                tree.Become( info.GetParamHolderValue( param.index ) );
                break; }
            case SubFunction:
              { const ParamSpec_SubFunction& param = *(const ParamSpec_SubFunction*) parampair.second;
                tree.SetOpcode( param.data.subfunc_opcode );
                for(unsigned a=0; a < param.data.param_count; ++a)
                {
                    CodeTree nparam;
                    SynthesizeParam( ParamSpec_Extract(param.data.param_list, a), nparam, info, true );
                    tree.AddParamMove(nparam);
                }
                if(param.data.restholder_index != 0)
                {
                    std::vector<CodeTree> trees
                        ( info.GetRestHolderValues( param.data.restholder_index ) );
                    tree.AddParamsMove(trees);
                    // ^note: this fails if the same restholder is synth'd twice
                    if(tree.GetParamCount() == 1)
                    {
                        /* Convert cMul <1> into <1> when <1> only contains one operand.
                         * This is redundant code; it is also done in ConstantFolding(),
                         * but it might be better for performance to do it here, too.
                         */
                        assert(tree.GetOpcode() == cAdd || tree.GetOpcode() == cMul
                            || tree.GetOpcode() == cMin || tree.GetOpcode() == cMax
                            || tree.GetOpcode() == cAnd || tree.GetOpcode() == cOr);
                        tree.Become(tree.GetParam(0));
                    }
                }
                if(inner)
                    tree.Rehash();
                break; }
        }
    }

    void SynthesizeRule(
        const Rule& rule,
        CodeTree& tree,
        MatchInfo& info)
    {
        switch(rule.ruletype)
        {
            case ProduceNewTree:
            {
                tree.DelParams();
                SynthesizeParam( ParamSpec_Extract(rule.repl_param_list, 0), tree, info, false );
                break;
            }
            case ReplaceParams:
            default:
            {
                /* Delete the matched parameters from the source tree */
                std::vector<unsigned> list = info.GetMatchedParamIndexes();
                std::sort(list.begin(), list.end());
                for(size_t a=list.size(); a-->0; )
                    tree.DelParam( list[a] );

                /* Synthesize the replacement params */
                for(unsigned a=0; a < rule.repl_param_count; ++a)
                {
                    CodeTree nparam;
                    SynthesizeParam( ParamSpec_Extract(rule.repl_param_list, a), nparam, info, true );
                    tree.AddParamMove(nparam);
                }
                break;
            }
        }
    }
}

#endif

#line 1 "fpoptimizer/fpoptimizer_optimize_debug.c"
// line removed for fpoptimizer.c: #include "fpoptimizer_grammar.h"
// line removed for fpoptimizer.c: #include "fpoptimizer_opcodename.h"
// line removed for fpoptimizer.c: #include "fpoptimizer_optimize.h"

#ifdef DEBUG_SUBSTITUTIONS

#include <sstream>
#include <cstring>

using namespace FUNCTIONPARSERTYPES;
using namespace FPoptimizer_Grammar;
using namespace FPoptimizer_CodeTree;
using namespace FPoptimizer_Optimize;

namespace FPoptimizer_Grammar
{
    void DumpMatch(const Rule& rule,
                   const CodeTree& tree,
                   const MatchInfo& info,
                   bool DidMatch,
                   std::ostream& o)
    {
        DumpMatch(rule,tree,info,DidMatch?"Found match":"Found mismatch",o);
    }

    void DumpMatch(const Rule& rule,
                   const CodeTree& tree,
                   const MatchInfo& info,
                   const char* whydump,
                   std::ostream& o)
    {
        static const char ParamHolderNames[][2] = {"%","&","x","y","z","a","b","c"};

        o << whydump
          << " (rule " << (&rule - grammar_rules) << ")"
          << ":\n"
            "  Pattern    : ";
        { ParamSpec tmp;
          tmp.first = SubFunction;
          ParamSpec_SubFunction tmp2;
          tmp2.data = rule.match_tree;
          tmp.second = (const void*) &tmp2;
          DumpParam(tmp, o);
        }
        o << "\n"
            "  Replacement: ";
        DumpParams(rule.repl_param_list, rule.repl_param_count, o);
        o << "\n";

        o <<
            "  Tree       : ";
        DumpTree(tree, o);
        o << "\n";
        if(!std::strcmp(whydump,"Found match")) DumpHashes(tree, o);

        for(size_t a=0; a<info.paramholder_matches.size(); ++a)
        {
            if(!info.paramholder_matches[a].IsDefined()) continue;
            o << "           " << ParamHolderNames[a] << " = ";
            DumpTree(info.paramholder_matches[a], o);
            o << "\n";
        }

        for(size_t b=0; b<info.restholder_matches.size(); ++b)
        {
            if(!info.restholder_matches[b].first) continue;
            for(size_t a=0; a<info.restholder_matches[b].second.size(); ++a)
            {
                o << "         <" << b << "> = ";
                DumpTree(info.restholder_matches[b].second[a], o);
                o << std::endl;
            }
        }
        o << std::flush;
    }

}

#endif

#line 1 "fpoptimizer/fpoptimizer_hash.c"
#include <list>
#include <algorithm>

// line removed for fpoptimizer.c: #include "fpoptimizer_codetree.h"
#include "fptypes.h"
// line removed for fpoptimizer.c: #include "crc32.h"

#ifdef FP_SUPPORT_OPTIMIZER

using namespace FUNCTIONPARSERTYPES;
//using namespace FPoptimizer_Grammar;

namespace
{
    bool MarkIncompletes(FPoptimizer_CodeTree::CodeTree& tree)
    {
        if(tree.Is_Incompletely_Hashed())
            return true;

        bool needs_rehash = false;
        for(size_t a=0; a<tree.GetParamCount(); ++a)
            needs_rehash |= MarkIncompletes(tree.GetParam(a));
        if(needs_rehash)
            tree.Mark_Incompletely_Hashed();
        return needs_rehash;
    }

    void FixIncompletes(FPoptimizer_CodeTree::CodeTree& tree)
    {
        if(tree.Is_Incompletely_Hashed())
        {
            for(size_t a=0; a<tree.GetParamCount(); ++a)
                FixIncompletes(tree.GetParam(a));
            tree.Rehash();
        }
    }
}

namespace FPoptimizer_CodeTree
{
    void CodeTree::Rehash(bool constantfolding)
    {
        if(constantfolding)
            ConstantFolding();
        data->Sort();
        data->Recalculate_Hash_NoRecursion();
    }

    void CodeTreeData::Recalculate_Hash_NoRecursion()
    {
        fphash_t NewHash = { Opcode * FPHASH_CONST(0x3A83A83A83A83A0),
                             Opcode * FPHASH_CONST(0x1131462E270012B)};
        Depth = 1;
        switch(Opcode)
        {
            case cImmed:
                if(Value != 0.0)
                {
                    crc32_t crc = crc32::calc( (const unsigned char*) &Value,
                                                sizeof(Value) );
                    NewHash.hash1 ^= crc | (fphash_value_t(crc) << FPHASH_CONST(32));
                    NewHash.hash2 += ((~fphash_value_t(crc)) * 3) ^ 1234567;
                }
                break; // no params
            case VarBegin:
                NewHash.hash1 ^= (Var<<24) | (Var>>24);
                NewHash.hash2 += (fphash_value_t(Var)*5) ^ 2345678;
                break; // no params
            case cFCall: case cPCall:
            {
                crc32_t crc = crc32::calc( (const unsigned char*) &Funcno, sizeof(Funcno) );
                NewHash.hash1 ^= (crc<<24) | (crc>>24);
                NewHash.hash2 += ((~fphash_value_t(crc)) * 7) ^ 3456789;
                /* passthru */
            }
            default:
            {
                size_t MaxChildDepth = 0;
                for(size_t a=0; a<Params.size(); ++a)
                {
                    if(Params[a].GetDepth() > MaxChildDepth)
                        MaxChildDepth = Params[a].GetDepth();

                    NewHash.hash1 += (1)*FPHASH_CONST(0x2492492492492492);
                    NewHash.hash1 *= FPHASH_CONST(1099511628211);
                    //assert(&*Params[a] != this);
                    NewHash.hash1 += Params[a].GetHash().hash1;

                    NewHash.hash2 += (3)*FPHASH_CONST(0x9ABCD801357);
                    NewHash.hash2 *= FPHASH_CONST(0xECADB912345);
                    NewHash.hash2 += (~Params[a].GetHash().hash1) ^ 4567890;
                }
                Depth += MaxChildDepth;
            }
        }
        if(Hash != NewHash)
        {
            Hash = NewHash;
            OptimizedUsing = 0;
        }
    }

    void CodeTree::FixIncompleteHashes()
    {
        MarkIncompletes(*this);
        FixIncompletes(*this);
    }
}

#endif

#line 1 "fpoptimizer/fpoptimizer_makebytecode.c"
#include <cmath>
#include <list>
#include <cassert>

// line removed for fpoptimizer.c: #include "fpoptimizer_codetree.h"
#include "fptypes.h"
// line removed for fpoptimizer.c: #include "fpoptimizer_consts.h"
// line removed for fpoptimizer.c: #include "fpoptimizer_bytecodesynth.h"
// line removed for fpoptimizer.c: #include "fpoptimizer_optimize.h"

#ifdef FP_SUPPORT_OPTIMIZER

using namespace FUNCTIONPARSERTYPES;
//using namespace FPoptimizer_Grammar;

namespace
{
    using namespace FPoptimizer_CodeTree;

    bool AssembleSequence(
                  const CodeTree& tree, long count,
                  const FPoptimizer_ByteCode::SequenceOpCode& sequencing,
                  FPoptimizer_ByteCode::ByteCodeSynth& synth,
                  size_t max_bytecode_grow_length);
}

namespace FPoptimizer_CodeTree
{
    void CodeTree::SynthesizeByteCode(
        std::vector<unsigned>& ByteCode,
        std::vector<double>&   Immed,
        size_t& stacktop_max)
    {
    #ifdef DEBUG_SUBSTITUTIONS
        std::cout << "Making bytecode for:\n";
        DumpTreeWithIndent(*this);
    #endif
        while(RecreateInversionsAndNegations())
        {
        #ifdef DEBUG_SUBSTITUTIONS
            std::cout << "One change issued, produced:\n";
            DumpTreeWithIndent(*this);
        #endif
            FixIncompleteHashes();
        }
    #ifdef DEBUG_SUBSTITUTIONS
        std::cout << "Actually synthesizing, after recreating inv/neg:\n";
        DumpTreeWithIndent(*this);
    #endif

        FPoptimizer_ByteCode::ByteCodeSynth synth;

        /* Then synthesize the actual expression */
        SynthesizeByteCode(synth, false);
        /* The "false" parameters tells SynthesizeByteCode
         * that at the outermost synthesizing level, it does
         * not matter if leftover temps are left in the stack.
         */
        synth.Pull(ByteCode, Immed, stacktop_max);
    }

    void CodeTree::SynthesizeByteCode(
        FPoptimizer_ByteCode::ByteCodeSynth& synth,
        bool MustPopTemps) const
    {
        // If the synth can already locate our operand in the stack,
        // never mind synthesizing it again, just dup it.
        if(synth.FindAndDup(*this))
        {
            return;
        }

        size_t n_subexpressions_synthesized = SynthCommonSubExpressions(synth);

        switch(GetOpcode())
        {
            case VarBegin:
                synth.PushVar(GetVar());
                break;
            case cImmed:
                synth.PushImmed(GetImmed());
                break;
            case cAdd:
            case cMul:
            case cMin:
            case cMax:
            case cAnd:
            case cOr:
            case cAbsAnd:
            case cAbsOr:
            {
                if(GetOpcode() == cMul) // Special treatment for cMul sequences
                {
                    // If the paramlist contains an Immed, and that Immed
                    // fits in a long-integer, try to synthesize it
                    // as add-sequences instead.
                    bool did_muli = false;
                    for(size_t a=0; a<GetParamCount(); ++a)
                    {
                        if(GetParam(a).IsLongIntegerImmed())
                        {
                            long value = GetParam(a).GetLongIntegerImmed();

                            CodeTree tmp(*this, CodeTree::CloneTag());
                            tmp.DelParam(a);
                            tmp.Rehash();
                            if(AssembleSequence(
                                tmp, value, FPoptimizer_ByteCode::AddSequence,
                                synth,
                                MAX_MULI_BYTECODE_LENGTH))
                            {
                                did_muli = true;
                                break;
                            }
                        }
                    }
                    if(did_muli)
                        break; // done
                }

                int n_stacked = 0;
                for(size_t a=0; a<GetParamCount(); ++a)
                {
                    GetParam(a).SynthesizeByteCode(synth);
                    ++n_stacked;

                    if(n_stacked > 1)
                    {
                        // Cumulate at the earliest opportunity.
                        synth.AddOperation(GetOpcode(), 2); // stack state: -2+1 = -1
                        n_stacked = n_stacked - 2 + 1;
                    }
                }
                if(n_stacked == 0)
                {
                    // Uh, we got an empty cAdd/cMul/whatever...
                    // Synthesize a default value.
                    // This should never happen.
                    switch(GetOpcode())
                    {
                        case cAdd:
                        case cOr:
                        case cAbsOr:
                            synth.PushImmed(0);
                            break;
                        case cMul:
                        case cAnd:
                        case cAbsAnd:
                            synth.PushImmed(1);
                            break;
                        case cMin:
                        case cMax:
                            //synth.PushImmed(NaN);
                            synth.PushImmed(0);
                            break;
                        default:
                            break;
                    }
                    ++n_stacked;
                }
                assert(n_stacked == 1);
                break;
            }
            case cPow:
            {
                const CodeTree& p0 = GetParam(0);
                const CodeTree& p1 = GetParam(1);

                if(!p1.IsLongIntegerImmed()
                || !AssembleSequence( /* Optimize integer exponents */
                        p0, p1.GetLongIntegerImmed(),
                        FPoptimizer_ByteCode::MulSequence,
                        synth,
                        MAX_POWI_BYTECODE_LENGTH)
                  )
                {
                    p0.SynthesizeByteCode(synth);
                    p1.SynthesizeByteCode(synth);
                    synth.AddOperation(GetOpcode(), 2); // Create a vanilla cPow.
                }
                break;
            }
            case cIf:
            case cAbsIf:
            {
                // Assume that the parameter count is 3 as it should.
                FPoptimizer_ByteCode::ByteCodeSynth::IfData ifdata;

                GetParam(0).SynthesizeByteCode(synth); // expression

                synth.SynthIfStep1(ifdata, GetOpcode());

                GetParam(1).SynthesizeByteCode(synth); // true branch

                synth.SynthIfStep2(ifdata);

                GetParam(2).SynthesizeByteCode(synth); // false branch

                synth.SynthIfStep3(ifdata);
                break;
            }
            case cFCall:
            case cPCall:
            {
                // Assume that the parameter count is as it should.
                for(size_t a=0; a<GetParamCount(); ++a)
                    GetParam(a).SynthesizeByteCode(synth);
                synth.AddOperation(GetOpcode(), (unsigned) GetParamCount());
                synth.AddOperation(GetFuncNo(), 0, 0);
                break;
            }
            default:
            {
                // Assume that the parameter count is as it should.
                for(size_t a=0; a<GetParamCount(); ++a)
                    GetParam(a).SynthesizeByteCode(synth);
                synth.AddOperation(GetOpcode(), (unsigned) GetParamCount());
                break;
            }
        }

        // Tell the synthesizer which tree was just produced in the stack
        synth.StackTopIs(*this);

        // If we added subexpressions, peel them off the stack now
        if(MustPopTemps && n_subexpressions_synthesized > 0)
        {
            size_t top = synth.GetStackTop();
            synth.DoPopNMov(top-1-n_subexpressions_synthesized, top-1);
        }
    }
}

namespace
{
    bool AssembleSequence(
        const CodeTree& tree, long count,
        const FPoptimizer_ByteCode::SequenceOpCode& sequencing,
        FPoptimizer_ByteCode::ByteCodeSynth& synth,
        size_t max_bytecode_grow_length)
    {
        if(count != 0)
        {
            FPoptimizer_ByteCode::ByteCodeSynth backup = synth;

            tree.SynthesizeByteCode(synth);

            // Ignore the size generated by subtree
            size_t bytecodesize_backup = synth.GetByteCodeSize();

            FPoptimizer_ByteCode::AssembleSequence(count, sequencing, synth);

            size_t bytecode_grow_amount = synth.GetByteCodeSize() - bytecodesize_backup;
            if(bytecode_grow_amount > max_bytecode_grow_length)
            {
                synth = backup;
                return false;
            }
            return true;
        }
        else
        {
            FPoptimizer_ByteCode::AssembleSequence(count, sequencing, synth);
            return true;
        }
    }
}

#endif

#line 1 "fpoptimizer/fpoptimizer_readbytecode.c"
#include <cmath>
#include <cassert>

// line removed for fpoptimizer.c: #include "fpoptimizer_codetree.h"
// line removed for fpoptimizer.c: #include "fpoptimizer_optimize.h"
// line removed for fpoptimizer.c: #include "fpoptimizer_opcodename.h"
// line removed for fpoptimizer.c: #include "fpoptimizer_grammar.h"
#include "fptypes.h"

// line removed for fpoptimizer.c: #include "fpoptimizer_consts.h"
#include "fparser.h"

#ifdef FP_SUPPORT_OPTIMIZER

using namespace FUNCTIONPARSERTYPES;
//using namespace FPoptimizer_Grammar;

namespace
{
    using namespace FPoptimizer_CodeTree;

    typedef std::vector<double> FactorStack;

    const struct PowiMuliType
    {
        unsigned opcode_square;
        unsigned opcode_cumulate;
        unsigned opcode_invert;
        unsigned opcode_half;
        unsigned opcode_invhalf;
    } iseq_powi = {cSqr,cMul,cInv,cSqrt,cRSqrt},
      iseq_muli = {~unsigned(0), cAdd,cNeg, ~unsigned(0),~unsigned(0) };

    double ParsePowiMuli(
        const PowiMuliType& opcodes,
        const std::vector<unsigned>& ByteCode, size_t& IP,
        size_t limit,
        size_t factor_stack_base,
        FactorStack& stack)
    {
        double result = 1;
        while(IP < limit)
        {
            if(ByteCode[IP] == opcodes.opcode_square)
            {
                if(!IsIntegerConst(result)) break;
                result *= 2;
                ++IP;
                continue;
            }
            if(ByteCode[IP] == opcodes.opcode_invert)
            {
                result = -result;
                ++IP;
                continue;
            }
            if(ByteCode[IP] == opcodes.opcode_half)
            {
                if(IsIntegerConst(result) && result > 0 && ((long)result) % 2 == 0)
                    break;
                result *= 0.5;
                ++IP;
                continue;
            }
            if(ByteCode[IP] == opcodes.opcode_invhalf)
            {
                if(IsIntegerConst(result) && result > 0 && ((long)result) % 2 == 0)
                    break;
                result *= -0.5;
                ++IP;
                continue;
            }

            size_t dup_fetch_pos = IP;
            double lhs = 1.0;

            if(ByteCode[IP] == cFetch)
            {
                unsigned index = ByteCode[++IP];
                if(index < factor_stack_base
                || size_t(index-factor_stack_base) >= stack.size())
                {
                    // It wasn't a powi-fetch after all
                    IP = dup_fetch_pos;
                    break;
                }
                lhs = stack[index - factor_stack_base];
                // Note: ^This assumes that cFetch of recentmost
                //        is always converted into cDup.
                goto dup_or_fetch;
            }
            if(ByteCode[IP] == cDup)
            {
                lhs = result;
                goto dup_or_fetch;

            dup_or_fetch:
                stack.push_back(result);
                ++IP;
                double subexponent = ParsePowiMuli
                    (opcodes,
                     ByteCode, IP, limit,
                     factor_stack_base, stack);
                if(IP >= limit || ByteCode[IP] != opcodes.opcode_cumulate)
                {
                    // It wasn't a powi-dup after all
                    IP = dup_fetch_pos;
                    break;
                }
                ++IP; // skip opcode_cumulate
                stack.pop_back();
                result += lhs*subexponent;
                continue;
            }
            break;
        }
        return result;
    }

    double ParsePowiSequence(const std::vector<unsigned>& ByteCode, size_t& IP,
                           size_t limit,
                           size_t factor_stack_base)
    {
        FactorStack stack;
        stack.push_back(1.0);
        return ParsePowiMuli(iseq_powi, ByteCode, IP, limit, factor_stack_base, stack);
    }

    double ParseMuliSequence(const std::vector<unsigned>& ByteCode, size_t& IP,
                           size_t limit,
                           size_t factor_stack_base)
    {
        FactorStack stack;
        stack.push_back(1.0);
        return ParsePowiMuli(iseq_muli, ByteCode, IP, limit, factor_stack_base, stack);
    }

    class CodeTreeParserData
    {
    public:
        explicit CodeTreeParserData(bool k_powi)
            : stack(), keep_powi(k_powi) { }

        void Eat(size_t nparams, OPCODE opcode)
        {
            CodeTree newnode;
            newnode.SetOpcode(opcode);

            std::vector<CodeTree> params = Pop(nparams);
            newnode.SetParamsMove(params);

            if(!keep_powi)
            switch(opcode)
            {
                //        asinh: log(x + sqrt(x*x + 1))
                //cAsinh [x] -> cLog (cAdd x (cPow (cAdd (cPow x 2) 1) 0.5))
                // Note: ^ Replacement function refers to x twice

                //        acosh: log(x + sqrt(x*x - 1))
                //cAcosh [x] -> cLog (cAdd x (cPow (cAdd (cPow x 2) -1) 0.5))

                //        atanh: log( (1+x) / (1-x)) / 2
                //cAtanh [x] -> cMul (cLog (cMul (cAdd 1 x) (cPow (cAdd 1 (cMul -1 x)) -1))) 0.5

                //        asin: atan2(x, sqrt(1-x*x))
                //cAsin[x] -> cAtan2 [x (cPow [(cAdd 1 (cMul (cPow [x 2] -1)) 0.5])]

                //        acos: atan2(sqrt(1-x*x), x)
                //cAcos[x] -> cAtan2 [(cPow [(cAdd 1 (cMul (cPow [x 2] -1)) 0.5]) x]

                //     The hyperbolic functions themselves are:
                //        sinh: (exp(x)-exp(-x)) / 2  = exp(-x) * (exp(2*x)-1) / 2
                //cSinh [x] -> cMul 0.5 (cPow [CONSTANT_EI x]) (cAdd [-1 (cPow [CONSTANT_2E x])])

                //        cosh: (exp(x)+exp(-x)) / 2  = exp(-x) * (exp(2*x)+1) / 2
                //        cosh(-x) = cosh(x)
                //cCosh [x] -> cMul 0.5 (cPow [CONSTANT_EI x]) (cAdd [ 1 (cPow [CONSTANT_2E x])])

                //        tanh: sinh/cosh = (exp(2*x)-1) / (exp(2*x)+1)
                //cTanh [x] -> (cMul (cAdd {(cPow [CONSTANT_2E x]) -1}) (cPow [(cAdd {(cPow [CONSTANT_2E x]) 1}) -1]))
                case cTanh:
                {
                    CodeTree sinh, cosh;
                    sinh.SetOpcode(cSinh); sinh.AddParam(newnode.GetParam(0)); sinh.Rehash();
                    cosh.SetOpcode(cCosh); cosh.AddParamMove(newnode.GetParam(0)); cosh.Rehash();
                    CodeTree pow;
                    pow.SetOpcode(cPow);
                    pow.AddParamMove(cosh);
                    pow.AddParam(CodeTree(-1.0));
                    pow.Rehash();
                    newnode.SetOpcode(cMul);
                    newnode.SetParamMove(0, sinh);
                    newnode.AddParamMove(pow);
                    break;
                }

                //        tan: sin/cos
                //cTan [x] -> (cMul (cSin [x]) (cPow [(cCos [x]) -1]))
                case cTan:
                {
                    CodeTree sin, cos;
                    sin.SetOpcode(cSin); sin.AddParam(newnode.GetParam(0)); sin.Rehash();
                    cos.SetOpcode(cCos); cos.AddParamMove(newnode.GetParam(0)); cos.Rehash();
                    CodeTree pow;
                    pow.SetOpcode(cPow);
                    pow.AddParamMove(cos);
                    pow.AddParam(CodeTree(-1.0));
                    pow.Rehash();
                    newnode.SetOpcode(cMul);
                    newnode.SetParamMove(0, sin);
                    newnode.AddParamMove(pow);
                    break;
                }

                case cPow:
                {
                    const CodeTree& p0 = newnode.GetParam(0);
                    const CodeTree& p1 = newnode.GetParam(1);
                    if(p1.GetOpcode() == cAdd)
                    {
                        // convert x^(a + b) into x^a * x^b just so that
                        // some optimizations can be run on it.
                        // For instance, exp(log(x)*-61.1 + log(z)*-59.1)
                        // won't be changed into exp(log(x*z)*-61.1)*z^2
                        // unless we do this.
                        std::vector<CodeTree> mulgroup(p1.GetParamCount());
                        for(size_t a=0; a<p1.GetParamCount(); ++a)
                        {
                            CodeTree pow;
                            pow.SetOpcode(cPow);
                            pow.AddParam(p0);
                            pow.AddParam(p1.GetParam(a));
                            pow.Rehash();
                            mulgroup[a].swap(pow);
                        }
                        newnode.SetOpcode(cMul);
                        newnode.SetParamsMove(mulgroup);
                    }
                    break;
                }

                // Should we change sin(x) into cos(pi/2-x)
                //               or cos(x) into sin(pi/2-x)?
                //                        note: cos(x-pi/2) = cos(pi/2-x) = sin(x)
                //                        note: sin(x-pi/2) = -sin(pi/2-x) = -cos(x)
                default: break;
            }

            newnode.Rehash(!keep_powi);
        /*
            using namespace FPoptimizer_Grammar;
            bool recurse = false;
            while(ApplyGrammar(pack.glist[0], newnode, recurse)) // intermediate
            { //std::cout << "Rerunning 1\n";
                newnode.FixIncompleteHashes();
                recurse = true;
            }
        */
            FindClone(newnode, false);
        #ifdef DEBUG_SUBSTITUTIONS
            std::cout << "POP " << nparams << ", " << FP_GetOpcodeName(opcode)
                      << "->" << FP_GetOpcodeName(newnode.GetOpcode())
                      << ": PUSH ";
            DumpTree(newnode);
            std::cout <<std::endl;
        #endif
            stack.push_back(newnode);
        }

        void EatFunc(size_t nparams, OPCODE opcode, unsigned funcno)
        {
            CodeTree newnode;
            newnode.SetFuncOpcode(opcode, funcno);
            std::vector<CodeTree> params = Pop(nparams);
            newnode.SetParamsMove(params);
            newnode.Rehash(false);
        #ifdef DEBUG_SUBSTITUTIONS
            std::cout << "POP " << nparams << ", PUSH ";
            DumpTree(newnode);
            std::cout << std::endl;
        #endif
            FindClone(newnode);
            stack.push_back(newnode);
        }

        void AddConst(double value)
        {
            CodeTree newnode(value);
            FindClone(newnode);
            Push(newnode);
        }

        void AddVar(unsigned varno)
        {
            CodeTree newnode(varno, CodeTree::VarTag());
            FindClone(newnode);
            Push(newnode);
        }

        void SwapLastTwoInStack()
        {
            stack[stack.size()-1].swap( stack[stack.size()-2] );
        }

        void Dup()
        {
            Fetch(stack.size()-1);
        }

        void Fetch(size_t which)
        {
            Push(stack[which]);
        }

        template<typename T>
        void Push(T tree)
        {
        #ifdef DEBUG_SUBSTITUTIONS
            std::cout << "PUSH ";
            DumpTree(tree);
            std::cout << std::endl;
        #endif
            stack.push_back(tree);
        }

        void PopNMov(size_t target, size_t source)
        {
            stack[target] = stack[source];
            stack.resize(target+1);
        }

        CodeTree PullResult()
        {
            clones.clear();
            CodeTree result(stack.back());
            stack.resize(stack.size()-1);
            return result;
        }
        std::vector<CodeTree> Pop(unsigned n_pop)
        {
            std::vector<CodeTree> result(n_pop);
            for(unsigned n=0; n<n_pop; ++n)
                result[n].swap(stack[stack.size()-n_pop+n]);
            stack.resize(stack.size()-n_pop);
            return result;
        }

        size_t GetStackTop() const { return stack.size(); }
    private:
        void FindClone(CodeTree& /*tree*/, bool /*recurse*/ = true)
        {
            // Disabled: Causes problems in optimization when
            // the same subtree is included in logical and non-logical
            // contexts: optimizations applied to the logical one will
            // mess up the non-logical one.
            return;
            /*
            std::multimap<fphash_t, CodeTree>::const_iterator
                i = clones.lower_bound(tree.GetHash());
            for(; i != clones.end() && i->first == tree.GetHash(); ++i)
            {
                if(i->second.IsIdenticalTo(tree))
                    tree.Become(i->second);
            }
            if(recurse)
                for(size_t a=0; a<tree.GetParamCount(); ++a)
                    FindClone(tree.GetParam(a));
            clones.insert(std::make_pair(tree.GetHash(), tree));
            */
        }
    private:
        std::vector<CodeTree> stack;
        std::multimap<fphash_t, CodeTree> clones;

        bool keep_powi;

    private:
        CodeTreeParserData(const CodeTreeParserData&);
        CodeTreeParserData& operator=(const CodeTreeParserData&);
    };

    struct IfInfo
    {
        CodeTree condition;
        CodeTree thenbranch;
        size_t endif_location;
    };
}

namespace FPoptimizer_CodeTree
{
    void CodeTree::GenerateFrom(
        const std::vector<unsigned>& ByteCode,
        const std::vector<double>& Immed,
        const FunctionParser::Data& fpdata,
        bool keep_powi)
    {
        std::vector<CodeTree> var_trees;
        var_trees.reserve(fpdata.numVariables);
        for(unsigned n=0; n<fpdata.numVariables; ++n)
        {
            var_trees.push_back( CodeTree(n+VarBegin, CodeTree::VarTag()) );
        }
        GenerateFrom(ByteCode,Immed,fpdata,var_trees,keep_powi);
    }

    void CodeTree::GenerateFrom(
        const std::vector<unsigned>& ByteCode,
        const std::vector<double>& Immed,
        const FunctionParser::Data& fpdata,
        const std::vector<CodeTree>& var_trees,
        bool keep_powi)
    {
        CodeTreeParserData sim(keep_powi);
        std::vector<IfInfo> if_stack;

        for(size_t IP=0, DP=0; ; ++IP)
        {
        after_powi:
            while(!if_stack.empty() &&
              (   // Normal If termination rule:
                  if_stack.back().endif_location == IP
                  // This rule matches when cJumps are threaded:
               || (IP < ByteCode.size() && ByteCode[IP] == cJump
                   && if_stack.back().thenbranch.IsDefined())
              ))
            {
                // The "else" of an "if" ends here
                CodeTree elsebranch = sim.PullResult();
                sim.Push(if_stack.back().condition);
                sim.Push(if_stack.back().thenbranch);
                sim.Push(elsebranch);
                sim.Eat(3, cIf);
                if_stack.pop_back();
            }
            if(IP >= ByteCode.size()) break;

            unsigned opcode = ByteCode[IP];
            if((opcode == cSqr || opcode == cDup
             || opcode == cInv || opcode == cNeg
             || opcode == cSqrt || opcode == cRSqrt
             || opcode == cFetch))
            {
                // Parse a powi sequence
                size_t was_ip = IP;
                double exponent = ParsePowiSequence(
                    ByteCode, IP, if_stack.empty() ? ByteCode.size() : if_stack.back().endif_location,
                    sim.GetStackTop()-1);
                if(exponent != 1.0)
                {
                    //std::cout << "Found exponent at " << was_ip << ": " << exponent << "\n";
                    sim.AddConst(exponent);
                    sim.Eat(2, cPow);
                    goto after_powi;
                }
                if(opcode == cDup
                || opcode == cFetch
                || opcode == cNeg)
                {
                    double factor = ParseMuliSequence(
                        ByteCode, IP, if_stack.empty() ? ByteCode.size() : if_stack.back().endif_location,
                        sim.GetStackTop()-1);
                    if(factor != 1.0)
                    {
                        //std::cout << "Found factor at " << was_ip << ": " << factor << "\n";
                        sim.AddConst(factor);
                        sim.Eat(2, cMul);
                        goto after_powi;
                    }
                }
                IP = was_ip;
            }
            if(OPCODE(opcode) >= VarBegin)
            {
                sim.Push(var_trees[opcode-VarBegin]);
            }
            else
            {
                switch( OPCODE(opcode) )
                {
                    // Specials
                    case cIf:
                    case cAbsIf:
                    {
                        if_stack.resize(if_stack.size() + 1);
                        CodeTree res( sim.PullResult() );
                        if_stack.back().condition.swap( res );
                        if_stack.back().endif_location = ByteCode.size();
                        IP += 2; // dp,sp for elsebranch are irrelevant.
                        continue;
                    }
                    case cJump:
                    {
                        CodeTree res( sim.PullResult() );
                        if_stack.back().thenbranch.swap( res );
                        if_stack.back().endif_location = ByteCode[IP+1]+1;
                        IP += 2;
                        continue;
                    }
                    case cImmed:
                        sim.AddConst(Immed[DP++]);
                        break;
                    case cDup:
                        sim.Dup();
                        break;
                    case cNop:
                        break;
                    case cFCall:
                    {
                        unsigned funcno = ByteCode[++IP];
                        assert(funcno < fpdata.FuncPtrs.size());
                        unsigned params = fpdata.FuncPtrs[funcno].params;
                        sim.EatFunc(params, OPCODE(opcode), funcno);
                        break;
                    }
                    case cPCall:
                    {
                        unsigned funcno = ByteCode[++IP];
                        assert(funcno < fpdata.FuncParsers.size());
                        const FunctionParserBase<double>& p =
                            *fpdata.FuncParsers[funcno].parserPtr;
                        unsigned params = fpdata.FuncParsers[funcno].params;

                        /* Inline the procedure call */
                        /* Works because cPCalls can never recurse */
                        std::vector<CodeTree> paramlist = sim.Pop(params);
                        CodeTree pcall_tree;
                        pcall_tree.GenerateFrom(p.data->ByteCode, p.data->Immed, *p.data,
                                                paramlist);
                        sim.Push(pcall_tree);
                        break;
                    }
                    // Unary operators requiring special attention
                    case cInv:  // already handled by powi_opt
                        //sim.Eat(1, cInv);
                        //break;
                        sim.AddConst(1);
                        sim.SwapLastTwoInStack();
                        sim.Eat(2, cDiv);
                        break;
                    case cNeg: // already handled by powi_opt
                        sim.Eat(1, cNeg);
                        break;
                        sim.AddConst(0);
                        sim.SwapLastTwoInStack();
                        sim.Eat(2, cSub);
                        break;
                    case cSqr: // already handled by powi_opt
                        //sim.Eat(1, cSqr);
                        //break;
                        sim.AddConst(2);
                        sim.Eat(2, cPow);
                        break;
                    // Unary functions requiring special attention
                    case cSqrt: // already handled by powi_opt
                        sim.AddConst(0.5);
                        sim.Eat(2, cPow);
                        break;
                    case cRSqrt: // already handled by powi_opt
                        sim.AddConst(-0.5);
                        sim.Eat(2, cPow);
                        break;
                    case cCbrt:
                        sim.AddConst(1.0 / 3.0);
                        sim.Eat(2, cPow);
                        break;
                    case cDeg:
                        sim.AddConst(CONSTANT_DR);
                        sim.Eat(2, cMul);
                        break;
                    case cRad:
                        sim.AddConst(CONSTANT_RD);
                        sim.Eat(2, cMul);
                        break;
                    case cExp:
                        if(keep_powi) goto default_function_handling;
                        sim.AddConst(CONSTANT_E);
                        sim.SwapLastTwoInStack();
                        sim.Eat(2, cPow);
                        break;
                    case cExp2: // from fpoptimizer
                        if(keep_powi) goto default_function_handling;
                        sim.AddConst(2.0);
                        sim.SwapLastTwoInStack();
                        sim.Eat(2, cPow);
                        break;
                    case cCot:
                        sim.Eat(1, cTan);
                        if(keep_powi) { sim.Eat(1, cInv); break; }
                        sim.AddConst(-1);
                        sim.Eat(2, cPow);
                        break;
                    case cCsc:
                        sim.Eat(1, cSin);
                        if(keep_powi) { sim.Eat(1, cInv); break; }
                        sim.AddConst(-1);
                        sim.Eat(2, cPow);
                        break;
                    case cSec:
                        sim.Eat(1, cCos);
                        if(keep_powi) { sim.Eat(1, cInv); break; }
                        sim.AddConst(-1);
                        sim.Eat(2, cPow);
                        break;
                    case cInt: // int(x) = floor(x + 0.5)
                    #ifndef __x86_64
                        if(keep_powi) { sim.Eat(1, cInt); break; }
                    #endif
                        sim.AddConst(0.5);
                        sim.Eat(2, cAdd);
                        sim.Eat(1, cFloor);
                        break;
                    case cLog10:
                        sim.Eat(1, cLog);
                        sim.AddConst(CONSTANT_L10I);
                        sim.Eat(2, cMul);
                        break;
                    case cLog2:
                        sim.Eat(1, cLog);
                        sim.AddConst(CONSTANT_L2I);
                        sim.Eat(2, cMul);
                        break;
                    case cLog2by: // x y     -> log(x)*CONSTANT_L2I*y
                        sim.SwapLastTwoInStack();   // y x
                        sim.Eat(1, cLog);           // y log(x)
                        sim.AddConst(CONSTANT_L2I); // y log(x) CONSTANT_L2I
                        sim.Eat(3, cMul);           // y*log(x)*CONSTANT_L2I
                        break;
                    //case cLog:
                    //    sim.Eat(1, cLog2);
                    //    sim.AddConst(CONSTANT_L2);
                    //    sim.Eat(2, cMul);
                    //    break;
                    // Binary operators requiring special attention
                    case cSub:
                        if(keep_powi) { sim.Eat(2, cSub); break; }
                        sim.AddConst(-1);
                        sim.Eat(2, cMul); // -x is x*-1
                        sim.Eat(2, cAdd); // Minus is negative adding
                        break;
                    case cRSub: // from fpoptimizer
                        sim.SwapLastTwoInStack();
                        if(keep_powi) { sim.Eat(2, cSub); break; }
                        sim.AddConst(-1);
                        sim.Eat(2, cMul); // -x is x*-1
                        sim.Eat(2, cAdd);
                        break;
                    case cDiv:
                        if(keep_powi) { sim.Eat(2, cDiv); break; }
                        sim.AddConst(-1);
                        sim.Eat(2, cPow); // 1/x is x^-1
                        sim.Eat(2, cMul); // Divide is inverse multiply
                        break;
                    case cRDiv: // from fpoptimizer
                        sim.SwapLastTwoInStack();
                        if(keep_powi) { sim.Eat(2, cDiv); break; }
                        sim.AddConst(-1);
                        sim.Eat(2, cPow); // 1/x is x^-1
                        sim.Eat(2, cMul); // Divide is inverse multiply
                        break;
                    // Binary operators not requiring special attention
                    case cAdd: case cMul:
                    case cMod: case cPow:
                    case cEqual: case cLess: case cGreater:
                    case cNEqual: case cLessOrEq: case cGreaterOrEq:
                    case cAnd: case cOr:
                    case cAbsAnd: case cAbsOr:
                        sim.Eat(2, OPCODE(opcode));
                        break;
                    // Unary operators not requiring special attention
                    case cNot:
                    case cNotNot: // from fpoptimizer
                    case cAbsNot:
                    case cAbsNotNot:
                        sim.Eat(1, OPCODE(opcode));
                        break;
                    // Special opcodes generated by fpoptimizer itself
                    case cFetch:
                        sim.Fetch(ByteCode[++IP]);
                        break;
                    case cPopNMov:
                    {
                        unsigned stackOffs_target = ByteCode[++IP];
                        unsigned stackOffs_source = ByteCode[++IP];
                        sim.PopNMov(stackOffs_target, stackOffs_source);
                        break;
                    }
                    // Other functions
#ifndef FP_DISABLE_EVAL
                    case cEval:
                    {
                        size_t paramcount = fpdata.numVariables;
                        sim.Eat(paramcount, OPCODE(opcode));
                        break;
                    }
#endif

                    default:
                    default_function_handling:;
                        unsigned funcno = opcode-cAbs;
                        assert(funcno < FUNC_AMOUNT);
                        const FuncDefinition& func = Functions[funcno];
                        sim.Eat(func.params, OPCODE(opcode));
                        break;
                }
            }
        }
        Become(sim.PullResult());
    #ifdef DEBUG_SUBSTITUTIONS
        std::cout << "Produced tree:\n";
        DumpTreeWithIndent(*this);
    #endif
    }
}

#endif

#line 1 "fpoptimizer/fpoptimizer_constantfolding.c"
// line removed for fpoptimizer.c: #include "fpoptimizer_codetree.h"
// line removed for fpoptimizer.c: #include "fpoptimizer_optimize.h"
// line removed for fpoptimizer.c: #include "fpoptimizer_consts.h"

#include <cmath> /* for CalculateResultBoundaries() */
#include <algorithm>

#include "fpconfig.h"
#include "fparser.h"
#include "fptypes.h"

#ifdef FP_SUPPORT_OPTIMIZER

using namespace FUNCTIONPARSERTYPES;
using namespace FPoptimizer_CodeTree;

#define FP_MUL_COMBINE_EXPONENTS

#ifdef _MSC_VER
#include <float.h>
#define isinf(x) (!_finite(x))
#endif

namespace
{
    bool IsLogicalTrueValue(const MinMaxTree& p, bool abs)
    {
        if(p.has_min && p.min >= 0.5) return true;
        if(!abs && p.has_max && p.max <= -0.5) return true;
        return false;
    }
    bool IsLogicalFalseValue(const MinMaxTree& p, bool abs)
    {
        if(abs)
            return p.has_max && p.max < 0.5;
        else
            return p.has_min && p.has_max
               && p.min > -0.5 && p.max < 0.5;
    }
    int GetLogicalValue(const MinMaxTree& p, bool abs)
    {
        if(IsLogicalTrueValue(p, abs)) return 1;
        if(IsLogicalFalseValue(p, abs)) return 0;
        return -1;
    }

    struct ComparisonSet /* For optimizing And, Or */
    {
        static const int Lt_Mask = 0x1; // 1=less
        static const int Eq_Mask = 0x2; // 2=equal
        static const int Le_Mask = 0x3; // 1+2 = Less or Equal
        static const int Gt_Mask = 0x4; // 4=greater
        static const int Ne_Mask = 0x5; // 4+1 = Greater or Less, i.e. Not equal
        static const int Ge_Mask = 0x6; // 4+2 = Greater or Equal
        static int Swap_Mask(int m) { return (m&Eq_Mask)
                                  | ((m&Lt_Mask) ? Gt_Mask : 0)
                                  | ((m&Gt_Mask) ? Lt_Mask : 0); }
        struct Comparison
        {
            CodeTree a;
            CodeTree b;
            int relationship;
        };
        std::vector<Comparison> relationships;
        struct Item
        {
            CodeTree value;
            bool negated;
        };
        std::vector<Item> plain_set;
        int const_offset;

        enum RelationshipResult
        {
            Ok,
            BecomeZero,
            BecomeOne,
            Suboptimal
        };
        enum ConditionType
        {
            cond_or,
            cond_and,
            cond_mul,
            cond_add
        };

        ComparisonSet():
            relationships(),
            plain_set(),
            const_offset(0)
        {
        }

        RelationshipResult AddItem(const CodeTree& a, bool negated, ConditionType type)
        {
            for(size_t c=0; c<plain_set.size(); ++c)
                if(plain_set[c].value.IsIdenticalTo(a))
                {
                    if(negated != plain_set[c].negated)
                    {
                        switch(type)
                        {
                            case cond_or:
                                return BecomeOne;
                            case cond_add:
                                plain_set.erase(plain_set.begin() + c);
                                const_offset += 1;
                                return Suboptimal;
                            case cond_and:
                            case cond_mul:
                                return BecomeZero;
                        }
                    }
                    return Suboptimal;
                }
            Item pole;
            pole.value   = a;
            pole.negated = negated;
            plain_set.push_back(pole);
            return Ok;
        }

        RelationshipResult AddRelationship(CodeTree a, CodeTree b, int reltype, ConditionType type)
        {
            switch(type)
            {
                case cond_or:
                    if(reltype == 7) return BecomeOne;
                    break;
                case cond_add:
                    if(reltype == 7) { const_offset += 1; return Suboptimal; }
                    break;
                case cond_and:
                case cond_mul:
                    if(reltype == 0) return BecomeZero;
                    break;
            }

            if(!(a.GetHash() < b.GetHash()))
            {
                a.swap(b);
                reltype = Swap_Mask(reltype);
            }

            for(size_t c=0; c<relationships.size(); ++c)
            {
                if(relationships[c].a.IsIdenticalTo(a)
                && relationships[c].b.IsIdenticalTo(b))
                {
                    switch(type)
                    {
                        case cond_or:
                        {
                            int newrel = relationships[c].relationship | reltype;
                            if(newrel == 7) return BecomeOne;
                            relationships[c].relationship = newrel;
                            break;
                        }
                        case cond_and:
                        case cond_mul:
                        {
                            int newrel = relationships[c].relationship & reltype;
                            if(newrel == 0) return BecomeZero;
                            relationships[c].relationship = newrel;
                            break;
                        }
                        case cond_add:
                        {
                            int newrel_or  = relationships[c].relationship | reltype;
                            int newrel_and = relationships[c].relationship & reltype;
                            if(newrel_or  == 5 // < + >
                            && newrel_and == 0)
                            {
                                // (x<y) + (x>y) = x!=y
                                relationships[c].relationship = Ne_Mask;
                                return Suboptimal;
                            }
                            if(newrel_or  == 7
                            && newrel_and == 0)
                            {
                                // (x<y) + (x>=y) = 1
                                // (x<=y) + (x>y) = 1
                                // (x=y) + (x!=y) = 1
                                const_offset += 1;
                                relationships.erase(relationships.begin()+c);
                                return Suboptimal;
                            }
                            if(newrel_or  == 7
                            && newrel_and == Eq_Mask)
                            {
                                // (x<=y) + (x>=y) = 1 + (x=y)
                                relationships[c].relationship = Eq_Mask;
                                const_offset += 1;
                                return Suboptimal;
                            }
                            continue;
                        }
                    }
                    return Suboptimal;
                }
            }
            Comparison comp;
            comp.a = a;
            comp.b = b;
            comp.relationship = reltype;
            relationships.push_back(comp);
            return Ok;
        }
    };

    struct CollectionSet /* For optimizing Add,  Mul */
    {
        struct Collection
        {
            CodeTree value;
            CodeTree factor;
            bool factor_needs_rehashing;

            Collection() : value(),factor(), factor_needs_rehashing(false) { }
            Collection(const CodeTree& v, const CodeTree& f)
                : value(v), factor(f), factor_needs_rehashing(false) { }
        };
        std::multimap<fphash_t, Collection> collections;

        enum CollectionResult
        {
            Ok,
            Suboptimal
        };

        typedef std::multimap<fphash_t, Collection>::iterator PositionType;

        PositionType FindIdenticalValueTo(const CodeTree& value)
        {
            fphash_t hash = value.GetHash();
            for(PositionType
                i = collections.lower_bound(hash);
                i != collections.end() && i->first == hash;
                ++i)
            {
                if(value.IsIdenticalTo(i->second.value))
                    return i;
            }
            return collections.end();
        }
        bool Found(const PositionType& b) { return b != collections.end(); }

        CollectionResult AddCollectionTo(const CodeTree& factor,
                                         const PositionType& into_which)
        {
            Collection& c = into_which->second;
            if(c.factor_needs_rehashing)
                c.factor.AddParam(factor);
            else
            {
                CodeTree add;
                add.SetOpcode(cAdd);
                add.AddParamMove(c.factor);
                add.AddParam(factor);
                c.factor.swap(add);
                c.factor_needs_rehashing = true;
            }
            return Suboptimal;
        }

        CollectionResult AddCollection(const CodeTree& value, const CodeTree& factor)
        {
            const fphash_t hash = value.GetHash();
            PositionType i = collections.lower_bound(hash);
            for(; i != collections.end() && i->first == hash; ++i)
            {
                if(i->second.value.IsIdenticalTo(value))
                    return AddCollectionTo(factor, i);
            }
            collections.insert(
                i,
                std::make_pair( hash, Collection(value, factor) ) );
            return Ok;
        }

        CollectionResult AddCollection(const CodeTree& a)
        {
            return AddCollection(a, CodeTree(1.0) );
        }
    };

    struct Select2ndRev
    {
        template<typename T>
        inline bool operator() (const T& a, const T& b) const
        {
            return a.second > b.second;
        }
    };
    struct Select1st
    {
        template<typename T1, typename T2>
        inline bool operator() (const std::pair<T1,T2>& a,
                                const std::pair<T1,T2>& b) const
        {
            return a.first < b.first;
        }

        template<typename T1, typename T2>
        inline bool operator() (const std::pair<T1,T2>& a, T1 b) const
        {
            return a.first < b;
        }

        template<typename T1, typename T2>
        inline bool operator() (T1 a, const std::pair<T1,T2>& b) const
        {
            return a < b.first;
        }
    };

    bool IsEvenIntegerConst(double v)
    {
        return IsIntegerConst(v) && ((long)v % 2) == 0;
    }

    struct ConstantExponentCollection
    {
        typedef std::pair<double, std::vector<CodeTree> > ExponentInfo;
        std::vector<ExponentInfo> data;

        void MoveToSet_Unique(double exponent, std::vector<CodeTree>& source_set)
        {
            data.push_back( std::pair<double, std::vector<CodeTree> >
                            (exponent, std::vector<CodeTree>() ) );
            data.back().second.swap(source_set);
        }
        void MoveToSet_NonUnique(double exponent, std::vector<CodeTree>& source_set)
        {
            std::vector<ExponentInfo>::iterator i
                = std::lower_bound(data.begin(), data.end(), exponent, Select1st());
            if(i != data.end() && i->first == exponent)
            {
                i->second.insert(i->second.end(), source_set.begin(), source_set.end());
            }
            else
            {
                //MoveToSet_Unique(exponent, source_set);
                data.insert(i,  std::pair<double, std::vector<CodeTree> >
                                (exponent, source_set) );
            }
        }

        bool Optimize()
        {
            /* TODO: Group them such that:
             *
             *      x^3 *         z^2 becomes (x*z)^2 * x^1
             *      x^3 * y^2.5 * z^2 becomes (x*z*y)^2 * y^0.5 * x^1
             *                    rather than (x*y*z)^2 * (x*y)^0.5 * x^0.5
             *
             *      x^4.5 * z^2.5     becomes (z * x)^2.5 * x^2
             *                        becomes (x*z*x)^2 * (z*x)^0.5
             *                        becomes (z*x*x*z*x)^0.5 * (z*x*x)^1.5 -- buzz, bad.
             *
             */
            bool changed = false;
            std::sort( data.begin(), data.end(), Select1st() );
        redo:
            /* Supposed algorithm:
             * For the smallest pair of data[] where the difference
             * between the two is a "neat value" (x*16 is positive integer),
             * do the combining as indicated above.
             */
            /*
             * NOTE: Hanged in Testbed test P44, looping the following
             *       (Var0 ^ 0.75) * ((1.5 * Var0) ^ 1.0)
             *     = (Var0 ^ 1.75) *  (1.5         ^ 1.0)
             *       Fixed by limiting to cases where (exp_a != 1.0).
             *
             * NOTE: Converting (x*z)^0.5 * x^16.5
             *              into x^17 * z^0.5
             *       is handled by code within CollectMulGroup().
             *       However, bacause it is prone for infinite looping,
             *       the use of "IsIdenticalTo(before)" is added at the
             *       end of ConstantFolding_MulGrouping().
             *
             *       This algorithm could make it into (x*z*x)^0.5 * x^16,
             *       but this is wrong, for it falsely includes x^evenint.. twice.
             */
            for(size_t a=0; a<data.size(); ++a)
            {
                double exp_a = data[a].first;
                if(FloatEqual(exp_a, 1.0)) continue;
                for(size_t b=a+1; b<data.size(); ++b)
                {
                    double exp_b = data[b].first;
                    double exp_diff = exp_b - exp_a;
                    if(exp_diff >= fabs(exp_a)) break;
                    double exp_diff_still_probable_integer = exp_diff * 16.0;
                    if(IsIntegerConst(exp_diff_still_probable_integer)
                    && !(IsIntegerConst(exp_b) && !IsIntegerConst(exp_diff))
                      )
                    {
                        /* When input is x^3 * z^2,
                         * exp_a = 2
                         * a_set = z
                         * exp_b = 3
                         * b_set = x
                         * exp_diff = 3-2 = 1
                         */
                        std::vector<CodeTree>& a_set = data[a].second;
                        std::vector<CodeTree>& b_set = data[b].second;
          #ifdef DEBUG_SUBSTITUTIONS
                        std::cout << "Before ConstantExponentCollection iteration:\n";
                        Dump(std::cout);
          #endif
                        if(IsIntegerConst(exp_b)
                        && IsEvenIntegerConst(exp_b)
                        //&& !IsEvenIntegerConst(exp_diff)
                        && !IsEvenIntegerConst(exp_diff+exp_a))
                        {
                            CodeTree tmp2;
                            tmp2.SetOpcode(cMul);
                            tmp2.SetParamsMove(b_set);
                            tmp2.Rehash();
                            CodeTree tmp;
                            tmp.SetOpcode(cAbs);
                            tmp.AddParamMove(tmp2);
                            tmp.Rehash();
                            b_set.resize(1);
                            b_set[0].swap(tmp);
                        }

                        a_set.insert(a_set.end(), b_set.begin(), b_set.end());

                        std::vector<CodeTree> b_copy = b_set;
                        data.erase(data.begin() + b);
                        MoveToSet_NonUnique(exp_diff, b_copy);
                        changed = true;

          #ifdef DEBUG_SUBSTITUTIONS
                        std::cout << "After ConstantExponentCollection iteration:\n";
                        Dump(std::cout);
          #endif
                        goto redo;
                    }
                }
            }
            return changed;
        }

    #ifdef DEBUG_SUBSTITUTIONS
        void Dump(std::ostream& out)
        {
            for(size_t a=0; a<data.size(); ++a)
            {
                out.precision(12);
                out << data[a].first << ": ";
                for(size_t b=0; b<data[a].second.size(); ++b)
                {
                    if(b > 0) out << '*';
                    DumpTree(data[a].second[b], out);
                }
                out << std::endl;
            }
        }
    #endif

    };

    struct RangeComparisonData
    {
        enum Decision
        {
            MakeFalse=0,
            MakeTrue=1,
            MakeNEqual=2,
            MakeEqual=3,
            MakeNotNotP0=4,
            MakeNotNotP1=5,
            MakeNotP0=6,
            MakeNotP1=7,
            Unchanged=8
        };
        enum WhatDoWhenCase
        {
            Never =0,
            Eq0   =1, // val==0
            Eq1   =2, // val==1
            Gt0Le1=3, // val>0 && val<=1
            Ge0Lt1=4  // val>=0 && val<1
        };
        Decision if_identical; // What to do when operands are identical
        Decision if_always[4]; // What to do if Always <, <=, >, >=
        struct { Decision what : 4; WhatDoWhenCase when : 4; }
            p0_logical_a, p1_logical_a,
            p0_logical_b, p1_logical_b;

        Decision Analyze(const CodeTree& a, const CodeTree& b) const
        {
            if(a.IsIdenticalTo(b))
                return if_identical;

            MinMaxTree p0 = a.CalculateResultBoundaries();
            MinMaxTree p1 = b.CalculateResultBoundaries();
            if(p0.has_max && p1.has_min)
            {
                if(p0.max <  p1.min && if_always[0] != Unchanged)
                    return if_always[0]; // p0 < p1
                if(p0.max <= p1.min && if_always[1] != Unchanged)
                    return if_always[1]; // p0 <= p1
            }
            if(p0.has_min && p1.has_max)
            {
                if(p0.min >  p1.max && if_always[2] != Unchanged)
                    return if_always[2]; // p0 > p1
                if(p0.min >= p1.max && if_always[3] != Unchanged)
                    return if_always[3]; // p0 >= p1
            }

            if(a.IsLogicalValue())
            {
                if(p0_logical_a.what != Unchanged)
                    if(TestCase(p0_logical_a.when, p1)) return p0_logical_a.what;
                if(p0_logical_b.what != Unchanged)
                    if(TestCase(p0_logical_b.when, p1)) return p0_logical_b.what;
            }
            if(b.IsLogicalValue())
            {
                if(p1_logical_a.what != Unchanged)
                    if(TestCase(p1_logical_a.when, p0)) return p1_logical_a.what;
                if(p1_logical_b.what != Unchanged)
                    if(TestCase(p1_logical_b.when, p0)) return p1_logical_b.what;
            }
            return Unchanged;
        }
        static bool TestCase(WhatDoWhenCase when, const MinMaxTree& p)
        {
            if(!p.has_min || !p.has_max) return false;
            switch(when)
            {
                case Eq0: return p.min==0.0 && p.max==p.min;
                case Eq1: return p.min==1.0 && p.max==p.max;
                case Gt0Le1: return p.min>0 && p.max<=1;
                case Ge0Lt1: return p.min>=0 && p.max<1;
                default:;
            }
            return false;
        }
    };

}

namespace FPoptimizer_CodeTree
{
    template<typename CondType> /* ComparisonSet::ConditionType */
    bool CodeTree::ConstantFolding_LogicCommon(CondType cond_type, bool is_logical)
    {
        bool should_regenerate = false;
        ComparisonSet comp;
        for(size_t a=0; a<GetParamCount(); ++a)
        {
            ComparisonSet::RelationshipResult change = ComparisonSet::Ok;
            const CodeTree& atree = GetParam(a);
            switch(atree.GetOpcode())
            {
                case cEqual:
                    change = comp.AddRelationship(atree.GetParam(0), atree.GetParam(1), ComparisonSet::Eq_Mask, cond_type);
                    break;
                case cNEqual:
                    change = comp.AddRelationship(atree.GetParam(0), atree.GetParam(1), ComparisonSet::Ne_Mask, cond_type);
                    break;
                case cLess:
                    change = comp.AddRelationship(atree.GetParam(0), atree.GetParam(1), ComparisonSet::Lt_Mask, cond_type);
                    break;
                case cLessOrEq:
                    change = comp.AddRelationship(atree.GetParam(0), atree.GetParam(1), ComparisonSet::Le_Mask, cond_type);
                    break;
                case cGreater:
                    change = comp.AddRelationship(atree.GetParam(0), atree.GetParam(1), ComparisonSet::Gt_Mask, cond_type);
                    break;
                case cGreaterOrEq:
                    change = comp.AddRelationship(atree.GetParam(0), atree.GetParam(1), ComparisonSet::Ge_Mask, cond_type);
                    break;
                case cNot:
                    change = comp.AddItem(atree.GetParam(0), true, cond_type);
                    break;
                case cNotNot:
                    change = comp.AddItem(atree.GetParam(0), false, cond_type);
                    break;
                default:
                    if(is_logical || atree.IsLogicalValue())
                        change = comp.AddItem(atree, false, cond_type);
            }
            switch(change)
            {
            ReplaceTreeWithZero:
                    data = new CodeTreeData(0.0);
                    return true;
            ReplaceTreeWithOne:
                    data = new CodeTreeData(1.0);
                    return true;
                case ComparisonSet::Ok: // ok
                    break;
                case ComparisonSet::BecomeZero: // whole set was invalidated
                    goto ReplaceTreeWithZero;
                case ComparisonSet::BecomeOne: // whole set was validated
                    goto ReplaceTreeWithOne;
                case ComparisonSet::Suboptimal: // something was changed
                    should_regenerate = true;
                    break;
            }
        }
        if(should_regenerate)
        {
          #ifdef DEBUG_SUBSTITUTIONS
            std::cout << "Before ConstantFolding_LogicCommon: "; DumpTree(*this);
            std::cout << "\n";
          #endif

            if(is_logical)
            {
                DelParams(); // delete all params
            }
            else
            {
                // Delete only logical params
                for(size_t a=GetParamCount(); a-- > 0; )
                {
                    const CodeTree& atree = GetParam(a);
                    if(atree.IsLogicalValue())
                        DelParam(a);
                }
            }

            for(size_t a=0; a<comp.plain_set.size(); ++a)
            {
                if(comp.plain_set[a].negated)
                {
                    CodeTree r;
                    r.SetOpcode(cNot);
                    r.AddParamMove(comp.plain_set[a].value);
                    r.Rehash();
                    AddParamMove(r);
                }
                else if(!is_logical)
                {
                    CodeTree r;
                    r.SetOpcode(cNotNot);
                    r.AddParamMove(comp.plain_set[a].value);
                    r.Rehash();
                    AddParamMove(r);
                }
                else
                    AddParamMove(comp.plain_set[a].value);
            }
            for(size_t a=0; a<comp.relationships.size(); ++a)
            {
                CodeTree r;
                r.SetOpcode(cNop); // dummy
                switch(comp.relationships[a].relationship)
                {
                    case ComparisonSet::Lt_Mask: r.SetOpcode( cLess ); break;
                    case ComparisonSet::Eq_Mask: r.SetOpcode( cEqual ); break;
                    case ComparisonSet::Gt_Mask: r.SetOpcode( cGreater ); break;
                    case ComparisonSet::Le_Mask: r.SetOpcode( cLessOrEq ); break;
                    case ComparisonSet::Ne_Mask: r.SetOpcode( cNEqual ); break;
                    case ComparisonSet::Ge_Mask: r.SetOpcode( cGreaterOrEq ); break;
                }
                r.AddParamMove(comp.relationships[a].a);
                r.AddParamMove(comp.relationships[a].b);
                r.Rehash();
                AddParamMove(r);
            }
            if(comp.const_offset != 0)
                AddParam( CodeTree( double(comp.const_offset) ) );
          #ifdef DEBUG_SUBSTITUTIONS
            std::cout << "After ConstantFolding_LogicCommon: "; DumpTree(*this);
            std::cout << "\n";
          #endif
            return true;
        }
        /*
        Note: One thing this does not yet do, is to detect chains
              such as x=y & y=z & x=z, which could be optimized
              to x=y & x=z.
        */
        return false;
    }

    bool CodeTree::ConstantFolding_AndLogic()
    {
        return ConstantFolding_LogicCommon( ComparisonSet::cond_and, true );
    }
    bool CodeTree::ConstantFolding_OrLogic()
    {
        return ConstantFolding_LogicCommon( ComparisonSet::cond_or, true );
    }
    bool CodeTree::ConstantFolding_AddLogicItems()
    {
        return ConstantFolding_LogicCommon( ComparisonSet::cond_add, false );
    }
    bool CodeTree::ConstantFolding_MulLogicItems()
    {
        return ConstantFolding_LogicCommon( ComparisonSet::cond_mul, false );
    }

    static CodeTree CollectMulGroup_Item(
        CodeTree& value,
        bool& has_highlevel_opcodes)
    {
        switch(value.GetOpcode())
        {
            case cPow:
            {
                CodeTree exponent = value.GetParam(1);
                value.Become( value.GetParam(0) );
                return exponent;
            }
            /* - disabled to avoid clashes with powi
            case cCbrt:
                value.Become( value.GetParam(0) );
                has_highlevel_opcodes = true;
                return CodeTree(1.0 / 3.0);
            case cSqrt:
                value.Become( value.GetParam(0) );
                has_highlevel_opcodes = true;
                return CodeTree(0.5);
            */
            case cRSqrt:
                value.Become( value.GetParam(0) );
                has_highlevel_opcodes = true;
                return CodeTree(-0.5);
            case cInv:
                value.Become( value.GetParam(0) );
                has_highlevel_opcodes = true;
                return CodeTree(-1.0);
            default: break;
        }
        return CodeTree(1.0);
    }

    static void CollectMulGroup(
        CollectionSet& mul, const CodeTree& tree, const CodeTree& factor,
        bool& should_regenerate,
        bool& has_highlevel_opcodes
    )
    {
        for(size_t a=0; a<tree.GetParamCount(); ++a)
        {
            CodeTree value(tree.GetParam(a));

            CodeTree exponent ( CollectMulGroup_Item(value, has_highlevel_opcodes) );

            if(!factor.IsImmed() || factor.GetImmed() != 1.0)
            {
                CodeTree new_exp;
                new_exp.SetOpcode(cMul);
                new_exp.AddParam( exponent );
                new_exp.AddParam( factor );
                new_exp.Rehash();
                exponent.swap( new_exp );
            }
        #if 0 /* FIXME: This does not work */
            if(value.GetOpcode() == cMul)
            {
                if(1)
                {
                    // Avoid erroneously converting
                    //          (x*z)^0.5 * z^2
                    // into     x^0.5 * z^2.5
                    // It should be x^0.5 * abs(z)^2.5, but this is not a good conversion.
                    bool exponent_is_even = exponent.IsImmed() && IsEvenIntegerConst(exponent.GetImmed());

                    for(size_t b=0; b<value.GetParamCount(); ++b)
                    {
                        bool tmp=false;
                        CodeTree val(value.GetParam(b));
                        CodeTree exp(CollectMulGroup_Item(val, tmp));
                        if(exponent_is_even
                        || (exp.IsImmed() && IsEvenIntegerConst(exp.GetImmed())))
                        {
                            CodeTree new_exp;
                            new_exp.SetOpcode(cMul);
                            new_exp.AddParam(exponent);
                            new_exp.AddParamMove(exp);
                            new_exp.ConstantFolding();
                            if(!new_exp.IsImmed() || !IsEvenIntegerConst(new_exp.GetImmed()))
                            {
                                goto cannot_adopt_mul;
                            }
                        }
                    }
                }
                CollectMulGroup(mul, value, exponent,
                                should_regenerate,
                                has_highlevel_opcodes);
            }
            else cannot_adopt_mul:
        #endif
            {
                if(mul.AddCollection(value, exponent) == CollectionSet::Suboptimal)
                    should_regenerate = true;
            }
        }
    }

    bool CodeTree::ConstantFolding_MulGrouping()
    {
        bool has_highlevel_opcodes = false;
        bool should_regenerate = false;
        CollectionSet mul;

        CollectMulGroup(mul, *this, CodeTree(1.0),
                        should_regenerate,
                        has_highlevel_opcodes);

        typedef std::pair<CodeTree/*exponent*/,
                          std::vector<CodeTree>/*base value (mul group)*/
                         > exponent_list;
        typedef std::multimap<fphash_t,/*exponent hash*/
                              exponent_list> exponent_map;
        exponent_map by_exponent;

        for(CollectionSet::PositionType
            j = mul.collections.begin();
            j != mul.collections.end();
            ++j)
        {
            CodeTree& value = j->second.value;
            CodeTree& exponent = j->second.factor;
            if(j->second.factor_needs_rehashing) exponent.Rehash();
            const fphash_t exponent_hash = exponent.GetHash();

            exponent_map::iterator i = by_exponent.lower_bound(exponent_hash);
            for(; i != by_exponent.end() && i->first == exponent_hash; ++i)
                if(i->second.first.IsIdenticalTo(exponent))
                {
                    if(!exponent.IsImmed() || !FloatEqual(exponent.GetImmed(), 1.0))
                        should_regenerate = true;
                    i->second.second.push_back(value);
                    goto skip_b;
                }
            by_exponent.insert(i, std::make_pair(exponent_hash,
                std::make_pair(exponent,
                               std::vector<CodeTree> (size_t(1), value)
                              )));
        skip_b:;
        }

    #ifdef FP_MUL_COMBINE_EXPONENTS
        ConstantExponentCollection by_float_exponent;
        for(exponent_map::iterator
            j,i = by_exponent.begin();
            i != by_exponent.end();
            i=j)
        {
            j=i; ++j;
            exponent_list& list = i->second;
            if(list.first.IsImmed())
            {
                double exponent = list.first.GetImmed();
                if(!(exponent == 0.0))
                    by_float_exponent.MoveToSet_Unique(exponent, list.second);
                by_exponent.erase(i);
            }
        }
        if(by_float_exponent.Optimize())
            should_regenerate = true;
    #endif

        if(should_regenerate)
        {
            CodeTree before = *this;
            before.CopyOnWrite();

          #ifdef DEBUG_SUBSTITUTIONS
            std::cout << "Before ConstantFolding_MulGrouping: "; DumpTree(before);
            std::cout << "\n";
          #endif
            DelParams();

            /* Group by exponents */
            /* First handle non-constant exponents */
            for(exponent_map::iterator
                i = by_exponent.begin();
                i != by_exponent.end();
                ++i)
            {
                exponent_list& list = i->second;
        #ifndef FP_MUL_COMBINE_EXPONENTS
                if(list.first.IsImmed())
                {
                    double exponent = list.first.GetImmed();
                    if(exponent == 0.0) continue;
                    if(FloatEqual(exponent, 1.0))
                    {
                        AddParamsMove(list.second);
                        continue;
                    }
                }
        #endif
                CodeTree mul;
                mul.SetOpcode(cMul);
                mul.SetParamsMove( list.second);
                mul.Rehash();

                if(has_highlevel_opcodes && list.first.IsImmed())
                {
                    if(list.first.GetImmed() == 1.0 / 3.0)
                    {
                        CodeTree cbrt;
                        cbrt.SetOpcode(cCbrt);
                        cbrt.AddParamMove(mul);
                        cbrt.Rehash();
                        AddParamMove(cbrt);
                        continue;
                    }
                    if(list.first.GetImmed() == 0.5)
                    {
                        CodeTree sqrt;
                        sqrt.SetOpcode(cSqrt);
                        sqrt.AddParamMove(mul);
                        sqrt.Rehash();
                        AddParamMove(sqrt);
                        continue;
                    }
                    if(list.first.GetImmed() == -0.5)
                    {
                        CodeTree rsqrt;
                        rsqrt.SetOpcode(cRSqrt);
                        rsqrt.AddParamMove(mul);
                        rsqrt.Rehash();
                        AddParamMove(rsqrt);
                        continue;
                    }
                    if(list.first.GetImmed() == -1.0)
                    {
                        CodeTree inv;
                        inv.SetOpcode(cInv);
                        inv.AddParamMove(mul);
                        inv.Rehash();
                        AddParamMove(inv);
                        continue;
                    }
                }
                CodeTree pow;
                pow.SetOpcode(cPow);
                pow.AddParamMove(mul);
                pow.AddParamMove( list.first );
                pow.Rehash();
                AddParamMove(pow);
            }
        #ifdef FP_MUL_COMBINE_EXPONENTS
            by_exponent.clear();
            /* Then handle constant exponents */
            for(size_t a=0; a<by_float_exponent.data.size(); ++a)
            {
                double exponent = by_float_exponent.data[a].first;
                if(FloatEqual(exponent, 1.0))
                {
                    AddParamsMove(by_float_exponent.data[a].second);
                    continue;
                }
                CodeTree mul;
                mul.SetOpcode(cMul);
                mul.SetParamsMove( by_float_exponent.data[a].second );
                mul.Rehash();
                CodeTree pow;
                pow.SetOpcode(cPow);
                pow.AddParamMove(mul);
                pow.AddParam( CodeTree( exponent ) );
                pow.Rehash();
                AddParamMove(pow);
            }
        #endif
          #ifdef DEBUG_SUBSTITUTIONS
            std::cout << "After ConstantFolding_MulGrouping: "; DumpTree(*this);
            std::cout << "\n";
          #endif
            // return true;
            return !IsIdenticalTo(before); // avoids infinite looping
        }
        return false;
    }

    bool CodeTree::ConstantFolding_AddGrouping()
    {
        bool should_regenerate = false;
        CollectionSet add;
        for(size_t a=0; a<GetParamCount(); ++a)
        {
            if(GetParam(a).GetOpcode() == cMul) continue;
            if(add.AddCollection(GetParam(a)) == CollectionSet::Suboptimal)
                should_regenerate = true;
            // This catches x + x and x - x
        }
        std::vector<bool> remaining ( GetParamCount() );
        size_t has_mulgroups_remaining = 0;
        for(size_t a=0; a<GetParamCount(); ++a)
        {
            const CodeTree& mulgroup = GetParam(a);
            if(mulgroup.GetOpcode() == cMul)
            {
                // This catches x + y*x*z, producing x*(1 + y*z)
                //
                // However we avoid changing 7 + 7*x into 7*(x+1),
                // because it may lead us into producing code such
                // as 20*x + 50*(x+1) + 10, which would be much
                // better expressed as 70*x + 60, and converting
                // back to that format would be needlessly hairy.
                for(size_t b=0; b<mulgroup.GetParamCount(); ++b)
                {
                    if(mulgroup.GetParam(b).IsImmed()) continue;
                    CollectionSet::PositionType c
                        = add.FindIdenticalValueTo(mulgroup.GetParam(b));
                    if(add.Found(c))
                    {
                        CodeTree tmp(mulgroup, CodeTree::CloneTag());
                        tmp.DelParam(b);
                        tmp.Rehash();
                        add.AddCollectionTo(tmp, c);
                        should_regenerate = true;
                        goto done_a;
                    }
                }
                remaining[a]  = true;
                has_mulgroups_remaining += 1;
            done_a:;
            }
        }

        if(has_mulgroups_remaining > 0)
        {
            if(has_mulgroups_remaining > 1) // is it possible to find a duplicate?
            {
                std::vector< std::pair<CodeTree, size_t> > occurance_counts;
                std::multimap<fphash_t, size_t> occurance_pos;
                bool found_dup = false;
                for(size_t a=0; a<GetParamCount(); ++a)
                    if(remaining[a])
                    {
                        // This catches x*a + x*b, producing x*(a+b)
                        for(size_t b=0; b<GetParam(a).GetParamCount(); ++b)
                        {
                            const CodeTree& p = GetParam(a).GetParam(b);
                            const fphash_t   p_hash = p.GetHash();
                            for(std::multimap<fphash_t, size_t>::const_iterator
                                i = occurance_pos.lower_bound(p_hash);
                                i != occurance_pos.end() && i->first == p_hash;
                                ++i)
                            {
                                if(occurance_counts[i->second].first.IsIdenticalTo(p))
                                {
                                    occurance_counts[i->second].second += 1;
                                    found_dup = true;
                                    goto found_mulgroup_item_dup;
                                }
                            }
                            occurance_counts.push_back(std::make_pair(p, size_t(1)));
                            occurance_pos.insert(std::make_pair(p_hash, occurance_counts.size()-1));
                        found_mulgroup_item_dup:;
                        }
                    }
                if(found_dup)
                {
                    // Find the "x" to group by
                    CodeTree group_by; { size_t max = 0;
                    for(size_t p=0; p<occurance_counts.size(); ++p)
                        if(occurance_counts[p].second <= 1)
                            occurance_counts[p].second = 0;
                        else
                        {
                            occurance_counts[p].second *= occurance_counts[p].first.GetDepth();
                            if(occurance_counts[p].second > max)
                                { group_by = occurance_counts[p].first; max = occurance_counts[p].second; }
                        } }
                    // Collect the items for adding in the group (a+b)
                    CodeTree group_add;
                    group_add.SetOpcode(cAdd);

        #ifdef DEBUG_SUBSTITUTIONS
                    std::cout << "Duplicate across some trees: ";
                    DumpTree(group_by);
                    std::cout << " in ";
                    DumpTree(*this);
                    std::cout << "\n";
        #endif
                    for(size_t a=0; a<GetParamCount(); ++a)
                        if(remaining[a])
                            for(size_t b=0; b<GetParam(a).GetParamCount(); ++b)
                                if(group_by.IsIdenticalTo(GetParam(a).GetParam(b)))
                                {
                                    CodeTree tmp(GetParam(a), CodeTree::CloneTag());
                                    tmp.DelParam(b);
                                    tmp.Rehash();
                                    group_add.AddParamMove(tmp);
                                    remaining[a] = false;
                                    break;
                                }
                    group_add.Rehash();
                    CodeTree group;
                    group.SetOpcode(cMul);
                    group.AddParamMove(group_by);
                    group.AddParamMove(group_add);
                    group.Rehash();
                    add.AddCollection(group);
                    should_regenerate = true;
                }
            }

            // all remaining mul-groups.
            for(size_t a=0; a<GetParamCount(); ++a)
                if(remaining[a])
                {
                    if(add.AddCollection(GetParam(a)) == CollectionSet::Suboptimal)
                        should_regenerate = true;
                }
        }

        if(should_regenerate)
        {
          #ifdef DEBUG_SUBSTITUTIONS
            std::cout << "Before ConstantFolding_AddGrouping: "; DumpTree(*this);
            std::cout << "\n";
          #endif
            DelParams();

            for(CollectionSet::PositionType
                j = add.collections.begin();
                j != add.collections.end();
                ++j)
            {
                CodeTree& value = j->second.value;
                CodeTree& coeff = j->second.factor;
                if(j->second.factor_needs_rehashing) coeff.Rehash();

                if(coeff.IsImmed())
                {
                    if(coeff.GetImmed() == 0.0)
                        continue;
                    if(FloatEqual(coeff.GetImmed(), 1.0))
                    {
                        AddParamMove(value);
                        continue;
                    }
                }
                CodeTree mul;
                mul.SetOpcode(cMul);
                mul.AddParamMove(value);
                mul.AddParamMove(coeff);
                mul.Rehash();
                AddParamMove(mul);
            }
          #ifdef DEBUG_SUBSTITUTIONS
            std::cout << "After ConstantFolding_AddGrouping: "; DumpTree(*this);
            std::cout << "\n";
          #endif
            return true;
        }
        return false;
    }

    bool CodeTree::ConstantFolding_IfOperations()
    {
        // If the If() condition begins with a cNot,
        // remove the cNot and swap the branches.
        for(;;)
        {
            if(GetParam(0).GetOpcode() == cNot)
            {
                SetOpcode(cIf);
                GetParam(0).Become( GetParam(0).GetParam(0) );
                GetParam(1).swap(GetParam(2));
            }
            else if(GetParam(0).GetOpcode() == cAbsNot)
            {
                SetOpcode(cAbsIf);
                GetParam(0).Become( GetParam(0).GetParam(0) );
                GetParam(1).swap(GetParam(2));
            }
            else break;
        }
        if(GetParam(0).GetOpcode() == cIf
        || GetParam(0).GetOpcode() == cAbsIf)
        {
            //     if(if(x, a,b), c,d)
            //  -> if(x, if(a, c,d), if(b, c,d))
            // when either a or b is constantly true/false
            CodeTree cond = GetParam(0);
            CodeTree truth_a;
            truth_a.SetOpcode(cond.GetOpcode() == cIf ? cNotNot : cAbsNotNot);
            truth_a.AddParam(cond.GetParam(1));
            truth_a.ConstantFolding();
            CodeTree truth_b;
            truth_b.SetOpcode(cond.GetOpcode() == cIf ? cNotNot : cAbsNotNot);
            truth_b.AddParam(cond.GetParam(2));
            truth_b.ConstantFolding();
            if(truth_a.IsImmed() || truth_b.IsImmed())
            {
                CodeTree then_tree;
                then_tree.SetOpcode(cond.GetOpcode());
                then_tree.AddParam(cond.GetParam(1));
                then_tree.AddParam(GetParam(1));
                then_tree.AddParam(GetParam(2));
                then_tree.Rehash();
                CodeTree else_tree;
                else_tree.SetOpcode(cond.GetOpcode());
                else_tree.AddParam(cond.GetParam(2));
                else_tree.AddParam(GetParam(1));
                else_tree.AddParam(GetParam(2));
                else_tree.Rehash();
                SetOpcode(cond.GetOpcode());
                SetParam(0, cond.GetParam(0));
                SetParamMove(1, then_tree);
                SetParamMove(2, else_tree);
                return true; // rerun cIf optimization
            }
        }
        if(GetParam(1).GetOpcode() == GetParam(2).GetOpcode()
        && (GetParam(1).GetOpcode() == cIf
         || GetParam(1).GetOpcode() == cAbsIf))
        {
            CodeTree& leaf1 = GetParam(1);
            CodeTree& leaf2 = GetParam(2);
            if(leaf1.GetParam(0).IsIdenticalTo(leaf2.GetParam(0))
            && (leaf1.GetParam(1).IsIdenticalTo(leaf2.GetParam(1))
             || leaf1.GetParam(2).IsIdenticalTo(leaf2.GetParam(2))))
            {
            //     if(x, if(y,a,b), if(y,c,d))
            // ->  if(y, if(x,a,c), if(x,b,d))
            // when either a,c are identical or b,d are identical
                CodeTree then_tree;
                then_tree.SetOpcode(GetOpcode());
                then_tree.AddParam(GetParam(0));
                then_tree.AddParam(leaf1.GetParam(1));
                then_tree.AddParam(leaf2.GetParam(1));
                then_tree.Rehash();
                CodeTree else_tree;
                else_tree.SetOpcode(GetOpcode());
                else_tree.AddParam(GetParam(0));
                else_tree.AddParam(leaf1.GetParam(2));
                else_tree.AddParam(leaf2.GetParam(2));
                else_tree.Rehash();
                SetOpcode(leaf1.GetOpcode());
                SetParam(0, leaf1.GetParam(0));
                SetParamMove(1, then_tree);
                SetParamMove(2, else_tree);
                return true; // rerun cIf optimization
            // cIf [x (cIf [y a z]) (cIf [y z b])] : (cXor x y) z (cIf[x a b])
            // ^ if only we had cXor opcode.
            }
            if(leaf1.GetParam(1).IsIdenticalTo(leaf2.GetParam(1))
            && leaf1.GetParam(2).IsIdenticalTo(leaf2.GetParam(2)))
            {
                //    if(x, if(y,a,b), if(z,a,b))
                // -> if( if(x, y,z), a,b)
                CodeTree cond_tree;
                cond_tree.SetOpcode(GetOpcode());
                cond_tree.AddParamMove(GetParam(0));
                cond_tree.AddParam(leaf1.GetParam(0));
                cond_tree.AddParam(leaf2.GetParam(0));
                cond_tree.Rehash();
                SetOpcode(leaf1.GetOpcode());
                SetParamMove(0, cond_tree);
                SetParam(2, leaf1.GetParam(2));
                SetParam(1, leaf1.GetParam(1));
                return true; // rerun cIf optimization
            }
            if(leaf1.GetParam(1).IsIdenticalTo(leaf2.GetParam(2))
            && leaf1.GetParam(2).IsIdenticalTo(leaf2.GetParam(1)))
            {
                //    if(x, if(y,a,b), if(z,b,a))
                // -> if( if(x, y,!z), a,b)
                CodeTree not_tree;
                not_tree.SetOpcode(leaf2.GetOpcode() == cIf ? cNot : cAbsNot);
                not_tree.AddParam(leaf2.GetParam(0));
                not_tree.Rehash();
                CodeTree cond_tree;
                cond_tree.SetOpcode(GetOpcode());
                cond_tree.AddParamMove(GetParam(0));
                cond_tree.AddParam(leaf1.GetParam(0));
                cond_tree.AddParamMove(not_tree);
                cond_tree.Rehash();
                SetOpcode(leaf1.GetOpcode());
                SetParamMove(0, cond_tree);
                SetParam(2, leaf1.GetParam(2));
                SetParam(1, leaf1.GetParam(1));
                return true; // rerun cIf optimization
            }
        }

        // If the sub-expression evaluates to approx. zero, yield param3.
        // If the sub-expression evaluates to approx. nonzero, yield param2.
        MinMaxTree p = GetParam(0).CalculateResultBoundaries();
        switch(GetLogicalValue(p, GetOpcode()==cAbsIf))
        {
            case 1: // true
                Become(GetParam(1));
                return true; // rerun optimization (opcode changed)
            case 0: // false
                Become(GetParam(2));
                return true; // rerun optimization (opcode changed)
            default: ;
        }

        CodeTree& branch1 = GetParam(1);
        CodeTree& branch2 = GetParam(2);

        if(branch1.IsIdenticalTo(branch2))
        {
            // If both branches of an If() are identical, the test becomes unnecessary
            Become(GetParam(1));
            return true; // rerun optimization (opcode changed)
        }

        const OPCODE op1 = branch1.GetOpcode();
        const OPCODE op2 = branch2.GetOpcode();
        if(op1 == op2)
        {
            // If both branches apply the same unary function to different values,
            // extract the function. E.g. if(x,sin(a),sin(b)) -> sin(if(x,a,b))
            if(branch1.GetParamCount() == 1)
            {
                CodeTree changed_if;
                changed_if.SetOpcode(GetOpcode());
                changed_if.AddParamMove(GetParam(0));
                changed_if.AddParam(branch1.GetParam(0));
                changed_if.AddParam(branch2.GetParam(0));
                changed_if.Rehash();
                SetOpcode(op1);
                DelParams();
                AddParamMove(changed_if);
                return true; // rerun optimization (opcode changed)
            }
            if(op1 == cAdd    || op1 == cMul
            || op1 == cAnd    || op1 == cOr
            || op1 == cAbsAnd || op1 == cAbsOr
            || op1 == cMin    || op1 == cMax)
            {
                // If the two groups contain one or more
                // identical values, extract them.
                std::vector<CodeTree> overlap;
                for(size_t a=branch1.GetParamCount(); a-- > 0; )
                {
                    for(size_t b=branch2.GetParamCount(); b-- > 0; )
                    {
                        if(branch1.GetParam(a).IsIdenticalTo(branch2.GetParam(b)))
                        {
                            if(overlap.empty()) { branch1.CopyOnWrite(); branch2.CopyOnWrite(); }
                            overlap.push_back(branch1.GetParam(a));
                            branch2.DelParam(b);
                            branch1.DelParam(a);
                            break;
                        }
                    }
                }
                if(!overlap.empty())
                {
                    branch1.Rehash();
                    branch2.Rehash();
                    CodeTree changed_if;
                    changed_if.SetOpcode(GetOpcode());
                    changed_if.SetParamsMove(GetParams());
                    changed_if.Rehash();
                    SetOpcode(op1);
                    SetParamsMove(overlap);
                    AddParamMove(changed_if);
                    return true; // rerun optimization (opcode changed)
                }
            }
        }
        // if(x, y+z, y) -> if(x, z,0)+y
        if(op1 == cAdd
        || op1 == cMul
        || (op1 == cAnd && branch2.IsLogicalValue())
        || (op1 == cOr  && branch2.IsLogicalValue())
          )
        {
            for(size_t a=branch1.GetParamCount(); a-- > 0; )
                if(branch1.GetParam(a).IsIdenticalTo(branch2))
                {
                    branch1.CopyOnWrite();
                    branch1.DelParam(a);
                    branch1.Rehash();
                    CodeTree branch2_backup = branch2;
                    branch2 = CodeTree( (op1==cAdd||op1==cOr) ? 0.0 : 1.0 );
                    CodeTree changed_if;
                    changed_if.SetOpcode(GetOpcode());
                    changed_if.SetParamsMove(GetParams());
                    changed_if.Rehash();
                    SetOpcode(op1);
                    AddParamMove(branch2_backup);
                    AddParamMove(changed_if);
                    return true; // rerun optimization (opcode changed)
                }
        }
        // if(x, y&z, !!y) -> if(x, z, 1) & y
        if((op1 == cAnd || op1 == cOr) && op2 == cNotNot)
        {
            CodeTree& branch2op = branch2.GetParam(0);
            for(size_t a=branch1.GetParamCount(); a-- > 0; )
                if(branch1.GetParam(a).IsIdenticalTo(branch2op))
                {
                    branch1.CopyOnWrite();
                    branch1.DelParam(a);
                    branch1.Rehash();
                    CodeTree branch2_backup = branch2op;
                    branch2 = CodeTree( (op1==cOr) ? 0.0 : 1.0 );
                    CodeTree changed_if;
                    changed_if.SetOpcode(GetOpcode());
                    changed_if.SetParamsMove(GetParams());
                    changed_if.Rehash();
                    SetOpcode(op1);
                    AddParamMove(branch2_backup);
                    AddParamMove(changed_if);
                    return true; // rerun optimization (opcode changed)
                }
        }
        // if(x, y, y+z) -> if(x, 0,z)+y
        if(op2 == cAdd
        || op2 == cMul
        || (op2 == cAnd && branch1.IsLogicalValue())
        || (op2 == cOr  && branch1.IsLogicalValue())
          )
        {
            for(size_t a=branch2.GetParamCount(); a-- > 0; )
                if(branch2.GetParam(a).IsIdenticalTo(branch1))
                {
                    branch2.CopyOnWrite();
                    branch2.DelParam(a);
                    branch2.Rehash();
                    CodeTree branch1_backup = branch1;
                    branch1 = CodeTree( (op2==cAdd||op2==cOr) ? 0.0 : 1.0 );
                    CodeTree changed_if;
                    changed_if.SetOpcode(GetOpcode());
                    changed_if.SetParamsMove(GetParams());
                    changed_if.Rehash();
                    SetOpcode(op2);
                    AddParamMove(branch1_backup);
                    AddParamMove(changed_if);
                    return true; // rerun optimization (opcode changed)
                }
        }
        // if(x, !!y, y&z) -> if(x, 1, z) & y
        if((op2 == cAnd || op2 == cOr) && op1 == cNotNot)
        {
            CodeTree& branch1op = branch1.GetParam(0);
            for(size_t a=branch2.GetParamCount(); a-- > 0; )
                if(branch2.GetParam(a).IsIdenticalTo(branch1op))
                {
                    branch2.CopyOnWrite();
                    branch2.DelParam(a);
                    branch2.Rehash();
                    CodeTree branch1_backup = branch1op;
                    branch1 = CodeTree( (op2==cOr) ? 0.0 : 1.0 );
                    CodeTree changed_if;
                    changed_if.SetOpcode(GetOpcode());
                    changed_if.SetParamsMove(GetParams());
                    changed_if.Rehash();
                    SetOpcode(op2);
                    AddParamMove(branch1_backup);
                    AddParamMove(changed_if);
                    return true; // rerun optimization (opcode changed)
                }
        }
        return false; // No changes
    }

    bool CodeTree::ConstantFolding_PowOperations()
    {
        if(GetParam(0).IsImmed()
        && GetParam(1).IsImmed())
        {
            double const_value = fp_pow(GetParam(0).GetImmed(),
                                        GetParam(1).GetImmed());
            data = new CodeTreeData(const_value);
            return false;
        }
        if(GetParam(1).IsImmed()
        && (float)GetParam(1).GetImmed() == 1.0)
        {
            // Conversion through a float type value gets rid of
            // awkward abs(x)^1 generated from exp(log(x^6)/6),
            // without sacrificing as much precision as FloatEqual() does.
            // x^1 = x
            Become(GetParam(0));
            return true; // rerun optimization (opcode changed)
        }
        if(GetParam(0).IsImmed()
        && (float)GetParam(0).GetImmed() == 1.0)
        {
            // 1^x = 1
            data = new CodeTreeData(1.0);
            return false;
        }

        // 5^(20*x) = (5^20)^x
        if(GetParam(0).IsImmed()
        && GetParam(1).GetOpcode() == cMul)
        {
            bool changes = false;
            double base_immed = GetParam(0).GetImmed();
            CodeTree mulgroup = GetParam(1);
            for(size_t a=mulgroup.GetParamCount(); a-->0; )
                if(mulgroup.GetParam(a).IsImmed())
                {
                    double imm = mulgroup.GetParam(a).GetImmed();
                    //if(imm >= 0.0)
                    {
                        double new_base_immed = fp_pow(base_immed, imm);
                      if(std::isinf(new_base_immed) || new_base_immed == 0.0)
                        {
                            // It produced an infinity. Do not change.
                            break;
                        }

                        if(!changes)
                        {
                            changes = true;
                            mulgroup.CopyOnWrite();
                        }
                        base_immed = new_base_immed;
                        mulgroup.DelParam(a);
                        break; //
                    }
                }
            if(changes)
            {
                mulgroup.Rehash();
            #ifdef DEBUG_SUBSTITUTIONS
                std::cout << "Before pow-mul change: "; DumpTree(*this);
                std::cout << "\n";
            #endif
                GetParam(0).Become(CodeTree(base_immed));
                GetParam(1).Become(mulgroup);
            #ifdef DEBUG_SUBSTITUTIONS
                std::cout << "After pow-mul change: "; DumpTree(*this);
                std::cout << "\n";
            #endif
            }
        }
        // (x*20)^2 = x^2 * 20^2
        if(GetParam(1).IsImmed()
        && GetParam(0).GetOpcode() == cMul)
        {
            double exponent_immed = GetParam(1).GetImmed();
            double factor_immed   = 1.0;
            bool changes = false;
            CodeTree& mulgroup = GetParam(0);
            for(size_t a=mulgroup.GetParamCount(); a-->0; )
                if(mulgroup.GetParam(a).IsImmed())
                {
                    double imm = mulgroup.GetParam(a).GetImmed();
                    //if(imm >= 0.0)
                    {
                        double new_factor_immed = fp_pow(imm, exponent_immed);
                        if(std::isinf(new_factor_immed) || new_factor_immed == 0.0)
                        {
                            // It produced an infinity. Do not change.
                            break;
                        }
                        if(!changes)
                        {
                            changes = true;
                            mulgroup.CopyOnWrite();
                        }
                        factor_immed *= new_factor_immed;
                        mulgroup.DelParam(a);
                        break; //
                    }
                }
            if(changes)
            {
                mulgroup.Rehash();
                CodeTree newpow;
                newpow.SetOpcode(cPow);
                newpow.SetParamsMove(GetParams());
                SetOpcode(cMul);
                AddParamMove(newpow);
                AddParam( CodeTree(factor_immed) );
                return true; // rerun optimization (opcode changed)
            }
        }

        // (x^3)^2 = x^6
        // NOTE: If 3 is even and 3*2 is not, x must be changed to abs(x).
        if(GetParam(0).GetOpcode() == cPow
        && GetParam(1).IsImmed()
        && GetParam(0).GetParam(1).IsImmed())
        {
            double a = GetParam(0).GetParam(1).GetImmed();
            double b = GetParam(1).GetImmed();
            double c = a * b; // new exponent
            if(IsEvenIntegerConst(a) // a is an even int?
            && !IsEvenIntegerConst(c)) // c is not?
            {
                CodeTree newbase;
                newbase.SetOpcode(cAbs);
                newbase.AddParam(GetParam(0).GetParam(0));
                newbase.Rehash();
                SetParamMove(0, newbase);
            }
            else
                SetParam(0, GetParam(0).GetParam(0));
            SetParam(1, CodeTree(c));
        }
        return false; // No changes that require a rerun
    }

    bool CodeTree::ConstantFolding_ComparisonOperations()
    {
        static const RangeComparisonData Data[6] =
        {
            // cEqual:
            // Case:      p0 == p1  Antonym: p0 != p1
            // Synonym:   p1 == p0  Antonym: p1 != p0
            { RangeComparisonData::MakeTrue,  // If identical: always true
              {RangeComparisonData::MakeFalse,  // If Always p0 < p1: always false
               RangeComparisonData::Unchanged,
               RangeComparisonData::MakeFalse,  // If Always p0 > p1: always false
               RangeComparisonData::Unchanged},
             // NotNot(p0) if p1==1    NotNot(p1) if p0==1
             //    Not(p0) if p1==0       Not(p1) if p0==0
              {RangeComparisonData::MakeNotNotP0, RangeComparisonData::Eq1},
              {RangeComparisonData::MakeNotNotP1, RangeComparisonData::Eq1},
              {RangeComparisonData::MakeNotP0, RangeComparisonData::Eq0},
              {RangeComparisonData::MakeNotP1, RangeComparisonData::Eq0}
            },
            // cNEqual:
            // Case:      p0 != p1  Antonym: p0 == p1
            // Synonym:   p1 != p0  Antonym: p1 == p0
            { RangeComparisonData::MakeFalse,  // If identical: always false
              {RangeComparisonData::MakeTrue,  // If Always p0 < p1: always true
               RangeComparisonData::Unchanged,
               RangeComparisonData::MakeTrue,  // If Always p0 > p1: always true
               RangeComparisonData::Unchanged},
             // NotNot(p0) if p1==0    NotNot(p1) if p0==0
             //    Not(p0) if p1==1       Not(p1) if p0==1
              {RangeComparisonData::MakeNotNotP0, RangeComparisonData::Eq0},
              {RangeComparisonData::MakeNotNotP1, RangeComparisonData::Eq0},
              {RangeComparisonData::MakeNotP0, RangeComparisonData::Eq1},
              {RangeComparisonData::MakeNotP1, RangeComparisonData::Eq1}
            },
            // cLess:
            // Case:      p0 < p1   Antonym: p0 >= p1
            // Synonym:   p1 > p0   Antonym: p1 <= p0
            { RangeComparisonData::MakeFalse,  // If identical: always false
              {RangeComparisonData::MakeTrue,  // If Always p0  < p1: always true
               RangeComparisonData::MakeNEqual,
               RangeComparisonData::MakeFalse, // If Always p0 > p1: always false
               RangeComparisonData::MakeFalse},// If Always p0 >= p1: always false
             // Not(p0)   if p1>0 & p1<=1    --   NotNot(p1) if p0>=0 & p0<1
              {RangeComparisonData::MakeNotP0,    RangeComparisonData::Gt0Le1},
              {RangeComparisonData::MakeNotNotP1, RangeComparisonData::Ge0Lt1},
              {RangeComparisonData::Unchanged, RangeComparisonData::Never},
              {RangeComparisonData::Unchanged, RangeComparisonData::Never}
            },
            // cLessOrEq:
            // Case:      p0 <= p1  Antonym: p0 > p1
            // Synonym:   p1 >= p0  Antonym: p1 < p0
            { RangeComparisonData::MakeTrue,   // If identical: always true
              {RangeComparisonData::Unchanged, // If Always p0  < p1: ?
               RangeComparisonData::MakeTrue,  // If Always p0 <= p1: always true
               RangeComparisonData::MakeFalse, // If Always p0  > p1: always false
               RangeComparisonData::MakeEqual},// If Never  p0  < p1:  use cEqual
             // Not(p0)    if p1>=0 & p1<1   --   NotNot(p1) if p0>0 & p0<=1
              {RangeComparisonData::MakeNotP0,    RangeComparisonData::Ge0Lt1},
              {RangeComparisonData::MakeNotNotP1, RangeComparisonData::Gt0Le1},
              {RangeComparisonData::Unchanged, RangeComparisonData::Never},
              {RangeComparisonData::Unchanged, RangeComparisonData::Never}
            },
            // cGreater:
            // Case:      p0 >  p1  Antonym: p0 <= p1
            // Synonym:   p1 <  p0  Antonym: p1 >= p0
            { RangeComparisonData::MakeFalse,  // If identical: always false
              {RangeComparisonData::MakeFalse, // If Always p0  < p1: always false
               RangeComparisonData::MakeFalse, // If Always p0 <= p1: always false
               RangeComparisonData::MakeTrue,  // If Always p0  > p1: always true
               RangeComparisonData::MakeNEqual},
             // NotNot(p0) if p1>=0 & p1<1   --   Not(p1)   if p0>0 & p0<=1
              {RangeComparisonData::MakeNotNotP0, RangeComparisonData::Ge0Lt1},
              {RangeComparisonData::MakeNotP1,    RangeComparisonData::Gt0Le1},
              {RangeComparisonData::Unchanged, RangeComparisonData::Never},
              {RangeComparisonData::Unchanged, RangeComparisonData::Never}
            },
            // cGreaterOrEq:
            // Case:      p0 >= p1  Antonym: p0 < p1
            // Synonym:   p1 <= p0  Antonym: p1 > p0
            { RangeComparisonData::MakeTrue,   // If identical: always true
              {RangeComparisonData::MakeFalse, // If Always p0  < p1: always false
               RangeComparisonData::MakeEqual, // If Always p0 >= p1: always true
               RangeComparisonData::Unchanged, // If always p0  > p1: ?
               RangeComparisonData::MakeTrue}, // If Never  p0  > p1:  use cEqual
             // NotNot(p0) if p1>0 & p1<=1   --   Not(p1)    if p0>=0 & p0<1
              {RangeComparisonData::MakeNotNotP0, RangeComparisonData::Gt0Le1},
              {RangeComparisonData::MakeNotP1,    RangeComparisonData::Ge0Lt1},
              {RangeComparisonData::Unchanged, RangeComparisonData::Never},
              {RangeComparisonData::Unchanged, RangeComparisonData::Never}
            }
        };
        switch(Data[GetOpcode()-cEqual].Analyze(GetParam(0), GetParam(1)))
        {
            case RangeComparisonData::MakeFalse:
                data = new CodeTreeData(0.0); return true;
            case RangeComparisonData::MakeTrue:
                data = new CodeTreeData(1.0); return true;
            case RangeComparisonData::MakeEqual:  SetOpcode(cEqual); return true;
            case RangeComparisonData::MakeNEqual: SetOpcode(cNEqual); return true;
            case RangeComparisonData::MakeNotNotP0: SetOpcode(cNotNot); DelParam(1); return true;
            case RangeComparisonData::MakeNotNotP1: SetOpcode(cNotNot); DelParam(0); return true;
            case RangeComparisonData::MakeNotP0: SetOpcode(cNot); DelParam(1); return true;
            case RangeComparisonData::MakeNotP1: SetOpcode(cNot); DelParam(0); return true;
            case RangeComparisonData::Unchanged:;
        }
        return false;
    }

    bool CodeTree::ConstantFolding_Assimilate()
    {
        /* If the list contains another list of the same kind, assimilate it */
        bool assimilated = false;
        for(size_t a=GetParamCount(); a-- > 0; )
            if(GetParam(a).GetOpcode() == GetOpcode())
            {
              #ifdef DEBUG_SUBSTITUTIONS
                if(!assimilated)
                {
                    std::cout << "Before assimilation: "; DumpTree(*this);
                    std::cout << "\n";
                    assimilated = true;
                }
              #endif
                // Assimilate its children and remove it
                AddParamsMove(GetParam(a).GetUniqueRef().GetParams(), a);
            }
      #ifdef DEBUG_SUBSTITUTIONS
        if(assimilated)
        {
            std::cout << "After assimilation:   "; DumpTree(*this);
            std::cout << "\n";
        }
      #endif
        return assimilated;
    }

    void CodeTree::ConstantFolding()
    {
    #ifdef DEBUG_SUBSTITUTIONS
        std::cout << "Runs ConstantFolding for: "; DumpTree(*this);
        std::cout << "\n";
    #endif
        using namespace std;
    redo:;

        // Insert here any hardcoded constant-folding optimizations
        // that you want to be done whenever a new subtree is generated.
        /* Not recursive. */

        double const_value = 1.0;
        if(GetOpcode() != cImmed)
        {
            MinMaxTree p = CalculateResultBoundaries();
            if(p.has_min && p.has_max && p.min == p.max)
            {
                // Replace us with this immed
                const_value = p.min;
                goto ReplaceTreeWithConstValue;
            }
        }

        if(false)
        {
            ReplaceTreeWithOne:
                const_value = 1.0;
                goto ReplaceTreeWithConstValue;
            ReplaceTreeWithZero:
                const_value = 0.0;
            ReplaceTreeWithConstValue:
              #ifdef DEBUG_SUBSTITUTIONS
                std::cout << "Replacing "; DumpTree(*this);
                if(IsImmed())
                    std::cout << "(" << std::hex
                              << *(const uint_least64_t*)&GetImmed()
                              << std::dec << ")";
                std::cout << " with const value " << const_value;
                std::cout << "(" << std::hex
                          << *(const uint_least64_t*)&const_value
                          << std::dec << ")";
                std::cout << "\n";
              #endif
                data = new CodeTreeData(const_value);
                return;
            ReplaceTreeWithParam0:
              #ifdef DEBUG_SUBSTITUTIONS
                std::cout << "Before replace: "; DumpTree(*this);
                std::cout << "\n";
              #endif
                Become(GetParam(0));
              #ifdef DEBUG_SUBSTITUTIONS
                std::cout << "After replace: "; DumpTree(*this);
                std::cout << "\n";
              #endif
                goto redo;
        }

        /* Constant folding */
        switch(GetOpcode())
        {
            case cImmed:
                break; // nothing to do
            case VarBegin:
                break; // nothing to do

            case cAnd:
            case cAbsAnd:
            {
                ConstantFolding_Assimilate();
                for(size_t a=GetParamCount(); a-- > 0; )
                    switch(GetLogicalValue(GetParam(a).CalculateResultBoundaries(),
                                           GetOpcode()==cAbsAnd))
                    {
                        case 0: goto ReplaceTreeWithZero;
                        case 1: DelParam(a); break; // x & y & 1 = x & y;  x & 1 = !!x
                        default: ;
                    }
                switch(GetParamCount())
                {
                    case 0: goto ReplaceTreeWithOne;
                    case 1: SetOpcode(GetOpcode()==cAnd ? cNotNot : cAbsNotNot); goto redo; // Replace self with the single operand
                    default: if(GetOpcode()==cAnd) if(ConstantFolding_AndLogic()) goto redo;
                }
                break;
            }
            case cOr:
            case cAbsOr:
            {
                ConstantFolding_Assimilate();
                for(size_t a=GetParamCount(); a-- > 0; )
                    switch(GetLogicalValue(GetParam(a).CalculateResultBoundaries(),
                                           GetOpcode()==cAbsOr))
                    {
                        case 1: goto ReplaceTreeWithOne;
                        case 0: DelParam(a); break;
                        default: ;
                    }
                switch(GetParamCount())
                {
                    case 0: goto ReplaceTreeWithZero;
                    case 1: SetOpcode(GetOpcode()==cOr ? cNotNot : cAbsNotNot); goto redo; // Replace self with the single operand
                    default: if(GetOpcode()==cOr) if(ConstantFolding_OrLogic()) goto redo;
                }
                break;
            }
            case cNot:
            case cAbsNot:
            {
                unsigned opposite = 0;
                switch(GetParam(0).GetOpcode())
                {
                    case cEqual:       opposite = cNEqual; break;
                    case cNEqual:      opposite = cEqual; break;
                    case cLess:        opposite = cGreaterOrEq; break;
                    case cGreater:     opposite = cLessOrEq; break;
                    case cLessOrEq:    opposite = cGreater; break;
                    case cGreaterOrEq: opposite = cLess; break;
                    //cNotNot already handled by grammar: @L cNotNot
                    case cNot:         opposite = cNotNot; break;
                    case cAbsNot:      opposite = cAbsNotNot; break;
                    case cAbsNotNot:   opposite = cAbsNot; break;
                    default: break;
                }
                if(opposite)
                {
                    SetOpcode(OPCODE(opposite));
                    SetParamsMove(GetParam(0).GetUniqueRef().GetParams());
                    goto redo;
                }

                // If the sub-expression evaluates to approx. zero, yield one.
                // If the sub-expression evaluates to approx. nonzero, yield zero.
                switch(GetLogicalValue(GetParam(0).CalculateResultBoundaries(),
                                       GetOpcode()==cAbsNot))
                {
                    case 1: goto ReplaceTreeWithZero;
                    case 0: goto ReplaceTreeWithOne;
                    default: ;
                }
                if(GetOpcode() == cNot && GetParam(0).IsAlwaysSigned(true))
                    SetOpcode(cAbsNot);

                if(GetParam(0).GetOpcode() == cIf
                || GetParam(0).GetOpcode() == cAbsIf)
                {
                    CodeTree iftree = GetParam(0);
                    const CodeTree& ifp1 = iftree.GetParam(1);
                    const CodeTree& ifp2 = iftree.GetParam(2);
                    if(ifp1.GetOpcode() == cNot
                    || ifp1.GetOpcode() == cAbsNot)
                    {
                        // cNot [(cIf [x (cNot[y]) z])] -> cIf [x (cNotNot[y]) (cNot[z])]
                        SetParam(0, iftree.GetParam(0)); // condition
                        CodeTree p1;
                        p1.SetOpcode(ifp1.GetOpcode()==cNot ? cNotNot : cAbsNotNot);
                        p1.AddParam(ifp1.GetParam(0));
                        p1.Rehash();
                        AddParamMove(p1);
                        CodeTree p2;
                        p2.SetOpcode(GetOpcode());
                        p2.AddParam(ifp2);
                        p2.Rehash();
                        AddParamMove(p2);
                        SetOpcode(iftree.GetOpcode());
                        goto redo;
                    }
                    if(ifp2.GetOpcode() == cNot
                    || ifp2.GetOpcode() == cAbsNot)
                    {
                        // cNot [(cIf [x y (cNot[z])])] -> cIf [x (cNot[y]) (cNotNot[z])]
                        SetParam(0, iftree.GetParam(0)); // condition
                        CodeTree p1;
                        p1.SetOpcode(GetOpcode());
                        p1.AddParam(ifp1);
                        p1.Rehash();
                        AddParamMove(p1);
                        CodeTree p2;
                        p2.SetOpcode(ifp2.GetOpcode()==cNot ? cNotNot : cAbsNotNot);
                        p2.AddParam(ifp2.GetParam(0));
                        p2.Rehash();
                        AddParamMove(p2);
                        SetOpcode(iftree.GetOpcode());
                        goto redo;
                    }
                }
                break;
            }
            case cNotNot:
            case cAbsNotNot:
            {
                // The function of cNotNot is to protect a logical value from
                // changing. If the parameter is already a logical value,
                // then the cNotNot opcode is redundant.
                if(GetParam(0).IsLogicalValue())
                    goto ReplaceTreeWithParam0;

                // If the sub-expression evaluates to approx. zero, yield zero.
                // If the sub-expression evaluates to approx. nonzero, yield one.
                switch(GetLogicalValue(GetParam(0).CalculateResultBoundaries(),
                                       GetOpcode()==cAbsNotNot))
                {
                    case 0: goto ReplaceTreeWithZero;
                    case 1: goto ReplaceTreeWithOne;
                    default: ;
                }
                if(GetOpcode() == cNotNot && GetParam(0).IsAlwaysSigned(true))
                    SetOpcode(cAbsNotNot);

                if(GetParam(0).GetOpcode() == cIf
                || GetParam(0).GetOpcode() == cAbsIf)
                {
                    CodeTree iftree = GetParam(0);
                    const CodeTree& ifp1 = iftree.GetParam(1);
                    const CodeTree& ifp2 = iftree.GetParam(2);
                    if(ifp1.GetOpcode() == cNot
                    || ifp1.GetOpcode() == cAbsNot)
                    {
                        // cNotNot [(cIf [x (cNot[y]) z])] -> cIf [x (cNot[y]) (cNotNot[z])]
                        SetParam(0, iftree.GetParam(0)); // condition
                        AddParam(ifp1);
                        CodeTree p2;
                        p2.SetOpcode(GetOpcode());
                        p2.AddParam(ifp2);
                        p2.Rehash();
                        AddParamMove(p2);
                        SetOpcode(iftree.GetOpcode());
                        goto redo;
                    }
                    if(ifp2.GetOpcode() == cNot
                    || ifp2.GetOpcode() == cAbsNot)
                    {
                        // cNotNot [(cIf [x y (cNot[z])])] -> cIf [x (cNotNot[y]) (cNot[z])]
                        SetParam(0, iftree.GetParam(0)); // condition
                        CodeTree p1;
                        p1.SetOpcode(GetOpcode());
                        p1.AddParam(ifp1);
                        p1.Rehash();
                        AddParamMove(p1);
                        AddParam(ifp2);
                        SetOpcode(iftree.GetOpcode());
                        goto redo;
                    }
                }
                break;
            }
            case cIf:
            case cAbsIf:
            {
                if(ConstantFolding_IfOperations())
                    goto redo;
                break;
            }
            case cMul:
            {
            NowWeAreMulGroup: ;
                ConstantFolding_Assimilate();
                // If one sub-expression evalutes to exact zero, yield zero.
                double immed_product = 1.0;
                size_t n_immeds = 0; bool needs_resynth=false;
                for(size_t a=0; a<GetParamCount(); ++a)
                {
                    if(!GetParam(a).IsImmed()) continue;
                    // ^ Only check constant values
                    double immed = GetParam(a).GetImmed();
                    if(immed == 0.0) goto ReplaceTreeWithZero;
                    immed_product *= immed; ++n_immeds;
                }
                // Merge immeds.
                if(n_immeds > 1 || (n_immeds == 1 && FloatEqual(immed_product, 1.0)))
                    needs_resynth = true;
                if(needs_resynth)
                {
                    // delete immeds and add new ones
                #ifdef DEBUG_SUBSTITUTIONS
                    std::cout << "cMul: Will add new immed " << immed_product << "\n";
                #endif
                    for(size_t a=GetParamCount(); a-->0; )
                        if(GetParam(a).IsImmed())
                        {
                        #ifdef DEBUG_SUBSTITUTIONS
                            std::cout << " - For that, deleting immed " << GetParam(a).GetImmed();
                            std::cout << "\n";
                        #endif
                            DelParam(a);
                        }
                    if(!FloatEqual(immed_product, 1.0))
                        AddParam( CodeTree(immed_product) );
                }
                switch(GetParamCount())
                {
                    case 0: goto ReplaceTreeWithOne;
                    case 1: goto ReplaceTreeWithParam0; // Replace self with the single operand
                    default:
                        if(ConstantFolding_MulGrouping()) goto redo;
                        if(ConstantFolding_MulLogicItems()) goto redo;
                }
                break;
            }
            case cAdd:
            {
                ConstantFolding_Assimilate();
                double immed_sum = 0.0;
                size_t n_immeds = 0; bool needs_resynth=false;
                for(size_t a=0; a<GetParamCount(); ++a)
                {
                    if(!GetParam(a).IsImmed()) continue;
                    // ^ Only check constant values
                    double immed = GetParam(a).GetImmed();
                    immed_sum += immed; ++n_immeds;
                }
                // Merge immeds.
                if(n_immeds > 1 || (n_immeds == 1 && immed_sum == 0.0))
                    needs_resynth = true;
                if(needs_resynth)
                {
                    // delete immeds and add new ones
                #ifdef DEBUG_SUBSTITUTIONS
                    std::cout << "cAdd: Will add new immed " << immed_sum << "\n";
                    std::cout << "In: "; DumpTree(*this);
                    std::cout << "\n";
                #endif
                    for(size_t a=GetParamCount(); a-->0; )
                        if(GetParam(a).IsImmed())
                        {
                        #ifdef DEBUG_SUBSTITUTIONS
                            std::cout << " - For that, deleting immed " << GetParam(a).GetImmed();
                            std::cout << "\n";
                        #endif
                            DelParam(a);
                        }
                    if(!(immed_sum == 0.0))
                        AddParam( CodeTree(immed_sum) );
                }
                switch(GetParamCount())
                {
                    case 0: goto ReplaceTreeWithZero;
                    case 1: goto ReplaceTreeWithParam0; // Replace self with the single operand
                    default:
                        if(ConstantFolding_AddGrouping()) goto redo;
                        if(ConstantFolding_AddLogicItems()) goto redo;
                }
                break;
            }
            case cMin:
            {
                ConstantFolding_Assimilate();
                /* Goal: If there is any pair of two operands, where
                 * their ranges form a disconnected set, i.e. as below:
                 *     xxxxx
                 *            yyyyyy
                 * Then remove the larger one.
                 *
                 * Algorithm: 1. figure out the smallest maximum of all operands.
                 *            2. eliminate all operands where their minimum is
                 *               larger than the selected maximum.
                 */
                size_t preserve=0;
                MinMaxTree smallest_maximum;
                for(size_t a=0; a<GetParamCount(); ++a)
                {
                    MinMaxTree p = GetParam(a).CalculateResultBoundaries();
                    if(p.has_max && (!smallest_maximum.has_max || p.max < smallest_maximum.max))
                    {
                        smallest_maximum.max = p.max;
                        smallest_maximum.has_max = true;
                        preserve=a;
                }   }
                if(smallest_maximum.has_max)
                    for(size_t a=GetParamCount(); a-- > 0; )
                    {
                        MinMaxTree p = GetParam(a).CalculateResultBoundaries();
                        if(p.has_min && a != preserve && p.min >= smallest_maximum.max)
                            DelParam(a);
                    }
                //fprintf(stderr, "Remains: %u\n", (unsigned)GetParamCount());
                if(GetParamCount() == 1)
                {
                    // Replace self with the single operand
                    goto ReplaceTreeWithParam0;
                }
                break;
            }
            case cMax:
            {
                ConstantFolding_Assimilate();
                /* Goal: If there is any pair of two operands, where
                 * their ranges form a disconnected set, i.e. as below:
                 *     xxxxx
                 *            yyyyyy
                 * Then remove the smaller one.
                 *
                 * Algorithm: 1. figure out the biggest minimum of all operands.
                 *            2. eliminate all operands where their maximum is
                 *               smaller than the selected minimum.
                 */
                size_t preserve=0;
                MinMaxTree biggest_minimum;
                for(size_t a=0; a<GetParamCount(); ++a)
                {
                    MinMaxTree p = GetParam(a).CalculateResultBoundaries();
                    if(p.has_min && (!biggest_minimum.has_min || p.min > biggest_minimum.min))
                    {
                        biggest_minimum.min = p.min;
                        biggest_minimum.has_min = true;
                        preserve=a;
                }   }
                if(biggest_minimum.has_min)
                {
                    //fprintf(stderr, "Removing all where max < %g\n", biggest_minimum.min);
                    for(size_t a=GetParamCount(); a-- > 0; )
                    {
                        MinMaxTree p = GetParam(a).CalculateResultBoundaries();
                        if(p.has_max && a != preserve && p.max < biggest_minimum.min)
                        {
                            //fprintf(stderr, "Removing %g\n", p.max);
                            DelParam(a);
                        }
                    }
                }
                //fprintf(stderr, "Remains: %u\n", (unsigned)GetParamCount());
                if(GetParamCount() == 1)
                {
                    // Replace self with the single operand
                    goto ReplaceTreeWithParam0;
                }
                break;
            }

            case cEqual:
            case cNEqual:
            case cLess:
            case cGreater:
            case cLessOrEq:
            case cGreaterOrEq:
                if(ConstantFolding_ComparisonOperations()) goto redo;
                // Any reversible functions:
                //   sin(x)  -> ASIN: Not doable, x can be cyclic
                //   asin(x) -> SIN: doable.
                //                   Invalid combinations are caught by
                //                   range-estimation. Threshold is at |pi/2|.
                //   acos(x) -> COS: doable.
                //                   Invalid combinations are caught by
                //                   range-estimation. Note that though
                //                   the range is contiguous, it is direction-flipped.
                //    log(x) -> EXP: no problem
                //   exp2, exp10: Converted to cPow, done by grammar.
                //   atan(x) -> TAN: doable.
                //                   Invalid combinations are caught by
                //                   range-estimation. Threshold is at |pi/2|.
                //   sinh(x) -> ASINH: no problem
                //   tanh(x) -> ATANH: no problem, but atanh is limited to -1..1
                //                     Invalid combinations are caught by
                //                     range-estimation, but the exact value
                //                     of 1.0 still needs checking, because
                //                     it involves infinity.
                if(GetParam(1).IsImmed())
                    switch(GetParam(0).GetOpcode())
                    {
                        case cAsin:
                            SetParam(0, GetParam(0).GetParam(0));
                            SetParam(1, CodeTree(fp_sin(GetParam(1).GetImmed())));
                            goto redo;
                        case cAcos:
                            // -1..+1 --> pi..0 (polarity-flipping)
                            SetParam(0, GetParam(0).GetParam(0));
                            SetParam(1, CodeTree(fp_cos(GetParam(1).GetImmed())));
                            SetOpcode( GetOpcode()==cLess ? cGreater
                                     : GetOpcode()==cLessOrEq ? cGreaterOrEq
                                     : GetOpcode()==cGreater ? cLess
                                     : GetOpcode()==cGreaterOrEq ? cLessOrEq
                                     : GetOpcode() );
                            goto redo;
                        case cAtan:
                            SetParam(0, GetParam(0).GetParam(0));
                            SetParam(1, CodeTree(fp_tan(GetParam(1).GetImmed())));
                            goto redo;
                        case cLog:
                            // Different logarithms have a constant-multiplication,
                            // which is no problem.
                            SetParam(0, GetParam(0).GetParam(0));
                            SetParam(1, CodeTree(fp_exp(GetParam(1).GetImmed())));
                            goto redo;
                        case cSinh:
                            SetParam(0, GetParam(0).GetParam(0));
                            SetParam(1, CodeTree(fp_asinh(GetParam(1).GetImmed())));
                            goto redo;
                        case cTanh:
                            if(fabs(GetParam(1).GetImmed()) < 1.0)
                            {
                                SetParam(0, GetParam(0).GetParam(0));
                                SetParam(1, CodeTree(fp_atanh(GetParam(1).GetImmed())));
                                goto redo;
                            }
                            break;
                        default: break;
                    }
                break;

            case cAbs:
            {
                /* If we know the operand is always positive, cAbs is redundant.
                 * If we know the operand is always negative, use actual negation.
                 */
                MinMaxTree p0 = GetParam(0).CalculateResultBoundaries();
                if(p0.has_min && p0.min >= 0.0)
                    goto ReplaceTreeWithParam0;
                if(p0.has_max && p0.max <= NEGATIVE_MAXIMUM)
                {
                    /* abs(negative) = negative*-1 */
                    SetOpcode(cMul);
                    AddParam( CodeTree(-1.0) );
                    /* The caller of ConstantFolding() will do Sort() and Rehash() next.
                     * Thus, no need to do it here. */
                    /* We were changed into a cMul group. Do cMul folding. */
                    goto NowWeAreMulGroup;
                }
                /* If the operand is a cMul group, find elements
                 * that are always positive and always negative,
                 * and move them out, e.g. abs(p*n*x*y) = p*(-n)*abs(x*y)
                 */
                if(GetParam(0).GetOpcode() == cMul)
                {
                    const CodeTree& p = GetParam(0);
                    std::vector<CodeTree> pos_set;
                    std::vector<CodeTree> neg_set;
                    for(size_t a=0; a<p.GetParamCount(); ++a)
                    {
                        p0 = p.GetParam(a).CalculateResultBoundaries();
                        if(p0.has_min && p0.min >= 0.0)
                            { pos_set.push_back(p.GetParam(a)); }
                        if(p0.has_max && p0.max <= NEGATIVE_MAXIMUM)
                            { neg_set.push_back(p.GetParam(a)); }
                    }
                #ifdef DEBUG_SUBSTITUTIONS
                    std::cout << "Abs: mul group has " << pos_set.size()
                              << " pos, " << neg_set.size() << "neg\n";
                #endif
                    if(!pos_set.empty() || !neg_set.empty())
                    {
                #ifdef DEBUG_SUBSTITUTIONS
                        std::cout << "AbsReplace-Before: ";
                        DumpTree(*this);
                        std::cout << "\n" << std::flush;
                        DumpHashes(*this, std::cout);
                #endif
                        CodeTree pclone;
                        pclone.SetOpcode(cMul);
                        for(size_t a=0; a<p.GetParamCount(); ++a)
                        {
                            p0 = p.GetParam(a).CalculateResultBoundaries();
                            if((p0.has_min && p0.min >= 0.0)
                            || (p0.has_max && p0.max <= NEGATIVE_MAXIMUM))
                                {/*pclone.DelParam(a);*/}
                            else
                                pclone.AddParam( p.GetParam(a) );
                            /* Here, p*n*x*y -> x*y.
                             * p is saved in pos_set[]
                             * n is saved in neg_set[]
                             */
                        }
                        pclone.Rehash();
                        CodeTree abs_mul;
                        abs_mul.SetOpcode(cAbs);
                        abs_mul.AddParamMove(pclone);
                        abs_mul.Rehash();
                        CodeTree mulgroup;
                        mulgroup.SetOpcode(cMul);
                        mulgroup.AddParamMove(abs_mul); // cAbs[whatever remains in p]
                        mulgroup.AddParamsMove(pos_set);
                        /* Now:
                         * mulgroup  = p * Abs(x*y)
                         */
                        if(!neg_set.empty())
                        {
                            if(neg_set.size() % 2)
                                mulgroup.AddParam( CodeTree(-1.0) );
                            mulgroup.AddParamsMove(neg_set);
                            /* Now:
                             * mulgroup = p * n * -1 * Abs(x*y)
                             */
                        }
                        Become(mulgroup);
                #ifdef DEBUG_SUBSTITUTIONS
                        std::cout << "AbsReplace-After: ";
                        DumpTree(*this, std::cout);
                        std::cout << "\n" << std::flush;
                        DumpHashes(*this, std::cout);
                #endif
                        /* We were changed into a cMul group. Do cMul folding. */
                        goto NowWeAreMulGroup;
                    }
                }
                break;
            }

            #define HANDLE_UNARY_CONST_FUNC(funcname) \
                if(GetParam(0).IsImmed()) \
                    { const_value = funcname(GetParam(0).GetImmed()); \
                      goto ReplaceTreeWithConstValue; }

            case cLog:
                HANDLE_UNARY_CONST_FUNC(log);
                if(GetParam(0).GetOpcode() == cPow)
                {
                    CodeTree pow = GetParam(0);
                    if(pow.GetParam(0).IsAlwaysSigned(true))  // log(posi ^ y) = y*log(posi)
                    {
                        pow.CopyOnWrite();
                        pow.SetOpcode(cLog);
                        SetOpcode(cMul);
                        AddParamMove(pow.GetParam(1));
                        pow.DelParam(1);
                        pow.Rehash();
                        SetParamMove(0, pow);
                        goto NowWeAreMulGroup;
                    }
                    if(pow.GetParam(1).IsAlwaysParity(false)) // log(x ^ even) = even*log(abs(x))
                    {
                        pow.CopyOnWrite();
                        CodeTree abs;
                        abs.SetOpcode(cAbs);
                        abs.AddParamMove(pow.GetParam(0));
                        abs.Rehash();
                        pow.SetOpcode(cLog);
                        SetOpcode(cMul);
                        pow.SetParamMove(0, abs);
                        AddParamMove(pow.GetParam(1));
                        pow.DelParam(1);
                        pow.Rehash();
                        SetParamMove(0, pow);
                        goto NowWeAreMulGroup;
                    }
                }
                else if(GetParam(0).GetOpcode() == cAbs)
                {
                    // log(abs(x^y)) = y*log(abs(x))
                    CodeTree pow = GetParam(0).GetParam(0);
                    if(pow.GetOpcode() == cPow)
                    {
                        pow.CopyOnWrite();
                        CodeTree abs;
                        abs.SetOpcode(cAbs);
                        abs.AddParamMove(pow.GetParam(0));
                        abs.Rehash();
                        pow.SetOpcode(cLog);
                        SetOpcode(cMul);
                        pow.SetParamMove(0, abs);
                        AddParamMove(pow.GetParam(1));
                        pow.DelParam(1);
                        pow.Rehash();
                        SetParamMove(0, pow);
                        goto NowWeAreMulGroup;
                    }
                }
                break;
            case cAcosh: HANDLE_UNARY_CONST_FUNC(fp_acosh); break;
            case cAsinh: HANDLE_UNARY_CONST_FUNC(fp_asinh); break;
            case cAtanh: HANDLE_UNARY_CONST_FUNC(fp_atanh); break;
            case cAcos: HANDLE_UNARY_CONST_FUNC(fp_acos); break;
            case cAsin: HANDLE_UNARY_CONST_FUNC(fp_asin); break;
            case cAtan: HANDLE_UNARY_CONST_FUNC(fp_atan); break;
            case cCosh: HANDLE_UNARY_CONST_FUNC(fp_cosh); break;
            case cSinh: HANDLE_UNARY_CONST_FUNC(fp_sinh); break;
            case cTanh: HANDLE_UNARY_CONST_FUNC(fp_tanh); break;
            case cSin: HANDLE_UNARY_CONST_FUNC(fp_sin); break;
            case cCos: HANDLE_UNARY_CONST_FUNC(fp_cos); break;
            case cTan: HANDLE_UNARY_CONST_FUNC(fp_tan); break;
            case cCeil: HANDLE_UNARY_CONST_FUNC(fp_ceil); break;
            case cTrunc: HANDLE_UNARY_CONST_FUNC(fp_trunc); break;
            case cFloor: HANDLE_UNARY_CONST_FUNC(fp_floor); break;
            case cCbrt: HANDLE_UNARY_CONST_FUNC(fp_cbrt); break; // converted into cPow x 0.33333
            case cSqrt: HANDLE_UNARY_CONST_FUNC(fp_sqrt); break; // converted into cPow x 0.5
            case cExp: HANDLE_UNARY_CONST_FUNC(fp_exp); break; // convered into cPow CONSTANT_E x
            case cInt: HANDLE_UNARY_CONST_FUNC(fp_int); break;
            case cLog2: HANDLE_UNARY_CONST_FUNC(fp_log2); break;
            case cLog10: HANDLE_UNARY_CONST_FUNC(fp_log10); break;

            case cLog2by:
                if(GetParam(0).IsImmed()
                && GetParam(1).IsImmed())
                    { const_value = fp_log2(GetParam(0).GetImmed()) * GetParam(1).GetImmed();
                      goto ReplaceTreeWithConstValue; }
                break;

            case cMod: /* Can more be done than this? */
                if(GetParam(0).IsImmed()
                && GetParam(1).IsImmed())
                    { const_value = fmod(GetParam(0).GetImmed(), GetParam(1).GetImmed());
                      goto ReplaceTreeWithConstValue; }
                break;

            case cAtan2:
            {
                /* Range based optimizations for (y,x):
                 * If y is +0 and x <= -0, +pi is returned
                 * If y is -0 and x <= -0, -pi is returned (assumed never happening)
                 * If y is +0 and x >= +0, +0 is returned
                 * If y is -0 and x >= +0, -0 is returned  (assumed never happening)
                 * If x is +-0 and y < 0, -pi/2 is returned
                 * If x is +-0 and y > 0, +pi/2 is returned
                 * Otherwise, perform constant folding when available
                 * If we know x <> 0, convert into atan(y / x)
                 *   TODO: Figure out whether the above step is wise
                 *         It allows e.g. atan2(6*x, 3*y) -> atan(2*x/y)
                 *         when we know y != 0
                 */
                MinMaxTree p0 = GetParam(0).CalculateResultBoundaries();
                MinMaxTree p1 = GetParam(1).CalculateResultBoundaries();
                if(GetParam(0).IsImmed()
                && FloatEqual(GetParam(0).GetImmed(), 0.0))   // y == 0
                {
                    if(p1.has_max && p1.max < 0)              // y == 0 && x < 0
                        { const_value = CONSTANT_PI; goto ReplaceTreeWithConstValue; }
                    if(p1.has_min && p1.min >= 0.0)           // y == 0 && x >= 0.0
                        { const_value = 0.0; goto ReplaceTreeWithConstValue; }
                }
                if(GetParam(1).IsImmed()
                && FloatEqual(GetParam(1).GetImmed(), 0.0))   // x == 0
                {
                    if(p0.has_max && p0.max < 0)              // y < 0 && x == 0
                        { const_value = -CONSTANT_PIHALF; goto ReplaceTreeWithConstValue; }
                    if(p0.has_min && p0.min > 0)              // y > 0 && x == 0
                        { const_value =  CONSTANT_PIHALF; goto ReplaceTreeWithConstValue; }
                }
                if(GetParam(0).IsImmed()
                && GetParam(1).IsImmed())
                    { const_value = atan2(GetParam(0).GetImmed(),
                                          GetParam(1).GetImmed());
                      goto ReplaceTreeWithConstValue; }
                if((p1.has_min && p1.min > 0.0)               // p1 != 0.0
                || (p1.has_max && p1.max < NEGATIVE_MAXIMUM)) // become atan(p0 / p1)
                {
                    CodeTree pow_tree;
                    pow_tree.SetOpcode(cPow);
                    pow_tree.AddParamMove(GetParam(1));
                    pow_tree.AddParam(CodeTree(-1.0));
                    pow_tree.Rehash();
                    CodeTree div_tree;
                    div_tree.SetOpcode(cMul);
                    div_tree.AddParamMove(GetParam(0));
                    div_tree.AddParamMove(pow_tree);
                    div_tree.Rehash();
                    SetOpcode(cAtan);
                    SetParamMove(0, div_tree);
                    DelParam(1);
                }
                break;
            }

            case cPow:
            {
                if(ConstantFolding_PowOperations()) goto redo;
                break;
            }

            /* The following opcodes are processed by GenerateFrom()
             * within fpoptimizer_bytecode_to_codetree.c and thus
             * they will never occur in the calling context for the
             * most of the parsing context. They may however occur
             * at the late phase, so we deal with them.
             */
            case cDiv: // converted into cPow y -1
                if(GetParam(0).IsImmed()
                && GetParam(1).IsImmed()
                && GetParam(1).GetImmed() != 0.0)
                    { const_value = GetParam(0).GetImmed() / GetParam(1).GetImmed();
                      goto ReplaceTreeWithConstValue; }
                break;
            case cInv: // converted into cPow y -1
                if(GetParam(0).IsImmed()
                && GetParam(0).GetImmed() != 0.0)
                    { const_value = 1.0 / GetParam(0).GetImmed();
                      goto ReplaceTreeWithConstValue; }
                // Note: Could use (mulgroup)^immed optimization from cPow
                break;
            case cSub: // converted into cMul y -1
                if(GetParam(0).IsImmed()
                && GetParam(1).IsImmed())
                    { const_value = GetParam(0).GetImmed() - GetParam(1).GetImmed();
                      goto ReplaceTreeWithConstValue; }
                break;
            case cNeg: // converted into cMul x -1
                if(GetParam(0).IsImmed())
                    { const_value = -GetParam(0).GetImmed();
                      goto ReplaceTreeWithConstValue; }
                break;
            case cRad: // converted into cMul x CONSTANT_RD
                if(GetParam(0).IsImmed())
                    { const_value = GetParam(0).GetImmed() * CONSTANT_RD;
                      goto ReplaceTreeWithConstValue; }
                break;
            case cDeg: // converted into cMul x CONSTANT_DR
                if(GetParam(0).IsImmed())
                    { const_value = GetParam(0).GetImmed() * CONSTANT_DR;
                      goto ReplaceTreeWithConstValue; }
                break;
            case cSqr: // converted into cMul x x
                if(GetParam(0).IsImmed())
                    { const_value = GetParam(0).GetImmed() * GetParam(0).GetImmed();
                      goto ReplaceTreeWithConstValue; }
                break;
            case cExp2: // converted into cPow 2.0 x
                HANDLE_UNARY_CONST_FUNC(fp_exp2); break;
            case cRSqrt: // converted into cPow x -0.5
                if(GetParam(0).IsImmed())
                    { const_value = 1.0 / sqrt(GetParam(0).GetImmed());
                      goto ReplaceTreeWithConstValue; }
                break;
            case cCot: // converted into cMul (cPow (cTan x) -1)
                if(GetParam(0).IsImmed())
                    { double tmp = tan(GetParam(0).GetImmed());
                      if(tmp != 0.0)
                      { const_value = 1.0 / tmp;
                        goto ReplaceTreeWithConstValue; } }
                break;
            case cSec: // converted into cMul (cPow (cCos x) -1)
                if(GetParam(0).IsImmed())
                    { double tmp = cos(GetParam(0).GetImmed());
                      if(tmp != 0.0)
                      { const_value = 1.0 / tmp;
                        goto ReplaceTreeWithConstValue; } }
                break;
            case cCsc: // converted into cMul (cPow (cSin x) -1)
                if(GetParam(0).IsImmed())
                    { double tmp = sin(GetParam(0).GetImmed());
                      if(tmp != 0.0)
                      { const_value = 1.0 / tmp;
                        goto ReplaceTreeWithConstValue; } }
                break;

            /* Opcodes that do not occur in the tree for other reasons */
            case cRDiv: // version of cDiv
            case cRSub: // version of cSub
            case cDup:
            case cFetch:
            case cPopNMov:
            case cNop:
            case cJump:
                break; /* Should never occur */

            /* Opcodes that we can't do anything about */
            case cPCall:
            case cFCall:
            case cEval:
                break;
        }
    }
}

#endif

#line 1 "fpoptimizer/fpoptimizer_rangeestimation.c"
// line removed for fpoptimizer.c: #include "fpoptimizer_codetree.h"
// line removed for fpoptimizer.c: #include "fpoptimizer_consts.h"

#include <cmath> /* for CalculateResultBoundaries() */

#ifdef FP_SUPPORT_OPTIMIZER

using namespace FUNCTIONPARSERTYPES;
using namespace FPoptimizer_CodeTree;

//#define DEBUG_SUBSTITUTIONS_extra_verbose

namespace FPoptimizer_CodeTree
{
    MinMaxTree CodeTree::CalculateResultBoundaries() const
#ifdef DEBUG_SUBSTITUTIONS_extra_verbose
    {
        MinMaxTree tmp = CalculateResultBoundaries_do();
        std::cout << "Estimated boundaries: ";
        if(tmp.has_min) std::cout << tmp.min; else std::cout << "-inf";
        std::cout << " .. ";
        if(tmp.has_max) std::cout << tmp.max; else std::cout << "+inf";
        std::cout << ": ";
        FPoptimizer_CodeTree::DumpTree(*this);
        std::cout << std::endl;
        return tmp;
    }
    MinMaxTree CodeTree::CalculateResultBoundaries_do() const
#endif
    {
        using namespace std;
        switch( GetOpcode() )
        {
            case cImmed:
                return MinMaxTree(GetImmed(), GetImmed()); // a definite value.
            case cAnd:
            case cAbsAnd:
            case cOr:
            case cAbsOr:
            case cNot:
            case cAbsNot:
            case cNotNot:
            case cAbsNotNot:
            case cEqual:
            case cNEqual:
            case cLess:
            case cLessOrEq:
            case cGreater:
            case cGreaterOrEq:
            {
                /* These operations always produce truth values (0 or 1) */
                /* Narrowing them down is a matter of performing Constant optimization */
                return MinMaxTree( 0.0, 1.0 );
            }
            case cAbs:
            {
                /* cAbs always produces a positive value */
                MinMaxTree m = GetParam(0).CalculateResultBoundaries();
                if(m.has_min && m.has_max)
                {
                    if(m.min < 0.0 && m.max >= 0.0) // ex. -10..+6 or -6..+10
                    {
                        /* -x..+y: spans across zero. min=0, max=greater of |x| and |y|. */
                        double tmp = -m.min; if(tmp > m.max) m.max = tmp;
                        m.min = 0.0; m.has_min = true;
                    }
                    else if(m.min < 0.0) // ex. -10..-4
                        { double tmp = m.max; m.max = -m.min; m.min = -tmp; }
                }
                else if(!m.has_min && m.has_max && m.max < 0.0) // ex. -inf..-10
                {
                    m.min = fabs(m.max); m.has_min = true; m.has_max = false;
                }
                else if(!m.has_max && m.has_min && m.min > 0.0) // ex. +10..+inf
                {
                    m.min = fabs(m.min); m.has_min = true; m.has_max = false;
                }
                else // ex. -inf..+inf, -inf..+10, -10..+inf
                {
                    // all of these cover -inf..0, 0..+inf, or both
                    m.min = 0.0; m.has_min = true; m.has_max = false;
                }
                return m;
            }

            case cLog: /* Defined for 0.0 < x <= inf */
            {
                MinMaxTree m = GetParam(0).CalculateResultBoundaries();
                if(m.has_min) { if(m.min < 0.0) m.has_min = false; else m.min = log(m.min); } // No boundaries
                if(m.has_max) { if(m.max < 0.0) m.has_max = false; else m.max = log(m.max); }
                return m;
            }

            case cLog2: /* Defined for 0.0 < x <= inf */
            {
                MinMaxTree m = GetParam(0).CalculateResultBoundaries();
                if(m.has_min) { if(m.min < 0.0) m.has_min = false; else m.min = fp_log2(m.min); } // No boundaries
                if(m.has_max) { if(m.max < 0.0) m.has_max = false; else m.max = fp_log2(m.max); }
                return m;
            }

            case cLog10: /* Defined for 0.0 < x <= inf */
            {
                MinMaxTree m = GetParam(0).CalculateResultBoundaries();
                if(m.has_min) { if(m.min < 0.0) m.has_min = false; else m.min = fp_log10(m.min); }
                if(m.has_max) { if(m.max < 0.0) m.has_max = false; else m.max = fp_log10(m.max); }
                return m;
            }

            case cAcosh: /* defined for             1.0 <  x <= inf */
            {
                MinMaxTree m = GetParam(0).CalculateResultBoundaries();
                if(m.has_min) { if(m.min <= 1.0) m.has_min = false; else m.min = fp_acosh(m.min); } // No boundaries
                if(m.has_max) { if(m.max <= 1.0) m.has_max = false; else m.max = fp_acosh(m.max); }
                return m;
            }
            case cAsinh: /* defined for all values -inf <= x <= inf */
            {
                MinMaxTree m = GetParam(0).CalculateResultBoundaries();
                if(m.has_min) m.min = fp_asinh(m.min); // No boundaries
                if(m.has_max) m.max = fp_asinh(m.max);
                return m;
            }
            case cAtanh: /* defined for all values -inf <= x <= inf */
            {
                MinMaxTree m = GetParam(0).CalculateResultBoundaries();
                if(m.has_min) m.min = fp_atanh(m.min); // No boundaries
                if(m.has_max) m.max = fp_atanh(m.max);
                return m;
            }
            case cAcos: /* defined for -1.0 <= x < 1, results within CONSTANT_PI..0 */
            {
                /* Somewhat complicated to narrow down from this */
                /* TODO: A resourceful programmer may add it later. */
                return MinMaxTree( 0.0, CONSTANT_PI );
            }
            case cAsin: /* defined for -1.0 <= x < 1, results within -CONSTANT_PIHALF..CONSTANT_PIHALF */
            {
                /* Somewhat complicated to narrow down from this */
                /* TODO: A resourceful programmer may add it later. */
                return MinMaxTree( -CONSTANT_PIHALF, CONSTANT_PIHALF );
            }
            case cAtan: /* defined for all values -inf <= x <= inf */
            {
                MinMaxTree m = GetParam(0).CalculateResultBoundaries();
                if(m.has_min) m.min = atan(m.min); else { m.min = -CONSTANT_PIHALF; m.has_min = true; }
                if(m.has_max) m.max = atan(m.max); else { m.max =  CONSTANT_PIHALF; m.has_max = true; }
                return m;
            }
            case cAtan2: /* too complicated to estimate */
            {
                MinMaxTree p0 = GetParam(0).CalculateResultBoundaries();
                MinMaxTree p1 = GetParam(1).CalculateResultBoundaries();
                if(GetParam(0).IsImmed()
                && FloatEqual(GetParam(0).GetImmed(), 0.0))   // y == 0
                {
                    // Either 0.0 or CONSTANT_PI
                    return MinMaxTree(0.0, CONSTANT_PI);
                }
                if(GetParam(1).IsImmed()
                && FloatEqual(GetParam(1).GetImmed(), 0.0))   // x == 0
                {
                    // EIther -CONSTANT_PIHALF or +CONSTANT_PIHALF
                    return MinMaxTree(-CONSTANT_PIHALF, CONSTANT_PIHALF);
                }
                // Anything else
                /* Somewhat complicated to narrow down from this */
                /* TODO: A resourceful programmer may add it later. */
                return MinMaxTree(-CONSTANT_PI, CONSTANT_PI);
            }

            case cSin:
            case cCos:
            {
                /* Could be narrowed down from here,
                 * but it's too complicated due to
                 * the cyclic nature of the function. */
                /* TODO: A resourceful programmer may add it later. */
                return MinMaxTree(-1.0, 1.0);
            }
            case cTan:
            {
                /* Could be narrowed down from here,
                 * but it's too complicated due to
                 * the cyclic nature of the function */
                /* TODO: A resourceful programmer may add it later. */
                return MinMaxTree(); // (CONSTANT_NEG_INF, CONSTANT_POS_INF);
            }

            case cCeil:
            {
                MinMaxTree m = GetParam(0).CalculateResultBoundaries();
                m.max = std::ceil(m.max); // ceil() may increase the value, may not decrease
                return m;
            }
            case cFloor:
            {
                MinMaxTree m = GetParam(0).CalculateResultBoundaries();
                m.min = std::floor(m.min); // floor() may decrease the value, may not increase
                return m;
            }
            case cTrunc:
            {
                MinMaxTree m = GetParam(0).CalculateResultBoundaries();
                m.min = std::floor(m.min); // trunc() may either increase or decrease the value
                m.max = std::ceil(m.max); // for safety, we assume both
                return m;
            }
            case cInt:
            {
                MinMaxTree m = GetParam(0).CalculateResultBoundaries();
                m.min = std::floor(m.min); // int() may either increase or decrease the value
                m.max = std::ceil(m.max); // for safety, we assume both
                return m;
            }
            case cSinh: /* defined for all values -inf <= x <= inf */
            {
                MinMaxTree m = GetParam(0).CalculateResultBoundaries();
                if(m.has_min) m.min = sinh(m.min); // No boundaries
                if(m.has_max) m.max = sinh(m.max);
                return m;
            }
            case cTanh: /* defined for all values -inf <= x <= inf, results within -1..1 */
            {
                MinMaxTree m = GetParam(0).CalculateResultBoundaries();
                if(m.has_min) m.min = tanh(m.min); else { m.has_min = true; m.min =-1.0; }
                if(m.has_max) m.max = tanh(m.max); else { m.has_max = true; m.max = 1.0; }
                return m;
            }
            case cCosh: /* defined for all values -inf <= x <= inf, results within 1..inf */
            {
                MinMaxTree m = GetParam(0).CalculateResultBoundaries();
                if(m.has_min)
                {
                    if(m.has_max) // max, min
                    {
                        if(m.min >= 0.0 && m.max >= 0.0) // +x .. +y
                            { m.min = cosh(m.min); m.max = cosh(m.max); }
                        else if(m.min < 0.0 && m.max >= 0.0) // -x .. +y
                            { double tmp = cosh(m.min); m.max = cosh(m.max);
                              if(tmp > m.max) m.max = tmp;
                              m.min = 1.0; }
                        else // -x .. -y
                            { m.min = cosh(m.min); m.max = cosh(m.max);
                              std::swap(m.min, m.max); }
                    }
                    else // min, no max
                    {
                        if(m.min >= 0.0) // 0..inf -> 1..inf
                            { m.has_max = false; m.min = cosh(m.min); }
                        else
                            { m.has_max = false; m.min = 1.0; } // Anything between 1..inf
                    }
                }
                else // no min
                {
                    m.has_min = true; m.min = 1.0; // always a lower boundary
                    if(m.has_max) // max, no min
                    {
                        m.min = cosh(m.max); // n..inf
                        m.has_max = false; // No upper boundary
                    }
                    else // no max, no min
                        m.has_max = false; // No upper boundary
                }
                return m;
            }

            case cIf:
            case cAbsIf:
            {
                // No guess which branch is chosen. Produce a spanning min & max.
                MinMaxTree res1 = GetParam(1).CalculateResultBoundaries();
                MinMaxTree res2 = GetParam(2).CalculateResultBoundaries();
                if(!res2.has_min) res1.has_min = false; else if(res1.has_min && res2.min < res1.min) res1.min = res2.min;
                if(!res2.has_max) res1.has_max = false; else if(res1.has_max && res2.max > res1.max) res1.max = res2.max;
                return res1;
            }

            case cMin:
            {
                bool has_unknown_min = false;
                bool has_unknown_max = false;

                MinMaxTree result;
                for(size_t a=0; a<GetParamCount(); ++a)
                {
                    MinMaxTree m = GetParam(a).CalculateResultBoundaries();
                    if(!m.has_min)
                        has_unknown_min = true;
                    else if(!result.has_min || m.min < result.min)
                        result.min = m.min;

                    if(!m.has_max)
                        has_unknown_max = true;
                    else if(!result.has_max || m.max < result.max)
                        result.max = m.max;
                }
                if(has_unknown_min) result.has_min = false;
                if(has_unknown_max) result.has_max = false;
                return result;
            }
            case cMax:
            {
                bool has_unknown_min = false;
                bool has_unknown_max = false;

                MinMaxTree result;
                for(size_t a=0; a<GetParamCount(); ++a)
                {
                    MinMaxTree m = GetParam(a).CalculateResultBoundaries();
                    if(!m.has_min)
                        has_unknown_min = true;
                    else if(!result.has_min || m.min > result.min)
                        result.min = m.min;

                    if(!m.has_max)
                        has_unknown_max = true;
                    else if(!result.has_max || m.max > result.max)
                        result.max = m.max;
                }
                if(has_unknown_min) result.has_min = false;
                if(has_unknown_max) result.has_max = false;
                return result;
            }
            case cAdd:
            {
                /* It's complicated. Follow the logic below. */
                /* Note: This also deals with the following opcodes:
                 *       cNeg, cSub, cRSub
                 */
                MinMaxTree result(0.0, 0.0);
                for(size_t a=0; a<GetParamCount(); ++a)
                {
                    MinMaxTree item = GetParam(a).CalculateResultBoundaries();

                    if(item.has_min) result.min += item.min;
                    else             result.has_min = false;
                    if(item.has_max) result.max += item.max;
                    else             result.has_max = false;

                    if(!result.has_min && !result.has_max) break; // hopeless
                }
                if(result.has_min && result.has_max
                && result.min > result.max) std::swap(result.min, result.max);
                return result;
            }
            case cMul:
            {
                /* It's very complicated. Follow the logic below. */
                struct Value
                {
                    enum ValueType { Finite, MinusInf, PlusInf };
                    ValueType valueType;
                    double value;

                    Value(ValueType t): valueType(t), value(0) {}
                    Value(double v): valueType(Finite), value(v) {}

                    bool isNegative() const
                    {
                        return valueType == MinusInf ||
                            (valueType == Finite && value < 0.0);
                    }

                    void operator*=(const Value& rhs)
                    {
                        if(valueType == Finite && rhs.valueType == Finite)
                            value *= rhs.value;
                        else
                            valueType = (isNegative() != rhs.isNegative() ?
                                         MinusInf : PlusInf);
                    }

                    bool operator<(const Value& rhs) const
                    {
                        return
                            (valueType == MinusInf && rhs.valueType != MinusInf) ||
                            (valueType == Finite &&
                             (rhs.valueType == PlusInf ||
                              (rhs.valueType == Finite && value < rhs.value)));
                    }
                };

                struct MultiplicationRange
                {
                    Value minValue, maxValue;

                    MultiplicationRange():
                        minValue(Value::PlusInf),
                        maxValue(Value::MinusInf) {}

                    void multiply(Value value1, const Value& value2)
                    {
                        value1 *= value2;
                        if(value1 < minValue) minValue = value1;
                        if(maxValue < value1) maxValue = value1;
                    }
                };

                MinMaxTree result(1.0, 1.0);
                for(size_t a=0; a<GetParamCount(); ++a)
                {
                    MinMaxTree item = GetParam(a).CalculateResultBoundaries();
                    if(!item.has_min && !item.has_max) return MinMaxTree(); // hopeless

                    Value minValue0 = result.has_min ? Value(result.min) : Value(Value::MinusInf);
                    Value maxValue0 = result.has_max ? Value(result.max) : Value(Value::PlusInf);
                    Value minValue1 = item.has_min ? Value(item.min) : Value(Value::MinusInf);
                    Value maxValue1 = item.has_max ? Value(item.max) : Value(Value::PlusInf);

                    MultiplicationRange range;
                    range.multiply(minValue0, minValue1);
                    range.multiply(minValue0, maxValue1);
                    range.multiply(maxValue0, minValue1);
                    range.multiply(maxValue0, maxValue1);

                    if(range.minValue.valueType == Value::Finite)
                        result.min = range.minValue.value;
                    else result.has_min = false;

                    if(range.maxValue.valueType == Value::Finite)
                        result.max = range.maxValue.value;
                    else result.has_max = false;

                    if(!result.has_min && !result.has_max) break; // hopeless
                }
                if(result.has_min && result.has_max
                && result.min > result.max) std::swap(result.min, result.max);
                return result;
            }
            case cMod:
            {
                /* TODO: The boundaries of modulo operator could be estimated better. */

                MinMaxTree x = GetParam(0).CalculateResultBoundaries();
                MinMaxTree y = GetParam(1).CalculateResultBoundaries();

                if(y.has_max)
                {
                    if(y.max >= 0.0)
                    {
                        if(!x.has_min || x.min < 0)
                            return MinMaxTree(-y.max, y.max);
                        else
                            return MinMaxTree(0.0, y.max);
                    }
                    else
                    {
                        if(!x.has_max || x.max >= 0)
                            return MinMaxTree(y.max, -y.max);
                        else
                            return MinMaxTree(y.max, NEGATIVE_MAXIMUM);
                    }
                }
                else
                    return MinMaxTree();
            }
            case cPow:
            {
                if(GetParam(1).IsImmed() && GetParam(1).GetImmed() == 0.0)
                {
                    // Note: This makes 0^0 evaluate into 1.
                    return MinMaxTree(1.0, 1.0); // x^0 = 1
                }
                if(GetParam(0).IsImmed() && GetParam(0).GetImmed() == 0.0)
                {
                    // Note: This makes 0^0 evaluate into 0.
                    return MinMaxTree(0.0, 0.0); // 0^x = 0
                }
                if(GetParam(0).IsImmed() && FloatEqual(GetParam(0).GetImmed(), 1.0))
                {
                    return MinMaxTree(1.0, 1.0); // 1^x = 1
                }
                if(GetParam(1).IsImmed()
                && GetParam(1).GetImmed() > 0
                && GetParam(1).IsAlwaysParity(false))
                {
                    // x ^ even_int_const always produces a non-negative value.
                    double exponent = GetParam(1).GetImmed();
                    MinMaxTree tmp = GetParam(0).CalculateResultBoundaries();
                    MinMaxTree result;
                    result.has_min = true;
                    result.min = 0;
                    if(tmp.has_min && tmp.min >= 0)
                        result.min = fp_pow(tmp.min, exponent);
                    else if(tmp.has_max && tmp.max <= 0)
                        result.min = fp_pow(tmp.max, exponent);

                    result.has_max = false;
                    if(tmp.has_min && tmp.has_max)
                    {
                        result.has_max = true;
                        result.max     = std::max(fabs(tmp.min), fabs(tmp.max));
                        result.max     = fp_pow(result.max, exponent);
                    }
                    return result;
                }

                MinMaxTree p0 = GetParam(0).CalculateResultBoundaries();
                MinMaxTree p1 = GetParam(1).CalculateResultBoundaries();
                TriTruthValue p0_positivity =
                    (p0.has_min && p0.min >= 0.0) ? IsAlways
                  : (p0.has_max && p0.max < 0.0 ? IsNever
                    : Unknown);
                TriTruthValue p1_evenness = GetParam(1).GetEvennessInfo();

                /* If param0 IsAlways, the return value is also IsAlways */
                /* If param1 is even, the return value is IsAlways */
                /* If param1 is odd, the return value is same as param0's */
                /* If param0 is negative and param1 is not integer,
                 * the return value is imaginary (assumed Unknown)
                 *
                 * Illustrated in this truth table:
                 *  P=positive, N=negative
                 *  E=even, O=odd, U=not integer
                 *  *=unknown, X=invalid (unknown), x=maybe invalid (unknown)
                 *
                 *   param1: PE PO P* NE NO N* PU NU *
                 * param0:
                 *   PE      P  P  P  P  P  P  P  P  P
                 *   PO      P  P  P  P  P  P  P  P  P
                 *   PU      P  P  P  P  P  P  P  P  P
                 *   P*      P  P  P  P  P  P  P  P  P
                 *   NE      P  N  *  P  N  *  X  X  x
                 *   NO      P  N  *  P  N  *  X  X  x
                 *   NU      P  N  *  P  N  *  X  X  x
                 *   N*      P  N  *  P  N  *  X  X  x
                 *   *       P  *  *  P  *  *  x  x  *
                 *
                 * Note: This also deals with the following opcodes:
                 *       cSqrt  (param0, PU) (x^0.5)
                 *       cRSqrt (param0, NU) (x^-0.5)
                 *       cExp   (PU, param1) (CONSTANT_E^x)
                 */
                TriTruthValue result_positivity = Unknown;
                switch(p0_positivity)
                {
                    case IsAlways:
                        // e.g.   5^x = positive.
                        result_positivity = IsAlways;
                        break;
                    case IsNever:
                    {
                        result_positivity = p1_evenness;
                        break;
                    }
                    default:
                        switch(p1_evenness)
                        {
                            case IsAlways:
                                // e.g. x^( 4) = positive
                                // e.g. x^(-4) = positive
                                result_positivity = IsAlways;
                                break;
                            case IsNever:
                                break;
                            case Unknown:
                            {
                                /* If p1 is const non-integer,
                                 * assume the result is positive
                                 * though it may be NaN instead.
                                 */
                                if(GetParam(1).IsImmed()
                                && GetParam(1).IsAlwaysInteger(false)
                                && GetParam(1).GetImmed() >= 0.0)
                                {
                                    result_positivity = IsAlways;
                                }
                                break;
                            }
                        }
                }
                switch(result_positivity)
                {
                    case IsAlways:
                    {
                        /* The result is always positive.
                         * Figure out whether we know the minimum value. */
                        double min = 0.0;
                        if(p0.has_min && p1.has_min)
                        {
                            min = pow(p0.min, p1.min);
                            if(p0.min < 0.0 && (!p1.has_max || p1.max >= 0.0) && min >= 0.0)
                                min = 0.0;
                        }
                        if(p0.has_min && p0.min >= 0.0 && p0.has_max && p1.has_max)
                        {
                            double max = pow(p0.max, p1.max);
                            if(min > max) std::swap(min, max);
                            return MinMaxTree(min, max);
                        }
                        return MinMaxTree(min, false);
                    }
                    case IsNever:
                    {
                        /* The result is always negative.
                         * TODO: Figure out whether we know the maximum value.
                         */
                        return MinMaxTree(false, NEGATIVE_MAXIMUM);
                    }
                    default:
                    {
                        /* It can be negative or positive.
                         * We know nothing about the boundaries. */
                        break;
                    }
                }
                break;
            }

            /* The following opcodes are processed by GenerateFrom()
             * within fpoptimizer_bytecode_to_codetree.c and thus
             * they will never occur in the calling context for the
             * most of the parsing context. They may however occur
             * at the late phase, so we deal with them.
             */
            case cNeg:
            {
                CodeTree tmp;
                tmp.SetOpcode(cMul);
                tmp.AddParam(CodeTree(-1.0));
                tmp.AddParam(GetParam(0));
                return tmp.CalculateResultBoundaries();
            }
            case cSub: // converted into cMul y -1
            {
                CodeTree tmp, tmp2;
                tmp2.SetOpcode(cNeg);
                tmp2.AddParam(GetParam(1));
                tmp.SetOpcode(cAdd);
                tmp.AddParam(GetParam(0));
                tmp.AddParamMove(tmp2);
                return tmp.CalculateResultBoundaries();
            }
            case cInv: // converted into cPow x -1
            {
                CodeTree tmp;
                tmp.SetOpcode(cPow);
                tmp.AddParam(GetParam(0));
                tmp.AddParam(CodeTree(-1.0));
                return tmp.CalculateResultBoundaries();
            }
            case cDiv: // converted into cPow y -1
            {
                CodeTree tmp, tmp2;
                tmp2.SetOpcode(cInv);
                tmp2.AddParam(GetParam(1));
                tmp.SetOpcode(cMul);
                tmp.AddParam(GetParam(0));
                tmp.AddParamMove(tmp2);
                return tmp.CalculateResultBoundaries();
            }
            case cRad: // converted into cMul x CONSTANT_RD
            {
                CodeTree tmp;
                tmp.SetOpcode(cMul);
                tmp.AddParam(GetParam(0));
                tmp.AddParam(CodeTree(CONSTANT_RD));
                return tmp.CalculateResultBoundaries();
            }
            case cDeg: // converted into cMul x CONSTANT_DR
            {
                CodeTree tmp;
                tmp.SetOpcode(cMul);
                tmp.AddParam(GetParam(0));
                tmp.AddParam(CodeTree(CONSTANT_DR));
                return tmp.CalculateResultBoundaries();
            }
            case cSqr: // converted into cMul x x    or cPow x 2
            {
                CodeTree tmp;
                tmp.SetOpcode(cPow);
                tmp.AddParam(GetParam(0));
                tmp.AddParam(CodeTree(2.0));
                return tmp.CalculateResultBoundaries();
            }
            case cExp: // converted into cPow CONSTANT_E x
            {
                CodeTree tmp;
                tmp.SetOpcode(cPow);
                tmp.AddParam(CodeTree(CONSTANT_E));
                tmp.AddParam(GetParam(0));
                return tmp.CalculateResultBoundaries();
            }
            case cExp2: // converted into cPow 2 x
            {
                CodeTree tmp;
                tmp.SetOpcode(cPow);
                tmp.AddParam(CodeTree(2.0));
                tmp.AddParam(GetParam(0));
                return tmp.CalculateResultBoundaries();
            }
            case cCbrt: // converted into cPow x 0.33333333
            {
                // However, contrary to x^(1/3), this allows
                // negative values for x, and produces those
                // as well.
                MinMaxTree result = GetParam(0).CalculateResultBoundaries();
                if(result.has_min) result.min = fp_cbrt(result.min);
                if(result.has_max) result.max = fp_cbrt(result.max);
                return result;
            }
            case cSqrt: // converted into cPow x 0.5
            {
                MinMaxTree result = GetParam(0).CalculateResultBoundaries();
                if(result.has_min) result.min = result.min < 0 ? 0 : fp_sqrt(result.min);
                if(result.has_max) result.max = result.max < 0 ? 0 : fp_sqrt(result.max);
                return result;
            }
            case cRSqrt: // converted into cPow x -0.5
            {
                CodeTree tmp;
                tmp.SetOpcode(cPow);
                tmp.AddParam(GetParam(0));
                tmp.AddParam(CodeTree(-0.5));
                return tmp.CalculateResultBoundaries();
            }
            case cLog2by: // converted into cMul y CONSTANT_L2I (cLog x)
            {
                CodeTree tmp, tmp2;
                tmp2.SetOpcode(cLog2);
                tmp2.AddParam(GetParam(0));
                tmp.SetOpcode(cMul);
                tmp.AddParamMove(tmp2);
                tmp.AddParam(GetParam(1));
                return tmp.CalculateResultBoundaries();
            }
            case cCot: // converted into 1 / cTan
            {
                CodeTree tmp, tmp2;
                tmp2.SetOpcode(cTan);
                tmp2.AddParam(GetParam(0));
                tmp.SetOpcode(cInv);
                tmp.AddParamMove(tmp2);
                return tmp.CalculateResultBoundaries();
            }
            case cSec: // converted into 1 / cCos
            {
                CodeTree tmp, tmp2;
                tmp2.SetOpcode(cCos);
                tmp2.AddParam(GetParam(0));
                tmp.SetOpcode(cInv);
                tmp.AddParamMove(tmp2);
                return tmp.CalculateResultBoundaries();
            }
            case cCsc: // converted into 1 / cSin
            {
                CodeTree tmp, tmp2;
                tmp2.SetOpcode(cSin);
                tmp2.AddParam(GetParam(0));
                tmp.SetOpcode(cInv);
                tmp.AddParamMove(tmp2);
                return tmp.CalculateResultBoundaries();
            }
            /* The following opcodes are processed by GenerateFrom()
             * within fpoptimizer_bytecode_to_codetree.c and thus
             * they will never occur in the calling context:
             */
                break; /* Should never occur */

            /* Opcodes that do not occur in the tree for other reasons */
            case cRDiv: // version of cDiv
            case cRSub: // version of cSub
            case cDup:
            case cFetch:
            case cPopNMov:
            case cNop:
            case cJump:
            case VarBegin:
                break; /* Should never occur */

            /* Opcodes that are completely unpredictable */
            case cPCall:
                break;
            case cFCall:
                break; // Cannot deduce
            case cEval:
                break; // Cannot deduce
        }
        return MinMaxTree(); /* Cannot deduce */
    }
}

#endif

#line 1 "fpoptimizer/fpoptimizer_transformations.c"
// line removed for fpoptimizer.c: #include "fpoptimizer_bytecodesynth.h"
// line removed for fpoptimizer.c: #include "fpoptimizer_codetree.h"
// line removed for fpoptimizer.c: #include "fpoptimizer_consts.h"

#ifdef FP_SUPPORT_OPTIMIZER

using namespace FUNCTIONPARSERTYPES;
//using namespace FPoptimizer_Grammar;

//#define DEBUG_POWI
//#define DEBUG_SUBSTITUTIONS_CSE

#if defined(__x86_64) || !defined(FP_SUPPORT_CBRT)
# define CBRT_IS_SLOW
#endif

namespace FPoptimizer_ByteCode
{
    extern const unsigned char powi_table[256];
}
namespace
{
    using namespace FPoptimizer_CodeTree;

    typedef
        std::multimap<fphash_t,  std::pair<size_t, CodeTree> >
        TreeCountType;

    void FindTreeCounts(TreeCountType& TreeCounts, const CodeTree& tree)
    {
        TreeCountType::iterator i = TreeCounts.lower_bound(tree.GetHash());
        bool found = false;
        for(; i != TreeCounts.end() && i->first == tree.GetHash(); ++i)
        {
            if(tree.IsIdenticalTo( i->second.second ) )
            {
                i->second.first += 1;
                found = true;
                break;
            }
        }
        if(!found)
        {
            TreeCounts.insert(i, std::make_pair(tree.GetHash(), std::make_pair(size_t(1), tree)));
        }

        for(size_t a=0; a<tree.GetParamCount(); ++a)
            FindTreeCounts(TreeCounts, tree.GetParam(a));
    }

    struct BalanceResultType
    {
        bool BalanceGood;
        bool FoundChild;
    };
    BalanceResultType IfBalanceGood(const CodeTree& root, const CodeTree& child)
    {
        if(root.IsIdenticalTo(child))
        {
            BalanceResultType result = {true,true};
            return result;
        }

        BalanceResultType result = {true,false};

        if(root.GetOpcode() == cIf
        || root.GetOpcode() == cAbsIf)
        {
            BalanceResultType cond    = IfBalanceGood(root.GetParam(0), child);
            BalanceResultType branch1 = IfBalanceGood(root.GetParam(1), child);
            BalanceResultType branch2 = IfBalanceGood(root.GetParam(2), child);

            if(cond.FoundChild || branch1.FoundChild || branch2.FoundChild)
                { result.FoundChild = true; }

            // balance is good if:
            //      branch1.found = branch2.found OR (cond.found AND cond.goodbalance)
            // AND  cond.goodbalance OR (branch1.found AND branch2.found)
            // AND  branch1.goodbalance OR (cond.found AND cond.goodbalance)
            // AND  branch2.goodbalance OR (cond.found AND cond.goodbalance)

            result.BalanceGood =
                (   (branch1.FoundChild == branch2.FoundChild)
                 || (cond.FoundChild && cond.BalanceGood) )
             && (cond.BalanceGood || (branch1.FoundChild && branch2.FoundChild))
             && (branch1.BalanceGood || (cond.FoundChild && cond.BalanceGood))
             && (branch2.BalanceGood || (cond.FoundChild && cond.BalanceGood));
        }
        else
        {
            bool has_bad_balance        = false;
            bool has_good_balance_found = false;

            // Balance is bad if one of the children has bad balance
            // Unless one of the children has good balance & found

            for(size_t b=root.GetParamCount(), a=0; a<b; ++a)
            {
                BalanceResultType tmp = IfBalanceGood(root.GetParam(a), child);
                if(tmp.FoundChild)
                    result.FoundChild = true;

                if(tmp.BalanceGood == false)
                    has_bad_balance = true;
                else if(tmp.FoundChild)
                    has_good_balance_found = true;

                // if the expression is
                //   if(x, sin(x), 0) + sin(x)
                // then sin(x) is a good subexpression
                // even though it occurs in unbalance.
            }
            if(has_bad_balance && !has_good_balance_found)
                result.BalanceGood = false;
        }
        return result;
    }

    bool IsOptimizableUsingPowi(long immed, long penalty = 0)
    {
        FPoptimizer_ByteCode::ByteCodeSynth synth;
        synth.PushVar(0);
        // Ignore the size generated by subtree
        size_t bytecodesize_backup = synth.GetByteCodeSize();
        FPoptimizer_ByteCode::AssembleSequence(immed, FPoptimizer_ByteCode::MulSequence, synth);

        size_t bytecode_grow_amount = synth.GetByteCodeSize() - bytecodesize_backup;

        return bytecode_grow_amount < size_t(MAX_POWI_BYTECODE_LENGTH - penalty);
    }

    void ChangeIntoRootChain(
        CodeTree& tree,
        bool inverted,
        long sqrt_count,
        long cbrt_count)
    {
        while(cbrt_count > 0)
        {
            CodeTree tmp;
            tmp.SetOpcode(cCbrt);
            tmp.AddParamMove(tree);
            tmp.Rehash();
            tree.swap(tmp);
            --cbrt_count;
        }
        while(sqrt_count > 0)
        {
            CodeTree tmp;
            tmp.SetOpcode(cSqrt);
            if(inverted)
            {
                tmp.SetOpcode(cRSqrt);
                inverted = false;
            }
            tmp.AddParamMove(tree);
            tmp.Rehash();
            tree.swap(tmp);
            --sqrt_count;
        }
        if(inverted)
        {
            CodeTree tmp;
            tmp.SetOpcode(cInv);
            tmp.AddParamMove(tree);
            tree.swap(tmp);
        }
    }

    double CalculatePowiFactorCost(long abs_int_exponent)
    {
        static std::map<long, double> cache;
        std::map<long,double>::iterator i = cache.lower_bound(abs_int_exponent);
        if(i != cache.end() && i->first == abs_int_exponent)
            return i->second;
        std::pair<long, double> result(abs_int_exponent, 0.0);
        double& cost = result.second;

        while(abs_int_exponent > 1)
        {
            int factor = 0;
            if(abs_int_exponent < 256)
            {
                factor = FPoptimizer_ByteCode::powi_table[abs_int_exponent];
                if(factor & 128) factor &= 127; else factor = 0;
                if(factor & 64) factor = -(factor&63) - 1;
            }
            if(factor)
            {
                cost += CalculatePowiFactorCost(factor);
                abs_int_exponent /= factor;
                continue;
            }
            if(!(abs_int_exponent & 1))
            {
                abs_int_exponent /= 2;
                cost += 3; // sqr
            }
            else
            {
                cost += 3.5; // dup+mul
                abs_int_exponent -= 1;
            }
        }

        cache.insert(i, result);
        return cost;
    }

    struct PowiResolver
    {
        /* Any exponentiation can be turned into one of these:
         *
         *   x^y  -> sqrt(x)^(y*2)         = x Sqrt       y*2  Pow
         *   x^y  -> cbrt(x)^(y*3)         = x Cbrt       y*3  Pow
         *   x^y  -> rsqrt(x)^(y*-2)       = x RSqrt     y*-2  Pow
         *   x^y  -> x^(y-1/2) * sqrt(x)   = x Sqrt   x y-0.5  Pow Mul
         *   x^y  -> x^(y-1/3) * cbrt(x)   = x Cbrt   x y-0.33 Pow Mul
         *   x^y  -> x^(y+1/2) * rsqrt(x)  = x Sqrt   x y+0.5  Pow Mul
         *   x^y  -> inv(x)^(-y)           = x Inv      -y     Pow
         *
         * These rules can be applied recursively.
         * The goal is to find the optimal chain of operations
         * that results in the least number of sqrt,cbrt operations;
         * an integer value of y, and that the integer is as close
         * to zero as possible.
         */
        static const unsigned MaxSep = 4;

        struct PowiResult
        {
            PowiResult() :
                n_int_sqrt(0),
                n_int_cbrt(0),
                resulting_exponent(0),
                sep_list() { }

            int n_int_sqrt;
            int n_int_cbrt;
            long resulting_exponent;
            int sep_list[MaxSep];
        };

        PowiResult CreatePowiResult(double exponent) const
        {
            static const double RootPowers[(1+4)*(1+3)] =
            {
                // (sqrt^n(x))
                1.0,
                1.0 / (2),
                1.0 / (2*2),
                1.0 / (2*2*2),
                1.0 / (2*2*2*2),
                // cbrt^1(sqrt^n(x))
                1.0 / (3),
                1.0 / (3*2),
                1.0 / (3*2*2),
                1.0 / (3*2*2*2),
                1.0 / (3*2*2*2*2),
                // cbrt^2(sqrt^n(x))
                1.0 / (3*3),
                1.0 / (3*3*2),
                1.0 / (3*3*2*2),
                1.0 / (3*3*2*2*2),
                1.0 / (3*3*2*2*2*2),
                // cbrt^3(sqrt^n(x))
                1.0 / (3*3*3),
                1.0 / (3*3*3*2),
                1.0 / (3*3*3*2*2),
                1.0 / (3*3*3*2*2*2),
                1.0 / (3*3*3*2*2*2*2)
            };

            PowiResult result;

            int best_factor = FindIntegerFactor(exponent);
            if(best_factor == 0)
            {
        #ifdef DEBUG_POWI
            printf("no factor found for %g\n", exponent);
        #endif
                return result; // Unoptimizable
            }

            double best_cost = EvaluateFactorCost(best_factor, 0, 0, 0)
                             + CalculatePowiFactorCost(long(exponent*best_factor));
            int s_count = 0;
            int c_count = 0;
            int mul_count = 0;

        #ifdef DEBUG_POWI
            printf("orig = %g\n", exponent);
            printf("plain factor = %d, cost %g\n", best_factor, best_cost);
        #endif

            for(unsigned n_s=0; n_s<MaxSep; ++n_s)
            {
                int best_selected_sep = 0;
                double best_sep_cost     = best_cost;
                int best_sep_factor   = best_factor;
                for(int s=1; s<5*4; ++s)
                {
#ifdef CBRT_IS_SLOW
                    if(s >= 5) break;
                    // When cbrt is implemented through exp and log,
                    // there is no advantage over exp(log()), so don't support it.
#endif
                    int n_sqrt = s%5;
                    int n_cbrt = s/5;
                    if(n_sqrt + n_cbrt > 4) continue;

                    double changed_exponent = exponent;
                    changed_exponent -= RootPowers[s];

                    int factor = FindIntegerFactor(changed_exponent);
                    if(factor != 0)
                    {
                        double cost = EvaluateFactorCost
                            (factor, s_count + n_sqrt, c_count + n_cbrt, mul_count + 1)
                          + CalculatePowiFactorCost(long(changed_exponent*factor));

        #ifdef DEBUG_POWI
                        printf("%d sqrt %d cbrt factor = %d, cost %g\n",
                            n_sqrt, n_cbrt, factor, cost);
        #endif
                        if(cost < best_sep_cost)
                        {
                            best_selected_sep = s;
                            best_sep_factor   = factor;
                            best_sep_cost     = cost;
                        }
                    }
                }
                if(!best_selected_sep) break;

                result.sep_list[n_s] = best_selected_sep;
                exponent -= RootPowers[best_selected_sep];
                s_count += best_selected_sep % 5;
                c_count += best_selected_sep / 5;
                best_cost   = best_sep_cost;
                best_factor = best_sep_factor;
                mul_count += 1;
            }

            result.resulting_exponent = (long) (exponent * best_factor + 0.5);
            while(best_factor % 2 == 0)
            {
                ++result.n_int_sqrt;
                best_factor /= 2;
            }
            while(best_factor % 3 == 0)
            {
                ++result.n_int_cbrt;
                best_factor /= 3;
            }
            return result;
        }

    private:
        // Find the integer that "value" must be multiplied
        // with to produce an integer...
        // Consisting of factors 2 and 3 only.
        bool MakesInteger(double value, int factor) const
        {
            double v = value * double(factor);
            double diff = fabs(v - (double)(long)(v+0.5));
            //printf("factor %d: v=%.20f, diff=%.20f\n", factor,v, diff);
            return diff < 1e-9;
        }
        int FindIntegerFactor(double value) const
        {
            int factor = (2*2*2*2);
#ifdef CBRT_IS_SLOW
            // When cbrt is implemented through exp and log,
            // there is no advantage over exp(log()), so don't support it.
#else
            factor *= (3*3*3);
#endif
            int result = 0;
            if(MakesInteger(value, factor))
            {
                result = factor;
                while((factor % 2) == 0 && MakesInteger(value, factor/2))
                    result = factor /= 2;
                while((factor % 3) == 0 && MakesInteger(value, factor/3))
                    result = factor /= 3;
            }
#ifdef CBRT_IS_SLOW
            if(result == 0)
            {
                /* Note: Even if we allow one cbrt,
                 *        cbrt(cbrt(x)) still gets turned into
                 *        exp(log(x)*0.111111)
                 *        which gives an error when x < 0...
                 *        should we use a special system here?
                 *        i.e. exp(log(-5)*y)
                 *      =      -exp(log(5)*y)
                 *        except when y is an even integer,
                 *      when  = exp(log(5)*y)
                 * We use a custom fp_pow() function
                 * in order to handle these situations.
                 */
                if(MakesInteger(value, 3)) return 3; // single cbrt opcode
            }
#endif
            return result;
        }

        int EvaluateFactorCost(int factor, int s, int c, int nmuls) const
        {
            const int sqrt_cost = 6;
#ifdef CBRT_IS_SLOW
            const int cbrt_cost = 25;
#else
            const int cbrt_cost = 8;
#endif
            int result = s * sqrt_cost + c * cbrt_cost;
            while(factor % 2 == 0) { factor /= 2; result += sqrt_cost; }
            while(factor % 3 == 0) { factor /= 3; result += cbrt_cost; }
            result += nmuls;
            return result;
        }
    };
}

namespace FPoptimizer_CodeTree
{
    bool CodeTree::RecreateInversionsAndNegations(bool prefer_base2)
    {
        bool changed = false;

        for(size_t a=0; a<GetParamCount(); ++a)
            if(GetParam(a).RecreateInversionsAndNegations(prefer_base2))
                changed = true;

        if(changed)
        {
        exit_changed:
            Mark_Incompletely_Hashed();
            return true;
        }

        switch(GetOpcode()) // Recreate inversions and negations
        {
            case cMul:
            {
                std::vector<CodeTree> div_params;
                CodeTree found_log2, found_log2by;

                if(true)
                {
                    /* This lengthy bit of code
                     * changes log2(x)^3 * 5
                     * to      log2by(x, 5^(1/3)) ^ 3
                     * which is better for runtime
                     * than    log2by(x,1)^3 * 5
                     */
                    bool found_log2_on_exponent = false;
                    double log2_exponent = 0;
                    for(size_t a = GetParamCount(); a-- > 0; )
                    {
                        const CodeTree& powgroup = GetParam(a);
                        if(powgroup.GetOpcode() == cPow
                        && powgroup.GetParam(0).GetOpcode() == cLog2
                        && powgroup.GetParam(1).IsImmed())
                        {
                            // Found log2 on exponent
                            found_log2_on_exponent = true;
                            log2_exponent = powgroup.GetParam(1).GetImmed();
                            break;
                        }
                    }
                    if(found_log2_on_exponent)
                    {
                        double immeds = 1.0;
                        for(size_t a = GetParamCount(); a-- > 0; )
                        {
                            const CodeTree& powgroup = GetParam(a);
                            if(powgroup.IsImmed())
                            {
                                immeds *= powgroup.GetImmed();
                                DelParam(a);
                            }
                        }
                        for(size_t a = GetParamCount(); a-- > 0; )
                        {
                            CodeTree& powgroup = GetParam(a);
                            if(powgroup.GetOpcode() == cPow
                            && powgroup.GetParam(0).GetOpcode() == cLog2
                            && powgroup.GetParam(1).IsImmed())
                            {
                                CodeTree& log2 = powgroup.GetParam(0);
                                log2.CopyOnWrite();
                                log2.SetOpcode(cLog2by);
                                log2.AddParam( CodeTree( fp_pow(immeds, 1.0 / log2_exponent) ) );
                                log2.Rehash();
                                break;
                            }
                        }
                    }
                }

                for(size_t a = GetParamCount(); a-- > 0; )
                {
                    const CodeTree& powgroup = GetParam(a);
                    if(powgroup.GetOpcode() == cPow
                    && powgroup.GetParam(1).IsImmed())
                    {
                        const CodeTree& exp_param = powgroup.GetParam(1);
                        double exponent = exp_param.GetImmed();
                        if(FloatEqual(exponent, -1.0))
                        {
                            CopyOnWrite();
                            div_params.push_back(GetParam(a).GetParam(0));
                            DelParam(a); // delete the pow group
                        }
                        else if(exponent < 0 && IsIntegerConst(exponent))
                        {
                            CodeTree edited_powgroup;
                            edited_powgroup.SetOpcode(cPow);
                            edited_powgroup.AddParam(powgroup.GetParam(0));
                            edited_powgroup.AddParam(CodeTree(-exponent));
                            edited_powgroup.Rehash();
                            div_params.push_back(edited_powgroup);
                            CopyOnWrite();
                            DelParam(a); // delete the pow group
                        }
                    }
                    else if(powgroup.GetOpcode() == cLog2 && !found_log2.IsDefined())
                    {
                        found_log2 = powgroup.GetParam(0);
                        CopyOnWrite();
                        DelParam(a);
                    }
                    else if(powgroup.GetOpcode() == cLog2by && !found_log2by.IsDefined())
                    {
                        found_log2by = powgroup;
                        CopyOnWrite();
                        DelParam(a);
                    }
                }
                if(!div_params.empty())
                {
                    changed = true;

                    CodeTree divgroup;
                    divgroup.SetOpcode(cMul);
                    divgroup.SetParamsMove(div_params);
                    divgroup.Rehash(); // will reduce to div_params[0] if only one item
                    CodeTree mulgroup;
                    mulgroup.SetOpcode(cMul);
                    mulgroup.SetParamsMove(GetParams());
                    mulgroup.Rehash(); // will reduce to 1.0 if none remained in this cMul
                    if(mulgroup.IsImmed() && FloatEqual(mulgroup.GetImmed(), 1.0))
                    {
                        SetOpcode(cInv);
                        AddParamMove(divgroup);
                    }
                    /*else if(mulgroup.IsImmed() && FloatEqual(mulgroup.GetImmed(), -1.0))
                    {
                        CodeTree invgroup;
                        invgroup.SetOpcode(cInv);
                        invgroup.AddParamMove(divgroup);
                        invgroup.Rehash();
                        SetOpcode(cNeg);
                        AddParamMove(invgroup);
                    }*/
                    else
                    {
                        if(mulgroup.GetDepth() >= divgroup.GetDepth())
                        {
                            SetOpcode(cDiv);
                            AddParamMove(mulgroup);
                            AddParamMove(divgroup);
                        }
                        else
                        {
                            SetOpcode(cRDiv);
                            AddParamMove(divgroup);
                            AddParamMove(mulgroup);
                        }
                    }
                }
                if(found_log2.IsDefined())
                {
                    CodeTree mulgroup;
                    mulgroup.SetOpcode(GetOpcode());
                    mulgroup.SetParamsMove(GetParams());
                    mulgroup.Rehash();
                    while(mulgroup.RecreateInversionsAndNegations(prefer_base2))
                        mulgroup.FixIncompleteHashes();
                    SetOpcode(cLog2by);
                    AddParamMove(found_log2);
                    AddParamMove(mulgroup);
                    changed = true;
                }
                if(found_log2by.IsDefined())
                {
                    CodeTree mulgroup;
                    mulgroup.SetOpcode(cMul);
                    mulgroup.AddParamMove(found_log2by.GetParam(1));
                    mulgroup.AddParamsMove(GetParams());
                    mulgroup.Rehash();
                    while(mulgroup.RecreateInversionsAndNegations(prefer_base2))
                        mulgroup.FixIncompleteHashes();
                    DelParams();
                    SetOpcode(cLog2by);
                    AddParamMove(found_log2by.GetParam(0));
                    AddParamMove(mulgroup);
                    changed = true;
                }
                break;
            }
            case cAdd:
            {
                std::vector<CodeTree> sub_params;

                for(size_t a = GetParamCount(); a-- > 0; )
                    if(GetParam(a).GetOpcode() == cMul)
                    {
                        bool is_signed = false; // if the mul group has a -1 constant...

                        CodeTree& mulgroup = GetParam(a);

                        for(size_t b=mulgroup.GetParamCount(); b-- > 0; )
                        {
                            if(mulgroup.GetParam(b).IsImmed())
                            {
                                double factor = mulgroup.GetParam(b).GetImmed();
                                if(FloatEqual(factor, -1.0))
                                {
                                    mulgroup.CopyOnWrite();
                                    mulgroup.DelParam(b);
                                    is_signed = !is_signed;
                                }
                                else if(FloatEqual(factor, -2.0))
                                {
                                    mulgroup.CopyOnWrite();
                                    mulgroup.DelParam(b);
                                    mulgroup.AddParam( CodeTree(2.0) );
                                    is_signed = !is_signed;
                                }
                            }
                        }
                        if(is_signed)
                        {
                            mulgroup.Rehash();
                            sub_params.push_back(mulgroup);
                            CopyOnWrite();
                            DelParam(a);
                        }
                    }
                    else if(GetParam(a).GetOpcode() == cDiv)
                    {
                        bool is_signed = false;
                        CodeTree& divgroup = GetParam(a);
                        if(divgroup.GetParam(0).IsImmed())
                        {
                            if(FloatEqual(divgroup.GetParam(0).GetImmed(), -1.0))
                            {
                                divgroup.CopyOnWrite();
                                divgroup.DelParam(0);
                                divgroup.SetOpcode(cInv);
                                is_signed = !is_signed;
                            }
                        }
                        if(is_signed)
                        {
                            divgroup.Rehash();
                            sub_params.push_back(divgroup);
                            CopyOnWrite();
                            DelParam(a);
                        }
                    }
                    else if(GetParam(a).GetOpcode() == cRDiv)
                    {
                        bool is_signed = false;
                        CodeTree& divgroup = GetParam(a);
                        if(divgroup.GetParam(1).IsImmed())
                        {
                            if(FloatEqual(divgroup.GetParam(1).GetImmed(), -1.0))
                            {
                                divgroup.CopyOnWrite();
                                divgroup.DelParam(1);
                                divgroup.SetOpcode(cInv);
                                is_signed = !is_signed;
                            }
                        }
                        if(is_signed)
                        {
                            divgroup.Rehash();
                            sub_params.push_back(divgroup);
                            CopyOnWrite();
                            DelParam(a);
                        }
                    }
                if(!sub_params.empty())
                {
                    CodeTree subgroup;
                    subgroup.SetOpcode(cAdd);
                    subgroup.SetParamsMove(sub_params);
                    subgroup.Rehash(); // will reduce to sub_params[0] if only one item
                    CodeTree addgroup;
                    addgroup.SetOpcode(cAdd);
                    addgroup.SetParamsMove(GetParams());
                    addgroup.Rehash(); // will reduce to 0.0 if none remained in this cAdd
                    if(addgroup.IsImmed() && FloatEqual(addgroup.GetImmed(), 0.0))
                    {
                        SetOpcode(cNeg);
                        AddParamMove(subgroup);
                    }
                    else
                    {
                        if(addgroup.GetDepth() == 1)
                        {
                            /* 5 - (x+y+z) is best expressed as rsub(x+y+z, 5);
                             * this has lowest stack usage.
                             * This is identified by addgroup having just one member.
                             */
                            SetOpcode(cRSub);
                            AddParamMove(subgroup);
                            AddParamMove(addgroup);
                        }
                        else if(subgroup.GetOpcode() == cAdd)
                        {
                            /* a+b-(x+y+z) is expressed as a+b-x-y-z.
                             * Making a long chain of cSubs is okay, because the
                             * cost of cSub is the same as the cost of cAdd.
                             * Thus we get the lowest stack usage.
                             * This approach cannot be used for cDiv.
                             */
                            SetOpcode(cSub);
                            AddParamMove(addgroup);
                            AddParamMove(subgroup.GetParam(0));
                            for(size_t a=1; a<subgroup.GetParamCount(); ++a)
                            {
                                CodeTree innersub;
                                innersub.SetOpcode(cSub);
                                innersub.SetParamsMove(GetParams());
                                innersub.Rehash(false);
                                //DelParams();
                                AddParamMove(innersub);
                                AddParamMove(subgroup.GetParam(a));
                            }
                        }
                        else
                        {
                            SetOpcode(cSub);
                            AddParamMove(addgroup);
                            AddParamMove(subgroup);
                        }
                    }
                }
                break;
            }
            case cPow:
            {
                const CodeTree& p0 = GetParam(0);
                const CodeTree& p1 = GetParam(1);
                if(p1.IsImmed())
                {
                    if(p1.GetImmed() != 0.0 && !p1.IsLongIntegerImmed())
                    {
                        PowiResolver::PowiResult
                            r = PowiResolver().CreatePowiResult(fabs(p1.GetImmed()));

                        if(r.resulting_exponent != 0)
                        {
                            bool signed_chain = false;

                            if(p1.GetImmed() < 0
                            && r.sep_list[0] == 0
                            && r.n_int_sqrt > 0)
                            {
                                // If one of the internal sqrts can be changed into rsqrt
                                signed_chain = true;
                            }

                        #ifdef DEBUG_POWI
                            printf("Will resolve powi %g as powi(chain(%d,%d),%ld)",
                                fabs(p1.GetImmed()),
                                r.n_int_sqrt,
                                r.n_int_cbrt,
                                r.resulting_exponent);
                            for(unsigned n=0; n<PowiResolver::MaxSep; ++n)
                            {
                                if(r.sep_list[n] == 0) break;
                                int n_sqrt = r.sep_list[n] % 5;
                                int n_cbrt = r.sep_list[n] / 5;
                                printf("*chain(%d,%d)", n_sqrt,n_cbrt);
                            }
                            printf("\n");
                        #endif

                            CodeTree source_tree = GetParam(0);

                            CodeTree pow_item = source_tree;
                            pow_item.CopyOnWrite();
                            ChangeIntoRootChain(pow_item,
                                signed_chain,
                                r.n_int_sqrt,
                                r.n_int_cbrt);
                            pow_item.Rehash();

                            CodeTree pow;
                            if(r.resulting_exponent != 1)
                            {
                                pow.SetOpcode(cPow);
                                pow.AddParamMove(pow_item);
                                pow.AddParam(CodeTree( double(r.resulting_exponent) ));
                            }
                            else
                                pow.swap(pow_item);

                            CodeTree mul;
                            mul.SetOpcode(cMul);
                            mul.AddParamMove(pow);

                            for(unsigned n=0; n<PowiResolver::MaxSep; ++n)
                            {
                                if(r.sep_list[n] == 0) break;
                                int n_sqrt = r.sep_list[n] % 5;
                                int n_cbrt = r.sep_list[n] / 5;

                                CodeTree mul_item = source_tree;
                                mul_item.CopyOnWrite();
                                ChangeIntoRootChain(mul_item, false, n_sqrt, n_cbrt);
                                mul_item.Rehash();
                                mul.AddParamMove(mul_item);
                            }

                            if(p1.GetImmed() < 0 && !signed_chain)
                            {
                                mul.Rehash();
                                SetOpcode(cInv);
                                SetParamMove(0, mul);
                                DelParam(1);
                            }
                            else
                            {
                                SetOpcode(cMul);
                                SetParamsMove(mul.GetParams());
                            }
                        #ifdef DEBUG_POWI
                            DumpTreeWithIndent(*this);
                        #endif
                            changed = true;
                            break;
                        }
                    }
                }
                if(GetOpcode() == cPow
                && (!p1.IsLongIntegerImmed()
                 || !IsOptimizableUsingPowi(p1.GetLongIntegerImmed())))
                {
                    if(p0.IsImmed() && p0.GetImmed() > 0.0)
                    {
                        // Convert into cExp or Exp2.
                        //    x^y = exp(log(x) * y) =
                        //    Can only be done when x is positive, though.
                        if(prefer_base2)
                        {
                            double mulvalue = fp_log2( p0.GetImmed() );
                            if(mulvalue == 1.0)
                            {
                                // exp2(1)^x becomes exp2(x)
                                DelParam(0);
                            }
                            else
                            {
                                // exp2(4)^x becomes exp2(4*x)
                                CodeTree exponent;
                                exponent.SetOpcode(cMul);
                                exponent.AddParam( CodeTree( mulvalue ) );
                                exponent.AddParam(p1);
                                exponent.Rehash();
                                SetParamMove(0, exponent);
                                DelParam(1);
                            }
                            SetOpcode(cExp2);
                            changed = true;
                        }
                        else
                        {
                            double mulvalue = std::log( p0.GetImmed() );
                            if(mulvalue == 1.0)
                            {
                                // exp(1)^x becomes exp(x)
                                DelParam(0);
                            }
                            else
                            {
                                // exp(4)^x becomes exp(4*x)
                                CodeTree exponent;
                                exponent.SetOpcode(cMul);
                                exponent.AddParam( CodeTree( mulvalue ) );
                                exponent.AddParam(p1);
                                exponent.Rehash();
                                SetParamMove(0, exponent);
                                DelParam(1);
                            }
                            SetOpcode(cExp);
                            changed = true;
                        }
                    }
                    else if(p0.IsAlwaysSigned(true))
                    {
                        if(prefer_base2)
                        {
                            CodeTree log;
                            log.SetOpcode(cLog2);
                            log.AddParam(p0);
                            log.Rehash();
                            CodeTree exponent;
                            exponent.SetOpcode(cMul);
                            exponent.AddParam(p1);
                            exponent.AddParamMove(log);
                            exponent.Rehash();
                            SetOpcode(cExp2);
                            SetParamMove(0, exponent);
                            DelParam(1);
                            changed = true;
                        }
                        else
                        {
                            CodeTree log;
                            log.SetOpcode(cLog);
                            log.AddParam(p0);
                            log.Rehash();
                            CodeTree exponent;
                            exponent.SetOpcode(cMul);
                            exponent.AddParam(p1);
                            exponent.AddParamMove(log);
                            exponent.Rehash();
                            SetOpcode(cExp);
                            SetParamMove(0, exponent);
                            DelParam(1);
                            changed = true;
                        }
                    }
                }
                break;
            }

            default: break;
        }

        if(changed)
            goto exit_changed;

        return changed;
    }

    bool ContainsOtherCandidates(
        const CodeTree& within,
        const CodeTree& tree,
        const FPoptimizer_ByteCode::ByteCodeSynth& synth,
        const TreeCountType& TreeCounts)
    {
        for(size_t b=tree.GetParamCount(), a=0; a<b; ++a)
        {
            const CodeTree& leaf = tree.GetParam(a);

            TreeCountType::iterator synth_it;
            for(TreeCountType::const_iterator
                i = TreeCounts.begin();
                i != TreeCounts.end();
                ++i)
            {
                if(i->first != leaf.GetHash())
                    continue;

                size_t          score     = i->second.first;
                const CodeTree& candidate = i->second.second;

                // It must not yet have been synthesized
                if(synth.Find(candidate))
                    continue;

                // And it must not be a simple expression
                // Because cImmed, VarBegin are faster than cFetch
                if(leaf.GetDepth() <= 1)
                    continue;

                // It must always occur at least twice
                if(score < 2)
                    continue;

                // And it must either appear on both sides
                // of a cIf, or neither
                if(IfBalanceGood(within, leaf).BalanceGood == false)
                    continue;

                return true;
            }
            if(ContainsOtherCandidates(within, leaf, synth, TreeCounts))
                return true;
        }
        return false;
    }

    size_t CodeTree::SynthCommonSubExpressions(
        FPoptimizer_ByteCode::ByteCodeSynth& synth) const
    {
        size_t stacktop_before = synth.GetStackTop();

        /* Find common subtrees */
        TreeCountType TreeCounts;
        FindTreeCounts(TreeCounts, *this);

    #ifdef DEBUG_SUBSTITUTIONS_CSE
        DumpHashes(*this);
    #endif

        /* Synthesize some of the most common ones */
        for(;;)
        {
            size_t best_score = 0;
            TreeCountType::iterator synth_it;
            for(TreeCountType::iterator
                j,i = TreeCounts.begin();
                i != TreeCounts.end();
                i=j)
            {
                j=i; ++j;

                size_t         score = i->second.first;
                const CodeTree& tree = i->second.second;

    #ifdef DEBUG_SUBSTITUTIONS_CSE
                std::cout << "Score " << score << ":\n";
                DumpTreeWithIndent(tree);
    #endif

                // It must not yet have been synthesized
                if(synth.Find(tree))
                {
                    TreeCounts.erase(i);
                    continue;
                }

                // And it must not be a simple expression
                // Because cImmed, VarBegin are faster than cFetch
                if(tree.GetDepth() <= 1)
                {
                    TreeCounts.erase(i);
                    continue;
                }

                // It must always occur at least twice
                if(score < 2)
                {
                    TreeCounts.erase(i);
                    continue;
                }

                // And it must either appear on both sides
                // of a cIf, or neither
                if(IfBalanceGood(*this, tree).BalanceGood == false)
                {
                    TreeCounts.erase(i);
                    continue;
                }

                // It must not contain other candidates
                if(ContainsOtherCandidates(*this, tree, synth, TreeCounts))
                {
                    // Don't erase it; it may be a proper candidate later
                    continue;
                }

                // Is a candidate.
                score *= tree.GetDepth();
                if(score > best_score)
                    { best_score = score; synth_it = i; }
            }

            if(best_score <= 0) break; // Didn't find anything.

            const CodeTree& tree = synth_it->second.second;
    #ifdef DEBUG_SUBSTITUTIONS_CSE
            std::cout << "Found Common Subexpression:"; DumpTree(tree); std::cout << "\n";
    #endif
            /* Synthesize the selected tree */
            tree.SynthesizeByteCode(synth, false);
            TreeCounts.erase(synth_it);
    #ifdef DEBUG_SUBSTITUTIONS_CSE
            std::cout << "Done with Common Subexpression:"; DumpTree(tree); std::cout << "\n";
    #endif
        }

        return synth.GetStackTop() - stacktop_before;
    }
}

#endif

#line 1 "fpoptimizer/fpoptimizer_main.c"
#include "fpconfig.h"
#include "fparser.h"
#include "fptypes.h"

// line removed for fpoptimizer.c: #include "fpoptimizer_codetree.h"
// line removed for fpoptimizer.c: #include "fpoptimizer_optimize.h"

using namespace FUNCTIONPARSERTYPES;

#ifdef FP_SUPPORT_OPTIMIZER
using namespace FPoptimizer_CodeTree;

template<>
void FunctionParserBase<double>::Optimize()
{
    typedef double Value_t;

    CopyOnWrite();

    //PrintByteCode(std::cout);

    CodeTree tree;
    tree.GenerateFrom(data->ByteCode, data->Immed, *data);

    FPoptimizer_Optimize::ApplyGrammars(tree);

    std::vector<unsigned> byteCode;
    std::vector<Value_t> immed;
    size_t stacktop_max = 0;
    tree.SynthesizeByteCode(byteCode, immed, stacktop_max);

    /*std::cout << std::flush;
    std::cerr << std::flush;
    fprintf(stderr, "Estimated stacktop %u\n", (unsigned)stacktop_max);
    fflush(stderr);*/

    if(data->StackSize != stacktop_max)
    {
        data->StackSize = stacktop_max; // note: gcc warning is meaningful
#ifndef FP_USE_THREAD_SAFE_EVAL
        data->Stack.resize(stacktop_max);
#endif
    }

    data->ByteCode.swap(byteCode);
    data->Immed.swap(immed);

    //PrintByteCode(std::cout);
}

template<>
void FunctionParserBase<float>::Optimize()
{}

template<>
void FunctionParserBase<long double>::Optimize()
{}

template<>
void FunctionParserBase<long>::Optimize()
{}

#ifdef FP_SUPPORT_MPFR_FLOAT_TYPE
template<>
void FunctionParserBase<MpfrFloat>::Optimize()
{}
#endif

#ifdef FP_SUPPORT_GMP_INT_TYPE
template<>
void FunctionParserBase<GmpInt>::Optimize()
{}
#endif


FUNCTIONPARSER_INSTANTIATE_TYPES

#endif

#line 1 "fpoptimizer/fpoptimizer_footer.txt"

#endif

