
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_Parser_h
#define MEDYAN_Parser_h

#include <algorithm>
#include <cctype> // isspace
#include <filesystem>
#include <fstream>
#include <iterator>
#include <list>
#include <stdexcept>
#include <string>
#include <sstream>
#include <variant>
#include <vector>

#include "common.h"
#include "SysParams.h"
#include "utility.h"

// A tokenizer for input files
//
// This extends the input config to s-expressions, which allows more
// extensibility, but it does not include the lisp syntax
//
// The input file will be treated as a list of s-expressions. For backwards
// compatibility, several special specifications exist:
//   - All symbols are parsed as strings.
//   - '#' and ';' both start a line comment.
//   - If a top level (not within any parentheses) token is not a parenthesis,
//     an implicit pair of parentheses will be added before the token and
//     before the next closest top level line break or end of input.
//   - Double quotation mark is parsed as string, but nested quotations,
//     escapings are not allowed. Line comment marker in a string stay as it
//     is.
//   - Evaluation is currently not supported.
//   - Most syntactic sugar is currently not supported.

namespace medyan {

struct ConfigFileToken {
    enum class Type {
        string,
        comment,
        parenthesisLeft,
        parenthesisRight,
        lineBreak,         // needed for implicit parenthesis
        unknown
    };

    Type        type = Type::string;
    std::string content;

    static std::string defaultContent(Type type) {
        switch(type) {
        case Type::string:           return "";
        case Type::comment:          return "";
        case Type::parenthesisLeft:  return "(";
        case Type::parenthesisRight: return ")";
        case Type::lineBreak:        return "\n";
        case Type::unknown:          return "";
        default:                     return "";
        }
    }

    static auto makeDefault(Type type) {
        return ConfigFileToken { type, defaultContent(type) };
    }
};

struct SExpr {
    using StringType = std::string;
    using ListType   = std::vector< SExpr >;

    std::variant<
        StringType,
        ListType
    > data;
};


// Tokenize and the input file.
// Comments are removed.
inline std::list< ConfigFileToken > tokenizeConfigFile(const std::string& str) {
    const auto len = str.length();
    const auto isSpecialChar = [](char x) {
        return std::isspace(x) ||
            x == '#' || x == ';' ||
            x == '(' || x == ')' ||
            x == '"';
    };

    std::list< ConfigFileToken > tokens;

    for(int i = 0; i < len; ++i) {
        const char a = str[i];
        if(a == '\n') {
            tokens.push_back(ConfigFileToken::makeDefault(ConfigFileToken::Type::lineBreak));
        }
        else if(std::isspace(a)) {
            // Do nothing
        }
        else if(a == '#' || a == ';') {
            int j = i + 1;
            while(j < len && str[j] != '\n') ++j;

            // Now j points to either the end of string, or the next line break
            tokens.push_back({ConfigFileToken::Type::comment, str.substr(i, j-i)});
            i = j - 1;
        }
        else if(a == '(') {
            tokens.push_back(ConfigFileToken::makeDefault(ConfigFileToken::Type::parenthesisLeft));
        }
        else if(a == ')') {
            tokens.push_back(ConfigFileToken::makeDefault(ConfigFileToken::Type::parenthesisRight));
        }
        else if(a == '"') {
            int j = i + 1;
            while(j < len && str[j] != '"') ++j;

            // Now j points to either end of string, or the next double quotation
            if(j < len) {
                tokens.push_back({ConfigFileToken::Type::string, str.substr(i+1, j-i-1)});
                i = j;
            } else {
                LOG(ERROR) << "Quotation marks do not match";
                throw std::runtime_error("Quotation marks do not match.");
            }
        }
        else {
            int j = i + 1;
            while(j < len && !isSpecialChar(str[j])) ++j;

            // Now j points to either the end of string, or the next special char
            tokens.push_back({ConfigFileToken::Type::string, str.substr(i, j-i)});
            i = j - 1;
        }
    }

    return tokens;
}

inline void outputConfigTokens(std::ostream& os, const std::list< ConfigFileToken >& tokens) {
    auto lastType = ConfigFileToken::Type::unknown;
    for(const auto& tok : tokens) {
        if(tok.type == ConfigFileToken::Type::string) {

            auto content = tok.content;
            // Strip double quotations
            content.erase(
                std::remove(content.begin(), content.end(), '"'),
                content.end()
            );

            // Check whether spaces exist
            const bool hasSpace = (
                std::find_if(content.begin(), content.end(), [](char x) { return std::isspace(x); })
                    != content.end());

            if(
                lastType == ConfigFileToken::Type::string ||
                lastType == ConfigFileToken::Type::parenthesisRight
            ) {
                os << ' ';
            }

            if(hasSpace) {
                os << '"' << content << '"';
            } else {
                os << content;
            }
        }
        else if(tok.type == ConfigFileToken::Type::comment) {
            if(
                lastType != ConfigFileToken::Type::lineBreak &&
                lastType != ConfigFileToken::Type::unknown
            ) {
                os << ' ';
            }
            os << tok.content;
        }
        else if(tok.type == ConfigFileToken::Type::parenthesisLeft) {
            if(
                lastType == ConfigFileToken::Type::string ||
                lastType == ConfigFileToken::Type::parenthesisRight
            ) {
                os << ' ';
            }

            os << tok.content;
        }
        else {
            os << tok.content;
        }

        lastType = tok.type;
    }
}

inline SExpr lexConfigTokensList(
    const std::list< ConfigFileToken >&           tokens,
    std::list< ConfigFileToken >::const_iterator& tokenIter,
    const int                                     depth,
    const bool                                    implicitParentheses
) {
    using namespace std;
    using VS = vector< SExpr >;

    SExpr se;
    se.data = VS{};

    while(tokenIter != tokens.cend()) {
        if(tokenIter->type == ConfigFileToken::Type::string) {
            // Add to the current list
            get<VS>(se.data).push_back(
                SExpr { tokenIter->content }
            );
            ++tokenIter;
        }
        else if(tokenIter->type == ConfigFileToken::Type::parenthesisLeft) {
            get<VS>(se.data).push_back(
                lexConfigTokensList(tokens, ++tokenIter, depth + 1, false)
            );
        }
        else if(tokenIter->type == ConfigFileToken::Type::parenthesisRight) {
            if(implicitParentheses) {
                LOG(ERROR) << "Unexpected ')' in implicit parentheses.";
                throw runtime_error("Unmatched parentheses");
            } else {
                // End current list
                ++tokenIter;
                return se;
            }
        }
        else if(tokenIter->type == ConfigFileToken::Type::lineBreak) {
            if(implicitParentheses) {
                // End current list
                ++tokenIter;
                return se;
            }
            else {
                ++tokenIter;
            }
        }
        else {
            ++tokenIter;
        }
    }

    if(implicitParentheses) {
        // Reaching the end
        return se;
    }
    else {
        LOG(ERROR) << "')' is expected, but not found.";
        throw runtime_error("Unmatched parentheses");
    }
}
inline SExpr lexConfigTokens(const list< ConfigFileToken >& tokens) {
    using namespace std;
    using VS = vector< SExpr >;

    SExpr se;
    se.data = VS{};

    auto tokenIter = tokens.cbegin();
    for(auto tokenIter = tokens.cbegin(); tokenIter != tokens.cend(); ) {

        if (tokenIter->type == ConfigFileToken::Type::string) {
            // Implicit parentheses assumed
            get<VS>(se.data).push_back(
                lexConfigTokensList(tokens, tokenIter, 1, true)
            );
        }
        else if (tokenIter->type == ConfigFileToken::Type::parenthesisLeft) {
            get<VS>(se.data).push_back(
                lexConfigTokensList(tokens, ++tokenIter, 1, false)
            );
        }
        else if (tokenIter->type == ConfigFileToken::Type::parenthesisRight) {
            LOG(ERROR) << "Unexpected ')'.";
            throw runtime_error("Unmatched parentheses");
        }
        else {
            ++tokenIter;
        }

    }

    return se;
}

} // namespace medyan

/// A general parser
/*!
 *  A parser object, when initialized, opens an input file. Upon destruction, it 
 *  closes the file.
 */
class Parser {
protected:
    fstream _inputFile; ///< input file being used
    
public:
    Parser(string inputFileName) {
        _inputFile.open(inputFileName);
        if(!_inputFile.is_open()) {
            cout << "There was an error parsing file " << inputFileName
                 << ". Exiting." << endl;
            exit(EXIT_FAILURE);
        }
        cout << "Loading file " << inputFileName << endl;
    }
    ~Parser() {_inputFile.close();}
};

/// To parse a system input file, initialized by the Controller.
struct SystemParser {
    
    //@{
    /// Parameter parser. Reads input directly into system parameters
    /// @note - does not check for correctness and consistency here.
    static MechParams    readMechParams(std::istream&);
    static ChemParams    readChemParams(std::istream&, const GeoParams&);
    static GeoParams     readGeoParams(std::istream&);
    static BoundParams   readBoundParams(std::istream&, const GeoParams&);
    static DyRateParams  readDyRateParams(std::istream&);
    static SpecialParams readSpecialParams(std::istream&);
    //@}
    
    //@{
    /// Algorithm parser
    static MechParams::MechanicsAlgorithm readMechanicsAlgorithm(std::istream&);
    static ChemParams::ChemistryAlgorithm readChemistryAlgorithm(std::istream&);
    //@}
    
    //@{
    /// Type parser
    static MechParams::MechanicsFFType     readMechanicsFFType(std::istream&);
    static DyRateParams::DynamicRateType   readDynamicRateType(std::istream&);
    static BoundParams::BoundaryType       readBoundaryType(std::istream&);
    static SpecialParams::SpecialSetupType readSpecialSetupType(std::istream&);
    //@}
    
    /// Read Filament information
    static FilamentSetup readFilamentSetup(std::istream&);
    
    /// Read Bubble information
    static BubbleSetup readBubbleSetup(std::istream&);
    
    /// Chemistry information
    static ChemParams::ChemistrySetup readChemistrySetup(std::istream&);
};

/// Used to parse initial Filament information, initialized by the Controller.
struct FilamentParser {
  
    /// Reads filament input file. Returns a vector of tuples containing
    /// filament type and positions (start and end points).
    /// @note - Does not check for coordinate correctness.
    static FilamentData readFilaments(std::istream&);
};

/// Used to parse initial Bubble information, initialized by the Controller.
struct BubbleParser {    
    /// Reads bubble input file. Returns a vector of tuples containing
    /// bubble type and position.
    /// @note - Does not check for coordinate correctness.
    static BubbleData readBubbles(std::istream&);
};


/// Used to parse all chemical information, initialized by the Controller.
struct ChemistryParser {
    /// Reads chemical reactions and species from input file. Returns a
    /// ChemistryData struct containing this data
    /// @note - this does not check for consistency and correctness, the only
    ///         sanity check here is that there are no duplicate species names.
    static ChemistryData readChemistryInput(std::istream&, const ChemParams&);
};


/// Used to parse pin positions if needed upon restart
class PinRestartParser: public Parser {
    
public:
    PinRestartParser(string inputFileName) : Parser(inputFileName) {}
    ~PinRestartParser() {}
    
    /// Reads pin positions from file, and sets filaments
    void resetPins();
};



#endif
