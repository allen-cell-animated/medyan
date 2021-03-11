
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
#include <charconv>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iterator>
#include <list>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <sstream>
#include <type_traits>
#include <variant>
#include <vector>

#include "common.h"
#include "SysParams.h"
#include "utility.h"
#include "Util/Math/Vec.hpp"
#include "Util/SExpr.hpp"

namespace medyan {

// Tokenizer and lexer for input files
//
// This extends the input config to s-expressions, which allows more extensibility, but it does not include the lisp syntax.
//
// The input file will be treated as a quoted list of s-expressions. For backwards compatibility and simplicity, several special specifications exist:
//   - All symbols are parsed as strings.
//   - '#' and ';' both start a line comment.
//   - If a top level (not within any parentheses) token is not a parenthesis, an implicit pair of parentheses will be added before the token and before the next closest top level line break or end of input.
//   - Double quotation mark is parsed as string, but nested quotations, escapings are not allowed. Line comment marker in a quoted string has no effect.
//   - cons that is not a list is currently not supported.
//   - Evaluation is currently not supported.
//   - Most syntactic sugar is currently not supported.

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

    static auto makeString(std::string s) {
        return ConfigFileToken { Type::string, std::move(s) };
    }
};


inline std::vector< std::string > getStringVector(const SExpr::ListType& sel) {
    // Precondition:
    //   - Every element in sel is of string type
    using namespace std;

    vector< string > res;
    res.reserve(sel.size());
    for(const auto& eachSE : sel) {
        res.push_back(get< SExpr::StringType >(eachSE.data));
    }
    return res;
}


inline bool isTokenSpecialChar(char x) {
    return std::isspace(x) ||
        x == '#' || x == ';' ||
        x == '(' || x == ')' ||
        x == '"';
}

// Tokenize and the input file.
inline std::list< ConfigFileToken > tokenizeConfigFile(const std::string& str) {
    const auto len = str.length();

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
                tokens.push_back(ConfigFileToken::makeString(str.substr(i+1, j-i-1)));
                i = j;
            } else {
                LOG(ERROR) << "Quotation marks do not match";
                throw std::runtime_error("Quotation marks do not match.");
            }
        }
        else {
            int j = i + 1;
            while(j < len && !isTokenSpecialChar(str[j])) ++j;

            // Now j points to either the end of string, or the next special char
            tokens.push_back(ConfigFileToken::makeString(str.substr(i, j-i)));
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

            // Check special characters exist
            const bool hasSpace = (
                std::find_if(content.begin(), content.end(), isTokenSpecialChar)
                    != content.end());

            if(
                lastType == ConfigFileToken::Type::string ||
                lastType == ConfigFileToken::Type::parenthesisRight
            ) {
                os << ' ';
            }

            if(hasSpace || content.empty()) {
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

// Build tokens from s-expression
// No comment or line break will be created
inline list< ConfigFileToken > buildConfigTokens(const SExpr& se) {
    using namespace std;
    using TokenList = list< ConfigFileToken >;

    TokenList tokens;

    struct TokenBuildVisitor {
        TokenList& theList;
        void operator()(const SExpr::StringType& str) const {
            theList.push_back(ConfigFileToken::makeString( str ));
        }
        void operator()(const SExpr::ListType& seList) const {
            theList.push_back(ConfigFileToken::makeDefault(ConfigFileToken::Type::parenthesisLeft));
            for(const auto& eachSE : seList) {
                visit(TokenBuildVisitor{theList}, eachSE.data);
            }
            theList.push_back(ConfigFileToken::makeDefault(ConfigFileToken::Type::parenthesisRight));
        }
    };

    visit(TokenBuildVisitor{ tokens }, se.data);
    return tokens;
}


// Key-value parser
//
// Treats an input s-expression as a list of key-value pairs and parse by
// matching keywords
//
// Features:
//   - Convert system input (s-expression) to params
//   - Build formatted tokens (for output) from params
template< typename Params >
struct KeyValueParser {

    using ParserFunc = std::function< void(Params&, const SExpr::ListType&) >;
    using TokenBuildFunc = std::function< std::list< ConfigFileToken >(const Params&) >;

    using ParserFuncDict = std::map< std::string, ParserFunc >;
    using TokenBuildList = std::vector< TokenBuildFunc >;

    // Two data structures are required for the parser to work.
    //   - A dictionary to match keywords to the appropriate function for
    //     parsing the arguments
    //   - A list which instructs how an input file should be generated from
    //     system params
    ParserFuncDict dict;
    TokenBuildList tokenBuildList;

    // Handy functions
    //---------------------------------

    // Parse or print a specific single-valued parameter.
    //
    // Template parameters
    //   - LocateParam:
    //     (Params&)       -> T&,       AND
    //     (const Params&) -> const T&
    template< typename LocateParam >
    void addSingleArg(
        std::string name,
        LocateParam&& locate
    ) {
        using namespace std;

        addStringArgs(
            name,
            [name, locate](Params& params, const vector<string>& lineVector) {
                auto& param = locate(params);
                using T = decay_t< decltype(param) >;

                if (lineVector.size() != 2) {
                    LOG(ERROR) << name << " must have exactly one value.";
                    throw runtime_error("Invalid argument.");
                }

                auto& arg = lineVector[1];

                if constexpr(is_integral_v< T > || is_floating_point_v< T >) {
                    const auto [p, ec] = from_chars(arg.data(), arg.data() + arg.size(), param);
                    if(ec != std::errc()) {
                        LOG(ERROR) << name << " argument invalid: " << arg;
                        throw runtime_error("Invalid argument.");
                    }
                }
                else if constexpr(is_same_v< T, string >) {
                    param = arg;
                }
            },
            [name, locate](const Params& params) {
                auto& param = locate(params);
                using T = decay_t< decltype(param) >;

                if constexpr(is_integral_v< T > || is_floating_point_v< T >) {
                    return vector<string> { to_string(param) };
                }
                else if constexpr(is_same_v< T, string >) {
                    return vector<string> { param };
                }
            }
        );
    }

    // Template parameters
    //   - FuncParse: void(Params&, const vector<string>&), including key
    //   - FuncBuild:
    //     vector<string>(const Params&), excluding key, OR
    //     vector<vector<string>>(const Params&), excluding key (for repeating items)
    template<
        typename FuncParse,
        typename FuncBuild
    >
    void addStringArgs(
        std::string name,
        FuncParse&& funcParse,
        FuncBuild&& funcBuild
    ) {
        using namespace std;

        return addStringArgsWithAliases(
            move(name),
            {},
            forward<FuncParse>(funcParse),
            forward<FuncBuild>(funcBuild)
        );
    }

    // With alias
    // Template parameters
    //   - FuncParse: void(Params&, const vector<string>&), including key
    //   - FuncBuild:
    //     vector<string>(const Param&), excluding key, OR
    //     vector<vector<string>>(const Param&), excluding key (for repeating items)
    template<
        typename FuncParse,
        typename FuncBuild
    >
    void addStringArgsWithAliases(
        std::string name,
        std::vector< std::string > aliases,
        FuncParse&& funcParse,
        FuncBuild&& funcBuild
    ) {
        using namespace std;

        return addArgsWithAliases(
            name,
            move(aliases),
            [funcParse] (Params& params, const SExpr::ListType& keyAndArgs) {
                funcParse(params, getStringVector(keyAndArgs));
            },
            [funcBuild, name] (const Params& params) {

                if constexpr(is_same_v< invoke_result_t< FuncBuild, const Params& >, vector<string> >) {

                    vector<string> tempResult = funcBuild(params);

                    list< ConfigFileToken > res;
                    for(auto& s : tempResult) {
                        res.push_back(ConfigFileToken::makeString(move(s)));
                    }

                    return res;

                } else {
                    vector<vector<string>> tempResult = funcBuild(params);

                    vector<list< ConfigFileToken >> res;
                    for(auto& eachVS : tempResult) {
                        res.emplace_back();
                        for(auto& s : eachVS) {
                            res.back().push_back(ConfigFileToken::makeString(move(s)));
                        }
                    }

                    return res;
                }
            }
        );
    }

    // With alias
    // Template parameters
    //   - FuncParse: (Params&, const SExpr::ListType&) -> void, including key
    //   - FuncBuild: (const Params&) -> list< ConfigFileToken >, excluding key
    template<
        typename FuncParse,
        typename FuncBuild
    >
    void addArgsWithAliases(
        std::string name,
        std::vector< std::string > aliases,
        FuncParse&& funcParse,
        FuncBuild&& funcBuild
    ) {
        using namespace std;

        const auto insertResult = dict.insert({ name, forward<FuncParse>(funcParse) });
        // Insert for aliases
        for(auto& alias : aliases) {
            dict.insert({ move(alias), insertResult.first->second });
        }

        tokenBuildList.push_back(
            [funcBuild, name] (const Params& params) {

                list< ConfigFileToken > res;

                if constexpr(
                    is_same_v< invoke_result_t< FuncBuild, const Params& >, list<ConfigFileToken> >
                ) {
                    // single entry
                    res.push_back(ConfigFileToken::makeDefault(ConfigFileToken::Type::parenthesisLeft));
                    res.push_back(ConfigFileToken::makeString(name));

                    res.splice(res.end(), funcBuild(params));

                    res.push_back(ConfigFileToken::makeDefault(ConfigFileToken::Type::parenthesisRight));
                    res.push_back(ConfigFileToken::makeDefault(ConfigFileToken::Type::lineBreak));

                } else {
                    // multiple entries
                    vector<list<ConfigFileToken>> tempResult = funcBuild(params);
                    for(auto& eachList : tempResult) {
                        res.push_back(ConfigFileToken::makeDefault(ConfigFileToken::Type::parenthesisLeft));
                        res.push_back(ConfigFileToken::makeString(name));

                        res.splice(res.end(), move(eachList));

                        res.push_back(ConfigFileToken::makeDefault(ConfigFileToken::Type::parenthesisRight));
                        res.push_back(ConfigFileToken::makeDefault(ConfigFileToken::Type::lineBreak));
                    }
                }

                return res;
            }
        );
    }

    void addComment(std::string comment) {
        // The comment must start with valid comment specifiers such as ';'
        tokenBuildList.push_back(
            [comment{ std::move(comment) }] (const Params&) {
                std::list< ConfigFileToken > res {
                    ConfigFileToken {
                        ConfigFileToken::Type::comment,
                        comment
                    }
                };
                res.push_back(ConfigFileToken::makeDefault(ConfigFileToken::Type::lineBreak));
                return res;
            }
        );
    }
    void addEmptyLine() {
        tokenBuildList.push_back(
            [] (const Params&) {
                return std::list< ConfigFileToken > {
                    ConfigFileToken::makeDefault(ConfigFileToken::Type::lineBreak)
                };
            }
        );
    }

};

// What to do when the key is not in the parser function dictionary.
enum class KeyValueParserUnknownKeyAction {
    ignore, warn, error
};

// Parse an s-expr as a key-value pair, using the parser func dictionary.
template< typename Params >
inline void parseKeyValue(
    Params&                                                params,
    const SExpr&                                           se,
    const typename KeyValueParser<Params>::ParserFuncDict& dict,
    KeyValueParserUnknownKeyAction                         unknownKeyAction = KeyValueParserUnknownKeyAction::ignore
) {
    using namespace std;
    using SES = SExpr::StringType;
    using SEL = SExpr::ListType;

    // se.data must be a list type, or an exception will be thrown
    // se contains the (key arg1 arg2 ...) data
    // The key must be a string type, or an exception will be thrown

    // Here the assignments are by value because of possible temporary values on the right hand side.
    const SES key = get<SES>(car(se).data);
    SEL keyAndArgs = get<SEL>(se.data);

    // Check against the parsing dictionary
    if(auto it = dict.find(key); it == dict.end()) {
        switch(unknownKeyAction) {
            case KeyValueParserUnknownKeyAction::ignore:
                break;
            case KeyValueParserUnknownKeyAction::warn:
                LOG(WARNING) << "In the input file, "
                    << key << " cannot be recognized.";
                break;
            case KeyValueParserUnknownKeyAction::error:
                LOG(ERROR) << "In the input file, "
                    << key << " cannot be recognized.";
                throw runtime_error("Unknown key in parser");
        }
    }
    else {
        // Execute the settings
        (it->second)(params, move(keyAndArgs));
    }
}

// Parse a list of s-expr as key-value pairs, using the parser func dictionary.
template< typename Params >
inline void parseKeyValueList(
    Params&                                                params,
    const SExpr&                                           se,
    const typename KeyValueParser<Params>::ParserFuncDict& dict,
    KeyValueParserUnknownKeyAction                         unknownKeyAction = KeyValueParserUnknownKeyAction::ignore
) {
    // se.data must be a list type, or an exception will be thrown
    for(const SExpr& eachList : get<SExpr::ListType>(se.data)) {
        parseKeyValue(params, eachList, dict, unknownKeyAction);
    }
}

// Parsing an s-expr as key-value pairs.
template< typename Params >
inline void parseKeyValueList(
    Params&                        params,
    const SExpr&                   se,
    const KeyValueParser<Params>&  parser,
    KeyValueParserUnknownKeyAction unknownKeyAction = KeyValueParserUnknownKeyAction::ignore
) {
    return parseKeyValueList(params, se, parser.dict, unknownKeyAction);
}

// Build the token list for key-value inputs.
template< typename Params >
inline auto buildTokens(
    const Params&                                          params,
    const typename KeyValueParser<Params>::TokenBuildList& tokenBuildList
) {
    std::list< ConfigFileToken > res;

    for(const auto& eachTokenBuild : tokenBuildList) {
        res.splice(res.end(), eachTokenBuild(params));
    }
    return res;
}
template< typename Params >
inline auto buildTokens(
    const Params&                 params,
    const KeyValueParser<Params>& parser
) {
    return buildTokens(params, parser.tokenBuildList);
}


//-----------------------------------------------------------------------------
// The actual parsers dealing with all medyan system parameters.
//-----------------------------------------------------------------------------

struct SystemParser {
    // Preferably, for one input file, only one parser is needed.
    // However, currently, some parsing rules require the results from parsed
    // params in the same file, and therefore multiple parsing passes are
    // needed.
    // So it is natural to write each set of parsing rules in a separate
    // parser.
    KeyValueParser< SimulConfig >
        headerParser,
        geoParser,
        boundParser,
        mechParser,
        chemParser,
        dyRateParser,
        initParser;

    SystemParser() {
        initInputHeader();
        initGeoParser();
        initBoundParser();
        initMechParser();
        initChemParser();
        initDyRateParser();
        initInitParser();
    }

    void parseInput(SimulConfig& conf, std::string input) const {
        const auto se =
            lexConfigTokens(
                tokenizeConfigFile(
                    std::move(input)));

        parseKeyValueList(conf, se, geoParser);
        geoPostProcessing(conf);

        parseKeyValueList(conf, se, boundParser);
        boundPostProcessing(conf);

        parseKeyValueList(conf, se, mechParser);

        parseKeyValueList(conf, se, chemParser);
        chemPostProcessing(conf);

        parseKeyValueList(conf, se, dyRateParser);

        parseKeyValueList(conf, se, initParser);
    }

    void outputInput(std::ostream& os, const SimulConfig& conf) const {
        std::list< ConfigFileToken > tokens;
        tokens.splice(tokens.end(), buildTokens(conf, headerParser));
        tokens.splice(tokens.end(), buildTokens(conf, geoParser));
        tokens.splice(tokens.end(), buildTokens(conf, boundParser));
        tokens.splice(tokens.end(), buildTokens(conf, mechParser));
        tokens.splice(tokens.end(), buildTokens(conf, chemParser));
        tokens.splice(tokens.end(), buildTokens(conf, dyRateParser));
        tokens.splice(tokens.end(), buildTokens(conf, initParser));

        outputConfigTokens(os, tokens);
    }

    // Add parsing rules
    void initInputHeader();
    void initGeoParser();
    void initBoundParser();
    void initMechParser();
    void initChemParser();
    void initDyRateParser();
    void initInitParser();

    // Post processing and validation
    void geoPostProcessing(SimulConfig&) const;
    void boundPostProcessing(SimulConfig&) const;
    void chemPostProcessing(SimulConfig&) const;

};

struct ChemistryParser {
    // Uses SimulConfig instead of ChemistryData, because the parsing uses
    // information of other parameters
    KeyValueParser< SimulConfig > chemDataParser;

    // State: this parser is not context-free
    // Remember species names to keep track of duplicate names
    std::set<std::string> allSpeciesNames;

    ChemistryParser() {
        initChemDataParser();
    }

    void parseInput(SimulConfig& conf, std::string input) {
        resetStates();

        const auto se =
            lexConfigTokens(
                tokenizeConfigFile(
                    std::move(input)));

        parseKeyValueList(conf, se, chemDataParser);
    }

    void outputInput(std::ostream& os, const SimulConfig& conf) const {
        outputConfigTokens(os, buildTokens(conf, chemDataParser));
    }

    // Add parsing rules
    void resetStates() {
        allSpeciesNames.clear();
    }
    void initChemDataParser();
};

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

/// Used to parse initial Filament information, initialized by the Controller.
struct FilamentParser {
  
    /// Reads filament input file. Returns a vector of tuples containing
    /// filament type and positions (start and end points).
    /// @note - Does not check for coordinate correctness.
    static FilamentData readFilaments(std::istream&);
};

/// Used to parse initial membrane vertex and neighbor information, initialized by the Controller.
struct MembraneParser {

    struct MembraneInfo {
        using coordinate_type = mathfunc::Vec< 3, floatingpoint >;
        std::vector< coordinate_type > vertexCoordinateList;
        std::vector< std::array< size_t, 3 > > triangleVertexIndexList;
    };

    /// Reads membrane vertex input file.
    /// @note - Does not check for coordinate correctness.
    static std::vector<MembraneInfo> readMembranes(std::istream&);
};

/// Used to parse initial Bubble information, initialized by the Controller.
struct BubbleParser {    
    /// Reads bubble input file. Returns a vector of tuples containing
    /// bubble type and position.
    /// @note - Does not check for coordinate correctness.
    static BubbleData readBubbles(std::istream&);
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
