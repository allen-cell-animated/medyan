#include <algorithm> // transform
#include <string>
#include <vector>

#include "catch2/catch.hpp"

#include "Util/Io/CmdParse.hpp"

namespace {

void parseForCommand(cmdparse::Command& cmd, std::vector< std::string > args) {
    using namespace std;
    size_t s = args.size();
    vector< char* > arr(s);
    for(size_t i = 0; i < s; ++i) arr[i] = &args[i][0];
    cmd.parse(s, arr.data());
}

} // namespace

TEST_CASE("Command line argument parsing", "[CmdParse]") {
    using namespace std;
    using namespace cmdparse;

    Command cmd("cmd", "desc");

    SECTION("Positional arguments") {
        SECTION("No positional argument") {
            CHECK_THROWS_AS(parseForCommand(cmd, {"cmd", "arg"}), ParsingError);
        }

        SECTION("1 positional argument") {
            double arg = 1.0;

            SECTION("Required argument") {
                cmd.addPosArgForVar("arg", "desc", true, arg);
                CHECK_THROWS_AS(parseForCommand(cmd, {"cmd"}), ValidationError);
            }
            SECTION("Optional argument") {
                cmd.addPosArgForVar("arg", "desc", false, arg);
                parseForCommand(cmd, {"cmd", "--", "-2.0"});
                CHECK(arg == Approx(-2.0));
            }
        }

        SECTION("Positional argument list") {
            vector< float > args;

            SECTION("Required argument") {
                cmd.addPosArgForVector("args", "desc", true, args);
                CHECK_THROWS_AS(parseForCommand(cmd, {"cmd"}), ValidationError);
            }
            SECTION("Optional argument") {
                cmd.addPosArgForVector("args", "desc", false, args);

                SECTION("Supply none") {
                    parseForCommand(cmd, {"cmd"});
                    CHECK(args.size() == 0);
                }
                SECTION("Supply multiple") {
                    parseForCommand(cmd, {"cmd", "1.0", "2"});
                    CHECK(args.size() == 2);
                    CHECK(args[1] == Approx(2.0f));
                }
            }
        }

        SECTION("Mixing positional arguments") {
            // TODO
        }

    }
}
