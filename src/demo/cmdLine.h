/**
 * @file cmdLine.h
 * @brief Command line option parsing
 * @author Pascal Monasse
 * 
 * Copyright (c) 2012-2017 Pascal Monasse
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CMDLINE_H
#define CMDLINE_H

#include <vector>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <string>
#include <sstream>
#include <cassert>

/// Base class for option/switch
class Option {
public:
    char c; ///< Option letter (eg 's' for -s), 0 if you want only a long name
    bool used; ///< Does the command line use that option?
    std::string longName; ///< Optional long name (eg "switch" for --switch)
    std::string desc; ///< Description
    std::string section; ///< Options can be sorted by section

    /// Constructor with short name/long name
    Option(char d, std::string name)
    : c(d), used(false), longName(name) {}
    virtual ~Option() {}
    /// Output the option identifier, short and/or long name
    virtual void print(std::ostream& str) const {
        if(c)
            str << '-' << c << (longName.empty()? "": ", ");
        if(! longName.empty())
            str << "--" << longName;
    }
    /// Set the description of an option, same as assigning field \a desc
    Option& doc(const std::string& description)
    { desc=description; return *this;}
    virtual bool check(int& argc, char* argv[])=0; ///< Option found at argv[0]?
    virtual Option* clone() const=0; ///< Copy
    /// Print value of option
    virtual std::string printValue() const { return std::string(); }
};

/// Option on/off is called a switch
class OptionSwitch : public Option {
public:
    /// Constructor with short name/long name (optional)
    OptionSwitch(char c, std::string name="")
    : Option(c,name) {}
    /// Find switch in argv[0]
    bool check(int& argc, char* argv[]) {
        if(std::string("-")+c==argv[0] ||
           (!longName.empty() && std::string("--")+longName==argv[0])) {
            used = true;
            std::rotate(argv, argv+1, argv+argc);
            argc -= 1;
            return true;
        } else if(std::string(argv[0]).find(std::string("-")+c)==0) {
            used = true; // Handle multiple switches in single option
            std::rotate(argv[0]+1, argv[0]+2,
                        argv[0]+std::string(argv[0]).size()+1);
            return true;
        }
        return false;
    }
    /// Copy
    Option* clone() const {
        return new OptionSwitch(*this);
    }
};

/// Option with an argument of type T, which must be readable by operator>>
template <class T>
class OptionField : public Option {
public:
    /// Constructor. The result with be stored in variable \c field.
    OptionField(char c, T& field, std::string name="")
    : Option(c,name), _field(field) {}
    /// Find option in argv[0] and argument in argv[1]. Throw an exception
    /// (type std::string) if the argument cannot be read.
    bool check(int& argc, char* argv[]) {
        std::string param; int arg=0;
        if(std::string("-")+c==argv[0] ||
           (!longName.empty() && std::string("--")+longName==argv[0])) {
            if(argc<=1)
                throw std::string("Option ")
                    +argv[0]+" requires argument";
            param=argv[1]; arg=2;
        } else if(std::string(argv[0]).find(std::string("-")+c)==0) {
            param=argv[0]+2; arg=1;
        } else if(!longName.empty() &&
                  std::string(argv[0]).find(std::string("--")+longName+'=')==0){
            size_t size=(std::string("--")+longName+'=').size();
            param=std::string(argv[0]).substr(size); arg=1;
        }
        if(arg>0) {
            if(! read_param(param))
                throw std::string("Unable to interpret ")
                    +param+" as argument of "+argv[0];
            used = true;
            std::rotate(argv, argv+arg, argv+argc);
            argc -= arg;
            return true;
        }
        return false;
    }
    /// Decode the string as template type T
    bool read_param(const std::string& param) {
        std::stringstream str(param); char unused;
        return !((str >> _field).fail() || !(str>>unused).fail());
    }
    /// Indicate that an argument is required
    void print(std::ostream& str) const {
        Option::print(str);
        str << (longName.empty()? ' ': '=') << "ARG";
    }
    std::string printValue() const {
        std::stringstream s;
        s << _field;
        return s.str();
    }
    /// Copy
    Option* clone() const {
        return new OptionField<T>(*this);
    }
private:
    T& _field; ///< Reference to variable where to store the value
};

/// Template specialization to declare a switch like an option, storing result
/// in the variable.
template <>
inline bool OptionField<bool>::check(int& argc, char* argv[]) {
    bool res = OptionSwitch(c,longName).check(argc, argv);
    if(res)
        _field = true;
    return res;
}

/// Specialisation for a switch, no argument to print.
template<>
void OptionField<bool>::print(std::ostream& str) const {
    Option::print(str);    
}

/// Specialisation for a switch, no argument to print.
template<>
std::string OptionField<bool>::printValue() const {
    return Option::printValue();
}

/// Template specialization to be able to take parameter including space.
/// Generic method would do >>_field (stops at space) and signal unused chars.
template <>
inline bool OptionField<std::string>::read_param(const std::string& param) {
    _field = param;
    return true;
}

/// New switch option
OptionSwitch make_switch(char c, std::string name="") {
    return OptionSwitch(c, name);
}

/// New option with argument.
template <class T>
OptionField<T> make_option(char c, T& field, std::string name="") {
    return OptionField<T>(c, field, name);
}

/// Utility function to sort options by section
static bool order_by_section(Option*const o1, Option*const o2) {
    return (o1->section < o2->section);
}

/// Command line parsing.
///
/// This class parses a command line into options and positional arguments. Its
/// purpose is similar to
///    - GNU getopt (https://www.gnu.org/software/libc/manual/html_node/Getopt.html)
///    - Boost.Program_options (http://www.boost.org/doc/libs/1_58_0/doc/html/program_options.html)
///
/// but should be simpler to use (portable, single header file). An example:
/// \code
/// #include "cmdLine.h"
/// int main(int argc, char* argv[]) {
///     CmdLine cmd;
///     bool help; int level=0;
///     cmd.add( make_option('h', help, "help").doc("usage info") );
///     cmd.add( make_option('l', level, "level").doc("integer level") );
///     try {
///         cmd.process(argc, argv);
///     } catch(std::string str) {
///         std::cerr << "Error: " << str << std::endl;
///         return 1;
///     }
///     // After calling process, positional args are argv[1]...argv[argc-1]
///     if(help)
///         std::cout << "Usage: " << argv[0] << '\n' << cmd;
///     std::cout << "Level argument is " << level << std::endl;
///     return 0;
/// }
/// \endcode
/// An instance of \c CmdLine is created, to which options are appended with
/// method add. Two options are proposed here, the first one invoked by "-h"
/// or "--help", the second one by "-l 2", '-l2', "--level 2", or "--level=2".
/// \c make_option is the preferred way to create an option, because it detects
/// the type of its variable and behaves appropriately: a \c bool indicates a
/// switch (option without argument), any other type a regular option requiring
/// a value. Non-standard types can be used, provided `operator<<` with
/// an `std::ostream` and `operator>>` with an `std::istream` are defined.
/// \note An option having no long name is created by omitting it:
/// <tt>make_option('l',help)</tt>. An option can have \a only a long name by putting
/// 0 as its \c char identifier: `make_option(0,help,"help")`.
/// \warning The single letter option must be a \c char, between single
/// quotes, not a \c char* in double quotes:
/// \code
/// make_option("h",help); // Does not compile: "h" (char*) must be 'h' (char)
/// \endcode
/// It is possible but not a good idea to have an option whose single
/// letter identifier is a digit. Indeed, the minus sign could be interpreted
/// either as the introduction of an option or as the sign of a number. The
/// method \c process tries the former, but failing that accepts it as a
/// positional argument if it can be interpreted as a number.
/// \code
/// bool b;
/// cmd.add( make_option('1',b) ); // Bad idea to use a digit as identifier
/// cmd.process(argc,argv);
/// std::cout << argc << std::endl;
/// 
/// $ prog -1
/// 1
/// $ prog -3
/// 2
/// \endcode
/// The first run identifies `-1` as the option, removing it from the arguments,
/// so only `argv[0]` remains; the second one interprets '-2' as a positional
/// argument, it remains in `argv[1]`.
///
/// To decode the command line, the method \c process is called. It extracts
/// from \c argc and \c argv the options, modifying these variables so that only
/// positional arguments are left. If an error occurs (unknown option, like "-k"
/// or missing value after "-l"), an exception, in the form of an
/// \c std::string, is raised. Our example catches such an exception to indicate
/// the problem to the user. The variable \c level has a default value of 0,
/// which is not changed if the corresponding option is not invoked. The
/// boolean variable \c help is set to \c true only when the switch is used.
///
/// The usage information is obtained by `std::cout << cmd`. This lists the
/// options and their description (set with the \c doc function). Several
/// fields of \c cmd control the display format:
///   - \c prefixDoc: an \c std::string to display as prefix before each option,
/// such as a few spaces or a tabulation.
///   - \c alignDoc: an integer indicating the column number where option
/// description should normally start if the length of option name allows.
///   - \c showDefaults: a boolean (\c true by default) controlling whether
/// default values of options should be displayed after their description. This
/// does not concern switches (options without arguments, controlled by a
/// variable of type \c bool).
///   - \c section: an \c std::string allowing to group options by sections, see
/// below.
///
/// \warning Notice that if \c showDefault is \c true, it is important that all
/// variables holding option values be initialized. Variables of type \c bool
/// are not concerned, since they are automatically set to \c false when used
/// in an option.
///
/// Options can be grouped by section, each holding a \c section (type
/// \c std::string) field. Instead of specifying the section for each option,
/// it is more convenient to set \c CmdLine::section, which affects the section
/// of each subsequent option added. The default section is the empty string.
/// If options are not all in the same single section, printing the object
/// \c CmdLine has these effects:
///    1. Options are sorted by section.
///    2. The name of each section is output before its options.
///    3. Section are printed by lexicographic order of their name.
///
/// To override such behavior and get back full control, the user can create
/// a new \c CmdLine filtering a specific section through the dedicated
/// constructor. For example:
/// \code
/// CmdLine cmd;
/// bool help; int level=0;
/// cmd.section = "General";
/// cmd.add( make_option('h', help, "help").doc("usage info") );
/// cmd.section = "Specific";
/// cmd.add( make_option('l', level, "level").doc("integer level") );
/// CmdLine cmdG(cmd, "General"), cmdS(cmd, "Specific");
/// std::cout << "General:\n" << cmdG << "Specific:\n" << cmdS;
/// \endcode
/// The last two lines could be compressed in
/// \code
/// std::cout << "General:\n"  << CmdLine(cmd,"General")
///           << "Specific:\n" << CmdLine(cmd,"Specific");
/// \endcode
class CmdLine {
    std::vector<Option*> opts;
public:
    std::string prefixDoc; ///< For example, a tabulation for each line of doc
    int alignDoc; ///< Column where option description starts
    bool showDefaults; ///< Show default values of options
    std::string section; ///< Section where to put next added options

    /// Constructor
    CmdLine(): alignDoc(0), showDefaults(true) {}
    /// Constructor copying only options from a given section
    CmdLine(const CmdLine& cmd, const std::string& sect) {
        *this = cmd;
        opts.clear();
        for(std::vector<Option*>::const_iterator it=cmd.opts.begin();
            it!=cmd.opts.end(); ++it)
            if((*it)->section == sect)
                opts.push_back( (*it)->clone() );
    }
    /// Destructor
    ~CmdLine() {
        std::vector<Option*>::iterator it=opts.begin();
        for(; it != opts.end(); ++it)
            delete *it;
    }
    /// Add an option
    void add(const Option& opt) {
        Option* c = opt.clone();
        c->section = section;
        opts.push_back(c);
    }
    /// Parse line acting as a filter removing options from the command line.
    void process(int& argc, char* argv[]) {
        std::vector<Option*>::iterator it=opts.begin();
        for(; it != opts.end(); ++it)
            (*it)->used = false;
        for(int i=1; i<argc;) {
            if(std::string("--")==argv[i]) { // "--" means stop option parsing
                std::rotate(argv+i, argv+i+1, argv+argc);
                -- argc;
                break;
            }
            bool found=false; // Find option
            for(it=opts.begin(); it != opts.end(); ++it) {
                int n = argc-i;
                found = (*it)->check(n, argv+i);
                if(found) {
                    argc = n+i;
                    break;
                }
            }
            if(! found) { // A negative number is not an option
                if(std::string(argv[i]).size()>1 && argv[i][0] == '-') {
                    std::istringstream str(argv[i]);
                    float v;
                    if(! (str>>v).eof())
                        throw std::string("Unrecognized option ")+argv[i];
                }
                ++i;
            }
        }
        // Order options by section, useful if print is called after
        std::stable_sort(opts.begin(), opts.end(), order_by_section);
    }
    /// Output options.
    void print(std::ostream& str) const {
        if(opts.empty()) return;
        std::string prevSection = opts.back()->section;
        bool showSection = (opts.front()->section != prevSection);
        std::vector<Option*>::const_iterator it=opts.begin();
        for(; it != opts.end(); ++it) {
            if(showSection && prevSection!=(*it)->section)
                str << (prevSection=(*it)->section) << std::endl;
            std::stringstream ss;
            ss << prefixDoc;
            (*it)->print(ss);
            ss << ' ';
            std::string d = ss.str();
            str << d;
            if((int)d.size() < alignDoc)
                std::fill_n(std::ostream_iterator<char>(str),
                            alignDoc-d.size(), ' ');
            str << (*it)->desc;
            if(showDefaults && !(*it)->printValue().empty())
                str << " (" << (*it)->printValue() << ')';
            str << std::endl;
        }
    }
    /// Was the option used in last parsing?
    bool used(char c) const {
        std::vector<Option*>::const_iterator it=opts.begin();
        for(; it != opts.end(); ++it)
            if((*it)->c == c)
                return (*it)->used;
        assert(false); // Called with non-existent option, probably a bug
        return false;
    }
};

/// Output possible options.
std::ostream& operator<<(std::ostream& str, const CmdLine& cmd) {
    cmd.print(str);
    return str;
}

#endif
