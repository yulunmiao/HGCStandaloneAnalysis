/* 
   Command line argument parser. Supports:
   - single arguments with any value
   - vector arguments with any value
   - single arguments with limited user-defined values
   - vector arguments with limited user-defined values
   - boolean flags

   All arguments must be preceded by `--`.

   Vector arguments are passed as follows:
   < command --vector_argument length val1 val2 val3 ... >
   where `length` refers to the number of `val` elements in the vector.
   For example: < plot.cc --energies 3 30 50 100 >

   The arguments must be defined in derived classes, by
   overriding the following methods:
   - set_args_options()
   - set_args_final()
   Do not forget to declare the command line arguments as private 
   members in the derived class, and their getters.

   See `src/InputParsePlotEGReso.cc` for an example.

   Author: Bruno Alves, CERN (bruno.alves_AT_cern.ch)
*/
#include<iostream>
#include<unordered_map>
#include<string>
#include<vector>
#include<algorithm>

class InputParserBase {
 public:
  explicit InputParserBase(int argc, char** argv): argc_(argc), argv_(argv) {}

  int argc_;
  char** argv_;

  virtual void set_args_options_() = 0;
  virtual void set_args_final_() = 0;
  
  std::unordered_map< std::string, std::vector<std::string> > valid_args_; //arguments with limited options
  std::unordered_map< std::string, std::vector<std::string> > valid_args_v_; //vector arguments with limited options
  std::vector<std::string> free_args_; //any argument allowed
  std::vector<std::string> free_args_v_; //any vector argument allowed
  std::vector<std::string> optional_args_; //optional arguments (boolean flags)

  std::vector<std::string> required_args_; //all required arguments (contains the above)

  std::unordered_map<std::string, std::string> chosen_args_;
  std::unordered_map<std::string, std::vector<std::string>> chosen_args_v_;

  void run();
  bool are_required_args_present_();
  bool are_args_supported_();
  void check_required_args_are_present_();
  bool is_arg_absent_(char*);
  bool help_();
  void parse_chosen_args_();
  void print_vector_elements_(const std::vector<std::string>&);
  void set_chosen_args_();
  bool was_help_invoked_();
};
