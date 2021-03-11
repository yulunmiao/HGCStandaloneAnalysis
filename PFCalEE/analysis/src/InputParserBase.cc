#include "../include/InputParserBase.h"

bool InputParserBase::are_required_args_present_()
{
  for(auto&& req: required_args_)
    {
      bool is_present = false;
      for(int iarg=0; iarg<argc_; ++iarg)
	{
	  if(std::string(argv_[iarg]) == req)
	    is_present = true;
	}
      if(!is_present) {
	std::cout << "You must specify the `" << req << "` argument." << std::endl;
	return 0;
      } 
    }
  return 1;
}

bool InputParserBase::are_args_supported_()
{
  for(int iarg=0; iarg<argc_; ++iarg) {
    if(is_arg_mistyped_(argv_[iarg]))
      {
	std::cout << "You mistyped " << std::string(argv_[iarg]) << "." << std::endl;
	return 0;
      }
    if(is_arg_absent_(argv_[iarg]))
      {
	std::cout << "Argument " << std::string(argv_[iarg]) << " was not recognized." << std::endl;
	std::cout << "The arguments currently supported are:" << std::endl;
	for(auto& elem : valid_args_)
	  std::cout << elem.first << std::endl;
	for(auto& elem : valid_args_v_)
	  std::cout << elem.first << std::endl;
	for(auto& elem : free_args_)
	  std::cout << elem << std::endl;
	for(auto& elem : free_args_v_)
	  std::cout << elem << std::endl;
	for(auto& elem : optional_args_)
	  std::cout << elem << std::endl;
	return 0;
      }
  }
  return 1;
}

void InputParserBase::check_required_args_are_present_()
{
  for(auto&& req: required_args_)
    {
      if( valid_args_.find(req) == valid_args_.end() and
	  valid_args_v_.find(req) == valid_args_v_.end() and
	  std::find(free_args_.begin(), free_args_.end(), req) == free_args_.end() and
	  std::find(free_args_v_.begin(), free_args_v_.end(), req) == free_args_v_.end() ) {
	std::cout << "You forgot to set optional arguments for `" << req << "`.";
	std::exit(1);
      }
    }
}

void InputParserBase::set_optional_args_to_false_()
{
  for(auto&& opt: optional_args_)
    chosen_args_[opt] = "false";
}

//returns true if the passed argument was not defined by the user
bool InputParserBase::is_arg_absent_(char* s)
{
  return ( std::string(s).find("--") != std::string::npos and
	   valid_args_.find(std::string(s)) == valid_args_.end() and
	   valid_args_v_.find(std::string(s)) == valid_args_v_.end() and
	   std::find(free_args_.begin(), free_args_.end(), std::string(s)) == free_args_.end() and
	   std::find(free_args_v_.begin(), free_args_v_.end(), std::string(s)) == free_args_v_.end() and
	   std::find(optional_args_.begin(), optional_args_.end(), std::string(s)) == optional_args_.end() );
}

bool InputParserBase::is_arg_mistyped_(char* s)
{
  if(std::string(s).find("--") == std::string::npos and s[0] == '-')
    return true;
  return false;
}


bool InputParserBase::help_()
{     
  if(argc_ == 1 or was_help_invoked_()) {
    std::cout << "=====================" << std::endl;
    std::cout << "You must specify the following:" << std::endl;
    for(auto& elem : valid_args_) {
      std::string elem2 = elem.first;
      elem2.erase(0,2);
      std::cout << elem2 + ": ";
      print_vector_elements_(elem.second);
    }
    for(auto& elem : valid_args_v_) {
      std::string elem2 = elem.first;
      elem2.erase(0,2);
      std::cout << elem2 + ": ";
      print_vector_elements_(elem.second);
    }
    for(std::string& elem : free_args_) {
      std::string elem2 = elem;
      elem2.erase(0,2);
      std::cout << elem2 + ": required, any choice allowed" << std::endl;
    }
    for(std::string& elem : free_args_v_) {
      std::string elem2 = elem;
      elem2.erase(0,2);
      std::cout << elem2 + ": vector required, any choice allowed" << std::endl;
    }
    for(std::string& elem : optional_args_) {
      std::string elem2 = elem;
      elem2.erase(0,2);    
      std::cout << elem2 + ": optional" << std::endl;
    }
    std::cout << "Note: vector arguments must obbey the "
      "`--arg length val1 val2 val3 ...` structure. " << std::endl;
    std::cout << "=====================" <<  std::endl;
    return 0;
  }
  return 1;
}

void InputParserBase::parse_chosen_args_()
{
  if( !help_() or !are_args_supported_() or !are_required_args_present_())
    std::exit(1);
  set_chosen_args_();
}

//convenience function which prints all the elements in a vector of strings to std::cout
void InputParserBase::print_vector_elements_(const std::vector<std::string>& v)
{
  for(auto& elem : v) {
    if(elem == v.back())
      std::cout << elem << std::endl;
    else
      std::cout << elem << " / ";
  }
}

void InputParserBase::run() {
  set_args_options_();
  check_required_args_are_present_();
  parse_chosen_args_();
  set_args_final_();
}

void InputParserBase::set_chosen_args_()
{
  for(int iarg=1; iarg<argc_; ++iarg)
    {
      std::string argvstr = std::string(argv_[iarg]);
      if( valid_args_.find(argvstr) != valid_args_.end() ) {
	if( std::find(valid_args_[argvstr].begin(), valid_args_[argvstr].end(), std::string(argv_[iarg+1])) == valid_args_[argvstr].end() )
	  {
	    std::cout << "The data type has to be one of the following: ";
	    print_vector_elements_(valid_args_[argvstr]);
	    std::exit(0);
	  }
	else
	  chosen_args_[argvstr] = std::string(argv_[iarg+1]);
      }
      else if(valid_args_v_.find(argvstr) != valid_args_v_.end() ) {
	int idx = iarg+1;
	for(int i(1); i<=std::stoi(std::string(argv_[idx])); ++i)
	  {
	    if( std::find(valid_args_v_[argvstr].begin(), valid_args_v_[argvstr].end(), std::string(argv_[idx+i])) == valid_args_v_[argvstr].end() )
	      {
		std::cout << "The vector data type has to be one of the following: ";
		print_vector_elements_(valid_args_v_[argvstr]);
		std::exit(0);
	      }
	    else
	      chosen_args_v_[argvstr].push_back( std::string(argv_[idx+i]) );
	  }
      }
      else if( std::find(free_args_.begin(), free_args_.end(), argvstr) != free_args_.end() )
	chosen_args_[argvstr] = std::string(argv_[iarg+1]);
      else if( std::find(free_args_v_.begin(), free_args_v_.end(), argvstr) != free_args_v_.end() )
	{
	  int idx = iarg+1;
	  for(int i(1); i<=std::stoi(argv_[idx]); ++i)
	    chosen_args_v_[argvstr].push_back( std::string(argv_[idx+i]) );
	}
      else if(std::find(optional_args_.begin(), optional_args_.end(), std::string(argv_[iarg])) != optional_args_.end()) {
	chosen_args_[argvstr] = "true";
      }
    }
}


/*
  Currently redundant, as the printout message appears as soon as one of the
  arguments is wrong.
*/
bool InputParserBase::was_help_invoked_()
{
  std::array<std::string,3> help_calls = {{"--help", "help", "-h"}};
  for(int iarg(0); iarg<argc_; ++iarg) {
    for(unsigned i(0); i<help_calls.size(); ++i)
      if(std::string(argv_[iarg]).find(help_calls[i]) != std::string::npos)
	return 1;
  }
  return 0;
}
