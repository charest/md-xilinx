#ifndef HOST_ARGS_HPP
#define HOST_ARGS_HPP

#include <string>

char* get_cmd_option(char ** begin, char ** end, const std::string & option);

bool cmd_option_exists(char** begin, char** end, const std::string& option);

#endif
