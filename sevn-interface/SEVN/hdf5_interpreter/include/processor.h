//
// Created by mario on 03/02/19.
//

#ifndef H5INTER_PROCESSOR_H
#define H5INTER_PROCESSOR_H

#include <iostream>
#include <map>

#include <algorithm>
#include <readline/readline.h>
#include <readline/history.h>
#include <command.h>

class CommandProcessor {

    public:

        typedef std::map <std::string, Command *> MapType;

        static CommandProcessor & Instance() {
            static CommandProcessor cp;
            return cp;
        }

        void Run() {
            char* line;
            std::string cmdline;

            for(;;) {

                line = readline("H5interpreter >>> ");

                try {
                    add_history(line);

                    cmdline = std::string(line); free(line);
                    std::vector<std::string> args;

                    Command *_command = Cleaner(cmdline, args);
                    _command->GetArguments(args);
                    _command->Execute();
                }
                catch( const std::string &msg ) {
                    std::cout << msg << "\n";
                }
            }
        }


        void Register( const std::string & name, Command * cmd, const int nargs) {
            cmdMap.insert( std::make_pair(name, cmd) );
            cmd->set_nargs(nargs);
        }


        Command * Cleaner(const std::string & cmd, std::vector<std::string> & args) const {

           // std::cout<<" insert command = "<<cmd<<std::endl;

            size_t pos1 = cmd.find('(');
            size_t pos2 = cmd.find(')');
            std::string cmd_cleaned = cmd, _args;

            cmd_cleaned.erase(remove_if(cmd_cleaned.begin(), cmd_cleaned.end(), ::isspace), cmd_cleaned.end());
            cmd_cleaned.erase(remove(cmd_cleaned.begin(), cmd_cleaned.end(), '\''), cmd_cleaned.end());
            cmd_cleaned.erase(remove(cmd_cleaned.begin(), cmd_cleaned.end(), '\"'), cmd_cleaned.end());


            if(pos1 != std::string::npos) {
                cmd_cleaned = cmd.substr(0, pos1);
                _args = cmd.substr(pos1 + 1, pos2 - pos1 - 1);

               // std::cout << " cleaned command = " << cmd_cleaned << std::endl;
               // std::cout << " aruments = " << _args << std::endl;

                //remove white spaces, remove ', remove "
                _args.erase(remove_if(_args.begin(), _args.end(), ::isspace), _args.end());
                _args.erase(remove(_args.begin(), _args.end(), '\''), _args.end());
                _args.erase(remove(_args.begin(), _args.end(), '\"'), _args.end());

                cmd_cleaned.erase(remove_if(cmd_cleaned.begin(), cmd_cleaned.end(), ::isspace), cmd_cleaned.end());
                cmd_cleaned.erase(remove(cmd_cleaned.begin(), cmd_cleaned.end(), '\''), cmd_cleaned.end());
                cmd_cleaned.erase(remove(cmd_cleaned.begin(), cmd_cleaned.end(), '\"'), cmd_cleaned.end());

              //  std::cout << " aruments cleaned = " << _args << std::endl;


                size_t pos = 0;
                std::string token;
                std::string delimiter = ",";
                while ((pos = _args.find(delimiter)) != std::string::npos) {
                    token = _args.substr(0, pos);
                    args.push_back(token);
                    _args.erase(0, pos + delimiter.length());
                }

                args.push_back(_args);
            }



            return LookupCommand(cmd_cleaned, args);




        }

        Command * LookupCommand( const std::string & cmd_cleaned, const std::vector<std::string> & args ) const {

            MapType::const_iterator it = cmdMap.find( cmd_cleaned );

            if ( it == cmdMap.end() ) {
                const std::string err = "Invalid command: [" + cmd_cleaned + "]";
                throw err;
            }


            if( it->second->get_nargs() != args.size() && it->second->get_nargs() != -1) {
                const std::string err = "Unexpected number of arguments for command [" + it->first + "]. Expected "
                                        + std::to_string(it->second->get_nargs()) + ", given " + std::to_string(args.size());
                throw err;
            }


            return it->second;

        }

    private:

        CommandProcessor() {}
        MapType cmdMap;
    };


#endif //H5INTER_PROCESSOR_H
