//
// Created by mario on 03/02/19.
//

#ifndef H5INTER_COMMAND_H
#define H5INTER_COMMAND_H

#include <vector>
#include <string>

class Command {

public:
    virtual ~Command() {}
    virtual void Execute(){}
    virtual void GetArguments(const std::vector<std::string> & _args) {args = _args;}

    void set_nargs(const int _nargs){ nargs = _nargs;}
    const int get_nargs() const {return nargs;}

protected:
    std::vector<std::string> args;
    int nargs;

};


#endif //H5INTER_COMMAND_H
