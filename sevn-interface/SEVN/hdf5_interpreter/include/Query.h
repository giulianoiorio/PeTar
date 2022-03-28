//
// Created by mario on 04/02/19.
//

#ifndef H5INTER_QUERY_H
#define H5INTER_QUERY_H

#include <iostream>
#include <processor.h>

class Header : public Command{

public:
    Header() {
        CommandProcessor::Instance().Register( "header", this, 0); //name, command, nargs
    }

    void Execute();
private:

};

class Dataset : public Command{

public:
    Dataset() {
        CommandProcessor::Instance().Register( "dataset", this, 0); //name, command, nargs
    }

    void Execute();
private:

};


#endif //H5INTER_QUERY_H
