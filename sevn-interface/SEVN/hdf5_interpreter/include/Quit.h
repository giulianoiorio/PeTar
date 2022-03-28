//
// Created by mario on 04/02/19.
//

#ifndef H5INTER_QUIT_H
#define H5INTER_QUIT_H

#include <processor.h>

class QuitCommand : public Command {
public:
    QuitCommand() {
        CommandProcessor::Instance().Register( "quit", this, 0 );
        CommandProcessor::Instance().Register( "exit", this, 0 );
        CommandProcessor::Instance().Register( "close", this, 0 );
    }

    void Execute() {
        exit(0);
    }

};

#endif //H5INTER_QUIT_H
