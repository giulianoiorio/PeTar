//
// Created by mario on 04/02/19.
//

#ifndef H5INTER_PLOT_H
#define H5INTER_PLOT_H

#include <processor.h>

class Plot : public Command{

public:
    Plot() {
        CommandProcessor::Instance().Register("plot", this, -1); //-1 means variable number of arguments

    }

    void Execute();

private:
    size_t xpos;
    std::vector<size_t> ypos;


};


#endif //H5INTER_PLOT_H
