//
// Created by mario on 03/02/19.
//


#include <processor.h>
#include <global.h>

//class-command list

LoadSStars _load;
LoadDataset _loaddataset;
QuitCommand _quit;
Plot _plot;
Header _header;
Dataset _dataset;

int main() {

    CommandProcessor::Instance().Run();

}