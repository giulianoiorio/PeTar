//
// Created by mario on 03/02/19.
//

#ifndef H5INTER_GLOBAL_H
#define H5INTER_GLOBAL_H

#include <LoadSStars.h>
#include <LoadDataset.h>
#include <Quit.h>
#include <Plot.h>
#include <Query.h>

extern LoadSStars _load;
extern LoadDataset _loaddataset;
extern QuitCommand _quit;
extern Plot _plot;
extern Header _header;
extern Dataset _dataset;

#endif //H5INTER_GLOBAL_H
