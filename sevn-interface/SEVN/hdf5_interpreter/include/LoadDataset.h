//
// Created by mario on 03/02/19.
//

#ifndef H5INTER_LOADDATASET_H
#define H5INTER_LOADDATASET_H

#include <processor.h>
#include <H5Cpp.h>

class LoadDataset : public Command{

public:
    LoadDataset() {

        loadedheader = false;
        loadeddata = false;
        dataset_name = "";
        CommandProcessor::Instance().Register("load_dataset", this, 2);
    }

    void Execute();

    std::vector<std::string> * get_allheader() {return &allheader;}
    std::vector<double> * get_alldata() {return &alldata;}
    inline bool hdrloaded() {return loadedheader;}
    inline bool dataloaded() {return loadeddata;}
    inline size_t get_xdim() {return dims_out[0];}
    inline size_t get_ydim() {return dims_out[1];}
    inline std::string get_name() {return dataset_name;}

private:
    hsize_t dims_out[2];
    bool loadedheader, loadeddata;
    std::vector<double> alldata;
    std::vector<std::string> allheader;
    std::string dataset_name;

};


#endif //H5INTER_LOADDATASET_H
