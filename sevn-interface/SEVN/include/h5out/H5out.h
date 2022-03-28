//
// Created by mario on 30/01/19.
//

#ifndef SEVN_H5OUT_H
#define SEVN_H5OUT_H

#include <H5Cpp.h>
#include <H5File.h>
#include <vector>
#include <string>
#include <utilities.h>
#include <sevnlog.h>
#include <errhand.h>
#include <zconf.h>

using sevnstd::SevnLogging;

namespace sevnstd {
    class H5out {

    public:
        H5out(const std::string &output_folder="output", const std::string &fname="stars.h5"){

            filename_h5 = utilities::gen_filename(output_folder,fname); //this should be OpenMP-safe
            DEBUG_LOG("f5name = "<<filename_h5);
            //utilities::wait();

            groupname_h5 = "/stars";
            header_name_h5 = groupname_h5 + "/header";
            header_printed = false;
        }

        void print(std::vector<std::string> &columns, std::vector<std::vector<double>> &printmatrix, const std::string &_name);
        void print_header(std::vector<std::string> &columns);
        inline void set_filename(std::string a) {filename_h5=a;};
        inline std::string get_filename() { return filename_h5;};

    private:

        std::string filename_h5, groupname_h5, header_name_h5;
        SevnLogging svlog;
        bool header_printed;

        H5::H5File file;
        H5::DataSpace dataspace;
        H5::DSetCreatPropList plist;
        H5::DataSet dataset;


        hsize_t dim[2];
        hsize_t chunk_dim[2];

        size_t col_size;


    };
}


#endif //SEVN_H5OUT_H
