//
// Created by mario on 30/01/19.
//

#include <H5out.h>
#include <lookup_and_phases.h>
#include <property.h>
#include <BinaryProperty.h>

void sevnstd::H5out::print(std::vector<std::string> &columns, std::vector<std::vector<double>> &printmatrix, const std::string &_name) {

    col_size = columns.size();
    print_header(columns);

    const std::string starname = groupname_h5 + "/" + _name; //name of the dataset

    const size_t xdim = dim[0] = printmatrix.size();
    const size_t ydim = dim[1] = col_size;
    chunk_dim[0] = xdim;
    chunk_dim[1] = ydim;

    if(printmatrix[0].size() != col_size)
        svlog.critical("Mismatch dimension number of columns, amd columns in printmatrix", __FILE__, __LINE__);


    //vector of vector is not contiguous in memory
    //HDF5 needs contiguous memory, thus we need to convert!
    double vec [xdim][ydim];

    for (size_t i = 0; i < xdim; i++) {
        for (size_t j = 0; j < ydim; j++) {
            vec[i][j] = printmatrix[i][j];
            std::cout << vec[i][j] << "  " << std::endl;
        }
    }


    // Create the data space for the dataset.
    dataspace = H5::DataSpace(2, dim);


    plist = H5::DSetCreatPropList();
    plist.setChunk(2, chunk_dim);
    plist.setDeflate(6);

    // Create the dataset
    std::cout<<" Create the dataset"<<std::endl;
    dataset = H5::DataSet(file.createDataSet(starname, H5::PredType::NATIVE_DOUBLE, dataspace, plist));

    dataset.write(vec, H5::PredType::NATIVE_DOUBLE);

}

void sevnstd::H5out::print_header(std::vector<std::string> &columns) {

    if(!header_printed) {


        file = H5::H5File(filename_h5, H5F_ACC_TRUNC);
        file.createGroup(groupname_h5);

        const size_t xdim = 1;
        const size_t ydim = col_size;

        dim[0] = xdim;
        dim[1] = ydim;
        chunk_dim[0] = xdim;
        chunk_dim[1] = ydim;


        dataspace = H5::DataSpace(2, dim);


        plist = H5::DSetCreatPropList();
        plist.setChunk(2, chunk_dim);
        plist.setDeflate(6);


//local type used to print the header line
        typedef struct HDR {
            char name[20];
        } HDR;
        HDR _hdr[ydim];


        hid_t string_type = H5Tcopy(H5T_C_S1);
        H5Tset_size(string_type, 20);

        dataset = H5::DataSet(file.createDataSet(header_name_h5, string_type, dataspace, plist));

        //HDF5 does not handle c++ vector<string>.... convert to contiguous char!
        for (size_t i = 0; i < ydim; i++) {
            auto it_star = Property::PrintMap.find(columns[i]);
            auto it_binary = Property::PrintMap.find(columns[i]);

            std::string value = columns[i];

            //Does the column you want to print exists in the Property or BinaryProperty map?
            if ( (it_star != Property::PrintMap.end()) || (it_binary != BinaryProperty::PrintMap.end()) ) {
                value.copy(&_hdr[i].name[0], value.size() + 1);
                _hdr[i].name[value.size()] = '\0'; //end of char*
            } else
                svlog.critical("I cannot recognize the parameter you want to print in the hdf5 header", __FILE__, __LINE__);
        }

        // Write data to dataset.
        dataset.write(_hdr, string_type);

        header_printed = true;
    }

}
