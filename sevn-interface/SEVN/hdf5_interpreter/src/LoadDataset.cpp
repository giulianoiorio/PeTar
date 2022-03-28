//
// Created by mario on 03/02/19.
//

#include <LoadSStars.h>
#include <LoadDataset.h>
#include <global.h>
//LoadDataset _loaddataset;

void LoadDataset::Execute() {


    try {

        //loadeddata = false;
        //loadedheader = false;
        //if the file is not loaded, please load it!!
        if(!_load.isloaded()) {
            const std::string err = "Please load the h5 file first (command: load(\"filename\") )";
            throw err;
        }

        H5::DataSet dataset = _load.get_h5file()->openDataSet(args[0] + "/" + args[1]);

        H5::DataSpace dataspace = dataset.getSpace();
        auto dataClass = dataset.getTypeClass();



        size_t byteSize;
        bool header = false;
        bool data = false;

        std::vector<double> data_out_data;


        if (dataClass == H5T_STRING) {
            //This means that I am reading the header file
            header = true;
            auto floatType = dataset.getStrType();
            byteSize = floatType.getSize();
        }

        else if (dataClass == H5T_FLOAT) {
            //This means that I am reading the data file
            data = true;
        }



        //get the dimensions
        int rank = dataspace.getSimpleExtentNdims();

        int ndims = dataspace.getSimpleExtentDims(dims_out, nullptr);
      //  std::cout << "rank " << rank << ", dimensions " <<
               //   (unsigned long) (dims_out[0]) << " x " <<
              //    (unsigned long) (dims_out[1]) << std::endl;

        const hsize_t DIM0 = dims_out[0];
        const hsize_t DIM1 = dims_out[1];





        int        numfilt;
        size_t     nelmts={1}, namelen={1};
        unsigned  flags, filter_info, cd_values[1];
        size_t idx;
        char       name[1];
        H5Z_filter_t filter_type;
        // Get the create property list of the dataset.
        H5::DSetCreatPropList plist(dataset.getCreatePlist());

        // Get the number of filters associated with the dataset.
        numfilt = plist.getNfilters();
       // std::cout << "Number of filters associated with dataset: " << numfilt << std::endl;

        for (idx=0; idx < numfilt; idx++) {
            nelmts = 0;

            filter_type = plist.getFilter(idx, flags, nelmts, cd_values, namelen, name , filter_info);

           // std::cout << "Filter Type: ";

#if 0
            switch (filter_type) {
                case H5Z_FILTER_DEFLATE:
                    //std::cout << "H5Z_FILTER_DEFLATE" << std::endl;
                    break;
                case H5Z_FILTER_SZIP:
                    //std::cout << "H5Z_FILTER_SZIP" << std::endl;
                    break;
                default:
                   // std::cout << "Other filter type included." << std::endl;
            }

#endif
        }



        if(data) {
            data_out_data.resize(DIM0*DIM1);
            dataset.read(&data_out_data[0], H5::PredType::NATIVE_DOUBLE);




            alldata.resize(DIM0*DIM1);

            for(int i = 0; i < DIM0; i++) {
                for (int j = 0; j < DIM1; j++) {
                    alldata[j + DIM1*i] = (data_out_data[j + DIM1*i]);
                }
            }

            loadeddata = true;
            dataset_name = args[0] + "/" + args[1];
        }
        else if (header) {
            //In case I am reading the header file

            hid_t string_type = H5Tcopy( H5T_C_S1 );
            H5Tset_size( string_type, byteSize );
            typedef char _hdrfile [byteSize];
            _hdrfile data_out_header[DIM1];
            dataset.read(data_out_header, string_type);


            allheader.resize(DIM1);

            for (int j = 0; j < DIM1; j++) {
                allheader[j]=data_out_header[j];
            }

            loadedheader = true;
        }
        else
            std::cout<<" Cannot recognize the type of read dataset";





    }
    catch( H5::FileIException &error )
    {
        error.printErrorStack();

    }
        // catch failure caused by the DataSet operations
    catch( H5::DataSetIException &error )
    {
        error.printErrorStack();

    }
        // catch failure caused by the DataSpace operations
    catch( H5::DataSpaceIException &error )
    {
        error.printErrorStack();

    }

}
