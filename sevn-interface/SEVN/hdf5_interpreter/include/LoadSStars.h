//
// Created by mario on 03/02/19.
//

#ifndef H5INTER_LOADSSTARS_H
#define H5INTER_LOADSSTARS_H

#include <command.h>
#include <processor.h>
#include <H5Cpp.h>

class LoadSStars : public Command {

public:

    bool fileloaded;

    LoadSStars() {
        fileloaded = false;
        CommandProcessor::Instance().Register( "load", this, 1); //name, command, nargs
    }

    void Execute() {

        try {

            fileloaded = false;

            H5::Exception::dontPrint();

            file.close();
            file.openFile(args[0], H5F_ACC_RDWR);

            group_names.erase(group_names.begin(), group_names.end());
            herr_t idx = H5Literate(file.getId(), H5_INDEX_NAME, H5_ITER_NATIVE, nullptr, file_info, &group_names);

            dataset_names.erase(dataset_names.begin(), dataset_names.end());
            dataset_names.resize(group_names.size());

            for(size_t i = 0; i < group_names.size(); i++) {

                H5::Group group(file.openGroup(group_names[i]));

                idx = H5Literate(group.getId(), H5_INDEX_NAME, H5_ITER_INC, nullptr, file_info, &dataset_names[i]);

                std::cout << " Group: [" << group_names[i] <<"]"<< std::endl;
                std::cout << " Datasets: [";
                for (size_t j = 0; j < dataset_names[i].size(); j++)
                    std::cout << dataset_names[i][j] <<";  ";
                std::cout<<"]"<<std::endl;

            }
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

        fileloaded = true;

    }

    inline const bool isloaded() const { return fileloaded;}
    inline H5::H5File* get_h5file() {return &file;}

//    inline const std::vector<std::string> get_gdatasets(int i) {return dataset_names[i];}
 //   inline const std::vector<std::string> get_groups() {return group_names;}
//    inline H5::H5File* get_h5file() {return file; };

private:
    std::vector<std::vector<std::string>> dataset_names;
    std::vector<std::string> group_names;
    H5::H5File file;

    static herr_t file_info(hid_t loc_id, const char *name, const H5L_info_t *linfo, void *opdata){


        H5O_info_t      infobuf;
        auto names=reinterpret_cast< std::vector<std::string>* >(opdata);

        herr_t status = H5Oget_info_by_name (loc_id, name, &infobuf, H5P_DEFAULT);

        switch (infobuf.type) {
            case H5O_TYPE_GROUP:
               // printf ("  Group: %s\n", name);
                names->push_back(name);
                break;
            case H5O_TYPE_DATASET:
               // printf ("  Dataset: %s\n", name);
                names->push_back(name);
                break;
            case H5O_TYPE_NAMED_DATATYPE:
               // printf ("  Datatype: %s\n", name);
                names->push_back(name);
                break;
            default:
                printf ( "  Unknown: %s\n", name);
              //  names->push_back(name);
        }


        return 0;


    }

};


#endif //H5INTER_LOADSSTARS_H
