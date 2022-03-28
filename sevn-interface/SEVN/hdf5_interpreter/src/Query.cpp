//
// Created by mario on 04/02/19.
//

#include <Query.h>
#include <global.h>


void Header::Execute() {
    std::vector<std::string> * hdr = _loaddataset.get_allheader();

    std::cout<<" Header: [";
    for(size_t i = 0; i < hdr->size(); i++)
        std::cout<<hdr->at(i)<<";  ";
    std::cout<<"]"<<std::endl;

}

void Dataset::Execute() {

    std::cout<<" Dataset: ["<<_loaddataset.get_name()<<"]"<<std::endl;

}
