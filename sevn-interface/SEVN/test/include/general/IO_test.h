//
// Created by iorio on 16/12/19.
//

#ifndef SEVN_IO_TEST_H
#define SEVN_IO_TEST_H

#include <IO.h>
#include <catch.hpp>
#include <string>
#include <sstream>
#include <fstream>
#include <cstdio>
#include <utility_test.h>

class IO_TEST : public IO {
public:
    using IO::IO;

    int load_test(int n, char** val) {
        IO::load(n, val);
        tablesloaded=0; // to be sura that the estimate could be iteratively call in benchmark test.
        return 0;
    }
};




TEST_CASE( "Test IO", "[testio]" ) {

    //REMOVE BUF FOR THE TEST
    std::streambuf *old = std::cout.rdbuf(); // <-- save the old buff.
    std::stringstream ss;
    std::cout.rdbuf(ss.rdbuf());

    std::string fname="list.tmp";
    std::ofstream list_tmp;

    char *argv[] = {"test", "-list", &fname[0], "-tables", "../../SEVNtracks_G","-tables_HE","../../SEVNtracks_pureHe"};
    int argc = sizeof(argv) / sizeof(char*) - 1;
    INFO("Narg input: "<<argc);

    IO_TEST io_test;


    list_tmp.open(&fname[0],  std::ios::out);

    std::vector<std::vector<double>> star_prop {
            {3.2, 0.02, 0.0},
            {25.0, 0.02, 0.0},
            {35.0,  0.02, 0.5}
    };

    int row=star_prop.size(), col=star_prop[0].size()+4;
    std::vector<std::vector<std::string>> M_star (row, std::vector<std::string>(col));
    std::ostringstream st_tmp;
    st_tmp.precision(3);

    for (std::size_t i=0; i<M_star.size(); ++i){
        std::size_t j;
        for (j=0; j<star_prop[0].size(); ++j){
            st_tmp << std::fixed << star_prop[i][j];
            M_star[i][j] = st_tmp.str();
            st_tmp.str(""); //To clear buffer
        }
        M_star[i][j] = "delayed";
        M_star[i][j+1] = "zams";
        M_star[i][j+2] = "end";
        M_star[i][j+3] = "all";

    }

    fill_tmp_file(list_tmp, M_star);


    SECTION( "IO LOAD STARS" ) {

            io_test.load_test(argc, argv);
            auto & smatrix = io_test.STARS_MATRIX;

            //CHECK DIMENSION
            REQUIRE(smatrix.size()==M_star.size());
            REQUIRE(smatrix[0].size()==M_star[0].size());

            //CHECK VALUES
            for (std::size_t i=0; i<star_prop.size(); ++i){
                std::size_t j;
                for (j=0; j<star_prop[0].size(); ++j){
                    REQUIRE(std::stod(smatrix[i][j])==star_prop[i][j]);
                }
                REQUIRE(smatrix[i][j]=="delayed");
                REQUIRE(smatrix[i][j+1]=="zams");
                REQUIRE(smatrix[i][j+2]=="end");
                REQUIRE(smatrix[i][j+3]=="all");
            }

    };


    std::cout.rdbuf (old);
    std::remove(&fname[0]);

}




#endif //SEVN_IO_TEST_H
