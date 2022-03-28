//
// Created by iorio on 21/01/20.
//

#ifndef SEVN_BINSTAR_TEST_H
#define SEVN_BINSTAR_TEST_H

#include <static_main.h>
#include <binstar.h>
#include <IO.h>
#include <fstream>
#include <sstream>
#include <lookup_and_phases.h>
#include <utility_test.h>


class Binstar_TEST : public Binstar {
public:
    using Binstar::Binstar;

};




TEST_CASE( "Test Binstar", "[testbinstar]" ) {


    //REMOVE BUF FOR THE TEST
    std::streambuf *old = std::cout.rdbuf(); // <-- save the old buff.
    std::stringstream ss;
    std::cout.rdbuf(ss.rdbuf());


    std::vector<std::vector<double>> star_prop1 {
            {10.2, 0.02, 0.0},
            {25.0, 0.02, 1.0},
            {35.0,  0.02, 5.0}
    };

    std::vector<std::vector<double>> star_prop2{
        {3.7, 0.02, 0.0},
        {5.0, 0.02, 2.0},
        {25.0,  0.02, 3.0}
    };

    std::vector<std::vector<double>> binary_prop {
            {300.0, 0.0},
            {200.0, 0.2},
            {150, 0.7}
    };

    std::string sn_type{"delayed"}, tini{"zams"}, tf{"end"}, dtout{"all"}, tstep{"0.01"};

    SECTION( "Binstar load new" ){

        int row = star_prop1.size();
        std::vector<std::string> col(Lookup::inputmapbin_nparam.at(InputBinaryOption::_new));
        std::vector<std::vector<std::string>> bstar_prop_new(row, col);


        for (std::vector<int>::size_type i = 0; i != bstar_prop_new.size(); i++){
            col.clear();
            for (auto sp: star_prop1[i]) col.push_back(std::to_string(sp));
            col.push_back(sn_type);
            col.push_back(tini);
            for (auto sp: star_prop2[i]) col.push_back(std::to_string(sp));
            col.push_back(sn_type);
            col.push_back(tini);
            for (auto bp: binary_prop[i]) col.push_back(std::to_string(bp));
            col.push_back(tf);
            col.push_back(dtout);
            bstar_prop_new[i]=col;
            };


        std::string fname="list_bin_new.tmp";
        std::ofstream list_tmp;
        list_tmp.open(&fname[0],  std::ios::out);
        std::ostringstream st_tmp;
        st_tmp.precision(3);
        fill_tmp_file(list_tmp, bstar_prop_new);
        list_tmp.close();
        char *argv[] = {"test", "-list", &fname[0], "-tables", "../../SEVNtracks_G","-tables_HE","../../SEVNtracks_pureHe", "-ibmode","new", "-wmode", "hurley_wind"};
        int argc = sizeof(argv) / sizeof(char*) - 1;
        INFO("Narg input: "<<argc);

        IO sevnio(argc, argv);
        std::remove(&fname[0]);

        //std::vector<Binstar> binaries;
        //binaries.reserve(sevnio.STARS_MATRIX.size());
        //size_t id=0;
        //binaries.emplace_back(Binstar(&sevnio, sevnio.STARS_MATRIX[0], id));
        for (size_t i = 0; i < sevnio.STARS_MATRIX.size(); i++){

            Binstar_TEST bin(&sevnio, sevnio.STARS_MATRIX[i], i);

            std::vector<double> zams_load=bin.get_zams(), Z_load=bin.get_Z();


            REQUIRE(star_prop1[i][0]==zams_load[0]);
            REQUIRE(star_prop2[i][0]==zams_load[1]);
            REQUIRE(star_prop1[i][1]==Z_load[0]);
            REQUIRE(star_prop2[i][1]==Z_load[1]);


        }
    };

    SECTION( "Binstar load legacy" ){

        int row = star_prop1.size();
        std::vector<std::string> col(Lookup::inputmapbin_nparam.at(InputBinaryOption::_new));
        std::vector<std::vector<std::string>> bstar_prop_new(row, col);


        for (std::vector<int>::size_type i = 0; i != bstar_prop_new.size(); i++){
            col.clear();
            col.push_back(std::to_string(star_prop1[i][0]));
            col.push_back(std::to_string(star_prop2[i][0]));
            col.push_back(std::to_string(star_prop1[i][1]));
            col.push_back(std::to_string(star_prop2[i][1]));
            col.push_back(std::to_string(star_prop1[i][2]));
            col.push_back(std::to_string(star_prop2[i][2]));
            col.push_back(std::to_string(binary_prop[i][0]));
            col.push_back(std::to_string(binary_prop[i][1]));
            col.push_back(tf);
            col.push_back(tini);
            col.push_back(tstep);
            col.push_back(sn_type);
            col.push_back(sn_type);
            col.push_back(dtout);
            bstar_prop_new[i]=col;
        };


        std::string fname="list_bin_new.tmp";
        std::ofstream list_tmp;
        list_tmp.open(&fname[0],  std::ios::out);
        std::ostringstream st_tmp;
        st_tmp.precision(3);
        fill_tmp_file(list_tmp, bstar_prop_new);
        list_tmp.close();

        char *argv[] = {"test", "-list", &fname[0], "-tables", "../../SEVNtracks_G","-tables_HE","../../SEVNtracks_pureHe", "-ibmode","legacy"};
        int argc = sizeof(argv) / sizeof(char*) - 1;
        INFO("Narg input: "<<argc);

        IO sevnio(argc, argv);
        std::remove(&fname[0]);

        //std::vector<Binstar> binaries;
        //binaries.reserve(sevnio.STARS_MATRIX.size());
        //size_t id=0;
        //binaries.emplace_back(Binstar(&sevnio, sevnio.STARS_MATRIX[0], id));
        for (size_t i = 0; i < sevnio.STARS_MATRIX.size(); i++){

            Binstar_TEST bin(&sevnio, sevnio.STARS_MATRIX[i], i);

            std::vector<double> zams_load=bin.get_zams(), Z_load=bin.get_Z();


            REQUIRE(star_prop1[i][0]==zams_load[0]);
            REQUIRE(star_prop2[i][0]==zams_load[1]);
            REQUIRE(star_prop1[i][1]==Z_load[0]);
            REQUIRE(star_prop2[i][1]==Z_load[1]);


        }
    };





    std::cout.rdbuf (old);


}

#endif //SEVN_BINSTAR_TEST_H
