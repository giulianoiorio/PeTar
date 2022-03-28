//
// Created by Giuliano Iorio on 2020-01-25.
//

#ifndef SEVN_UTILITY_TEST_H
#define SEVN_UTILITY_TEST_H


template<typename T>
void fill_tmp_file(std::ofstream & file, const std::vector<std::vector<T>> & M){

    if (!file.is_open()) file.open("list.tmp",  std::ios::out);

    std::string sep = "    ";
    for (auto row : M) {
        for (auto col: row)
            file << col << sep;
        file<<std::endl;
    }
    file.close();
};


#endif //SEVN_UTILITY_TEST_H
