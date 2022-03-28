//
// Created by mario on 04/02/19.
//

#include <Plot.h>
#include <iostream>
#include <vector>
#include <global.h>
#include <Python.h>
#include <plot_matplotlib.h>

void Plot::Execute() {
    try {

        if (!_loaddataset.hdrloaded()) {
            const std::string err = "Please load the header dataset first (command: load_dataset(stars, header) )";
            throw err;
        }

        if(!_loaddataset.dataloaded()){
            const std::string err = "Please load the dataset first (command: load_dataset(stars, ''datasetname'') )";
            throw err;
        }


        std::string x = args[0];
        std::vector<std::string> y;
        y.resize(args.size()-1);

        for(size_t i = 1; i < args.size(); i++)
            y[i-1] = args[i];

        ypos.resize(y.size());

        //check if x and y are both in the header... and return their position

        std::vector<std::string> *hdr = _loaddataset.get_allheader();
        std::vector<double> *data = _loaddataset.get_alldata();
        size_t dimy = _loaddataset.get_ydim();
        size_t dimx = _loaddataset.get_xdim();

        size_t pos = 0;
        pos = find(hdr->begin(), hdr->end(), x) - hdr->begin();


        if (pos >= hdr->size()) {
            std::string err = "Impossible to get the specified quantitiy [" + x + "]" + "Possible values = ";
	for(int i = 0; i < hdr->size(); i++)
		err += "[" + hdr->at(i) + "]";
	
            throw err;
        }

        xpos = pos;


        for(size_t k = 0; k < y.size(); k++) {
            pos = 0;
            pos = find(hdr->begin(), hdr->end(), y[k]) - hdr->begin();


            if (pos >= hdr->size()) {
                std::string err = "Impossible to get the specified quantitiy [" + y[k] + "]" + "Possible values = ";
		for(int i = 0; i < hdr->size(); i++)
                err += "[" + hdr->at(i) + "]";
		
                throw err;
            }

            ypos[k] = pos;
        }


        Py_Initialize(); // Initialize Python
        std::vector<double> X;
        std::vector<double> Y;

        X.resize(dimx);
        Y.resize(dimx);

        // Plot using matplotlib
        plot_matplotlib * plot = new plot_matplotlib();

        for (size_t i = 0; i < dimx; i++) {
            X[i] = data->at(i * dimy + xpos);
            Y[i] = data->at(i * dimy + ypos[0]);
        }
        plot->add_somedata(&X, &Y, plot->get_marker());

        for(size_t k = 1; k < ypos.size(); k++) {

            for (size_t i = 0; i < dimx; i++) {
                Y[i] = data->at(i * dimy + ypos[k]);
            }

            plot->add_somedata(&Y, plot->get_marker());
        }

        plot->set_range_auto();
        plot->set_xlabel(x);
        if(y.size()==1)
            plot->set_ylabel(y[0]);
        else
            plot->set_ylabel("Various quantities");
        plot->show();

        plot->reset_counters();

    }
    catch( const std::string &msg ) {
        std::cout << msg << "\n";
    }
}
