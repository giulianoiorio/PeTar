// Don't forget the includes
#include <iostream>
#include <Python.h>
#include "plot_matplotlib.h"
#include <vector>

using namespace std;

int main(){

// ...

Py_Initialize(); // Initialize Python
std::vector<double> X;
std::vector<double> Y;

for(int i = 0; i < 100; i++){
	X.push_back(i);
	Y.push_back(i*i);
}


// Plot using matplotlib
plot_matplotlib * plot = new plot_matplotlib();
plot->add_somedata(&X,&Y);
plot->set_range_auto();
plot->set_xlabel("Time");
plot->set_ylabel("Error");
plot->show();

}
