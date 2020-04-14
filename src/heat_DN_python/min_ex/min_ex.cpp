#include "mpi.h"
#include "Python.h"
#include "iostream"

// compile via ""mpic++ -I /usr/include/python3.6m/ min_ex.cpp -lpython3.6m"", or similar, depending what python version you are using
int main(int argc, char *argv[]){
    MPI_Init(&argc, &argv);
    int ID_SELF;
    MPI_Comm_rank(MPI_COMM_WORLD, &ID_SELF);
    
    Py_Initialize(); // Python initialize
    
    PyRun_SimpleString("import sys; import os; sys.path.append(os.getcwd())"); // required to find python file in same folder
    PyObject *pFile = PyImport_Import(PyUnicode_FromString("py_min")); // load/execute python file
    
    std::cout << "all done " << ID_SELF << std::endl;
    Py_Finalize(); // Python finalize
    std::cout << "py finalized " << ID_SELF << std::endl;
    //MPI_Finalize();
    return 0;
}
