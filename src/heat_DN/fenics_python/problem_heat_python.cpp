/*
Authors: Peter Meisrimel
September 2018
*/

#ifndef PROBLEM_HEAT_PY_CPP_
#define PROBLEM_HEAT_PY_CPP_

#include "Python.h"
#include "problem_heat_python.h"
#include "iostream"
#include "cassert"
#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL cool_ARRAY_API
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"

Problem_heat_python::Problem_heat_python(int gridsize, double a, double g, const char* class_name){
    other_init_done = false;
    _length = gridsize + 2;
    
//     this appears to change the paths to the current directory? without it, the import fails
//    PyRun_SimpleString("import sys; import os");
//    PyRun_SimpleString("sys.path.append(os.getcwd())");
    PyRun_SimpleString("import sys");
    PyRun_SimpleString("sys.path.append('/home/peter/asynch-WFR/src/heat_DN/fenics_python/')");
//    PyRun_SimpleString("os.environ['OMP_NUM_THREADS'] = '1'");
    
    // import file, resp. python script
//    PyObject *pFile = PyImport_Import(PyUnicode_FromString("heat_py"));
    PyObject *pFile = PyImport_Import(PyUnicode_FromString("heat_fenics_grad"));
    assert(pFile != NULL);
    
    PyObject *pDict = PyModule_GetDict(pFile);
    // get class name to initialize from dictionary
    PyObject *pClass_name = PyDict_GetItemString(pDict, class_name);
    
    // initialize values for constructor class
    PyObject *pClassCall = PyTuple_New(3);
    PyTuple_SetItem(pClassCall, 0, PyLong_FromLong((long)gridsize));
    PyTuple_SetItem(pClassCall, 1, PyFloat_FromDouble(a));
    PyTuple_SetItem(pClassCall, 2, PyFloat_FromDouble(g));
    
    // initialize class object
    pClass_obj = PyObject_CallObject(pClass_name, pClassCall);
}

void Problem_heat_python::create_checkpoint(){
    PyObject_CallMethodObjArgs(pClass_obj, PyUnicode_FromString("create_checkpoint"), NULL);
}

void Problem_heat_python::reset_to_checkpoint(){
    PyObject_CallMethodObjArgs(pClass_obj, PyUnicode_FromString("reset"), NULL);
}

void Problem_heat_python::get_u0(double* u_out){
    PyObject *pOutput = PyObject_CallMethodObjArgs(pClass_obj, PyUnicode_FromString("get_u0"), NULL);
    Py_XINCREF(pOutput);
    PyArrayObject *pOutput_arr = reinterpret_cast<PyArrayObject*>(pOutput);
    Py_XINCREF(pOutput_arr);
    double * u_out_local = reinterpret_cast<double*>(PyArray_DATA(pOutput_arr));
    for(int i = 0; i < _length; i++)
        u_out[i] = u_out_local[i];
    
    Py_XDECREF(pOutput);
    Py_XDECREF(pOutput_arr);
}

void Problem_heat_python_D::do_step(double t, double dt, double *u_out, Waveform *WF_in){
    Pdt = PyFloat_FromDouble(dt);
    
    npy_intp dims[1]{_length_other};
    double *cInput = new double[_length_other];
    WF_in->eval(t + dt, cInput);
    PyObject *pInput = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, reinterpret_cast<void*>(cInput));
    Py_XINCREF(pInput);
    
    PyObject *pOutput = PyObject_CallMethodObjArgs(pClass_obj, PyUnicode_FromString("do_step"), Pdt, pInput, NULL);
    Py_XINCREF(pOutput);
    PyArrayObject *pOutput_arr = reinterpret_cast<PyArrayObject*>(pOutput);
    Py_XINCREF(pOutput_arr);
    double * u_out_local;
    u_out_local = reinterpret_cast<double*>(PyArray_DATA(pOutput_arr));
    for(int i = 0; i < _length; i++)
        u_out[i] = u_out_local[i];
        
    Py_XDECREF(pInput);
    Py_XDECREF(pOutput);
    Py_XDECREF(pOutput_arr);
}

void Problem_heat_python_N::do_step(double t, double dt, double *u_out, Waveform *WF_in){
    Pdt = PyFloat_FromDouble(dt);
    npy_intp dims[1]{_length_other};
    
    double *cInput_flux_old = new double[_length_other];
    WF_in->eval(t, cInput_flux_old);
    PyObject *pInput_old = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, reinterpret_cast<void*>(cInput_flux_old));
    Py_XINCREF(pInput_old);
    
    double *cInput_flux_new = new double[_length_other];
    WF_in->eval(t + dt, cInput_flux_new);
    PyObject *pInput_new = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, reinterpret_cast<void*>(cInput_flux_new));
    Py_XINCREF(pInput_new);
    
    PyObject *pOutput = PyObject_CallMethodObjArgs(pClass_obj, PyUnicode_FromString("do_step"), Pdt, pInput_old, pInput_new, NULL);
    Py_XINCREF(pOutput);
    PyArrayObject *pOutput_arr = reinterpret_cast<PyArrayObject*>(pOutput);
    Py_XINCREF(pOutput_arr);
    double * u_out_local;
    u_out_local = reinterpret_cast<double*>(PyArray_DATA(pOutput_arr));
    for(int i = 0; i < _length; i++)
        u_out[i] = u_out_local[i];
    
    Py_XDECREF(pInput_new);
    Py_XDECREF(pOutput);
    Py_XDECREF(pOutput_arr);
}

#endif //PROBLEM_HEAT_PY_CPP_