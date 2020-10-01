/*
Authors: Peter Meisrimel
December 2018
*/

#ifndef WAVEFORM_H_
#define WAVEFORM_H_

/*
A general purpose class for organization and handling of discrete Waveforms.
Important note: This class does not allocate the memory for the data and is done inside the WFR methods. The rationale behind this is the need to allocate memory using MPI_win_allocate for One-sided communication.
*/
class Waveform{
protected:
    int _n;      // number of points in time
    int _length; // size of vector
    double * _times;
    double * _data;
public:
    Waveform(int, int, double *, double *);           // n, length, times, pointer to data
    Waveform(int, int, double *, double *, double *); // same plus vector for initialization

    void get_time_full      (double *);         // times_vec (out)
    void set_time_full      (double *);         // times_vec
    void time_shift         (double);           // time
    // evaluate waveform at a given time t, uses linear interpolation or extrapolation
    virtual void eval               (double, double *); // t, vec (out)
    virtual void set                (int, double *);    // idx, vec
    virtual void set_last           (double *);         // vec
    virtual void get                (int, double *);    // idx, vec
    virtual void get_last           (double *);         // vec (out)
    virtual void get_all            (Waveform *);       // copy data vector
    virtual void init_by_last       ();                 // set(get_last) for all indices
    virtual void relax_by           (double, Waveform *); // w, WF_old
    virtual void relax_by (double, double*); //testing
    virtual void relax_by_single    (double, int, double *); // w, idx, WF_old
    virtual void copy(double*); // testing 

// perform relaxation, (w, idx_old, new, out)

    // get pointer to according start of vector
    double * operator[] (int i){ return _data + i*_length;};
    double get_time(int idx){ return _times[idx];};
    // squared 2-norm of difference between last vector and input vector
    double get_err_norm_sq_last (double *); // vec(in)
    double get_norm_sq_last ();
    void print();

    double* get_data_p()          { return _data;};
    void   set_data_p(double * p) { _data = p;};
};

void WF_swap_data_pointers(Waveform *, Waveform *);
#endif // WAVEFORM_H_
