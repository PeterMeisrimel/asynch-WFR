/*
Authors: Peter Meisrimel
April 2019
*/

#ifndef WAVEFORM_LOCKING_H_
#define WAVEFORM_LOCKING_H_

#include "waveform.h"
#include "mpi.h"

/*
Special instance of the Waveform class to include locking around the associated window object for evaluations
*/
class Waveform_locking : public Waveform{
private:
    MPI_Win * _win;
    int ID_SELF;
public:
    // same as basic ones plus pointer to window
    Waveform_locking(int, int, double *, double *, MPI_Win *, int);
    Waveform_locking(int, int, double *, double *, double *, MPI_Win *, int);
    // overwrite all functions that have data accesses, as they now require MPI_Win_locks
	void eval               (double, double *); // t, vec (out)
	void set                (int, double *);    // idx, vec
	void set_last           (double *);         // vec
    void get_last           (double *);         // vec (out)
    void init_by_last       ();                 // set(get_last) for all indices

    // squared 2-norm of difference between last vector and input vector
    double get_err_norm_sq_last (double *); // vec(in)
};

void WF_swap_data_pointers(Waveform *, Waveform *);
#endif // WAVEFORM_LOCKING_H_
