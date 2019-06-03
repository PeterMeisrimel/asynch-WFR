/*
Authors: Peter Meisrimel
April 2019
*/

#ifndef WAVEFORM_LOCKING_CPP_
#define WAVEFORM_LOCKING_CPP_

#include "waveform_locking.h"

Waveform_locking::Waveform_locking(int n, int length, double * t, double * data, MPI_Win * win_data, int ID): Waveform(n, length, t, data){
    _win = win_data;
    ID_SELF = ID;
}

Waveform_locking::Waveform_locking(int n, int length, double * t, double * data, double * initvec, MPI_Win * win_data, int ID): Waveform(n, length, t, data){
    _win = win_data;
    ID_SELF = ID;
    for(int i = 0; i < _n; i++)
        // locks in set function
        set(i, initvec);
}

// linear interpolation for values in the times vector, linear extrapolation for times outside
void Waveform_locking::eval(double t, double * out){
	int idx;
	// find smallest idx such that _times[idx] > t (yes, strict)
	for(idx = 1; idx < _n; idx++)
		if(_times[idx] > t)
			break;
    if(idx == _n)
        idx--;
    idx--;
    double t_fac = (t - _times[idx])/(_times[idx+1] - _times[idx]);
    MPI_Win_lock(MPI_LOCK_SHARED, ID_SELF, 0, *_win);
    for(int i = 0; i < _length; i++)
        out[i] = (1 - t_fac) * _data[idx * _length + i] + t_fac * _data[(idx + 1) * _length + i];
    MPI_Win_unlock(ID_SELF, *_win);
}

void Waveform_locking::set(int idx, double * vec){
    MPI_Win_lock(MPI_LOCK_EXCLUSIVE, ID_SELF, 0, *_win);
    for(int i = 0; i < _length; i++)
        _data[idx*_length + i] = vec[i];
    MPI_Win_unlock(ID_SELF, *_win);
}

void Waveform_locking::set_last(double * vec){
    set(_n - 1, vec); // can be removed?
}

void Waveform_locking::get_last(double * out){
    MPI_Win_lock(MPI_LOCK_SHARED, ID_SELF, 0, *_win);
    for(int i = 0; i < _length; i++)
        out[i] = _data[(_n - 1)*_length + i];
    MPI_Win_unlock(ID_SELF, *_win);
}

void Waveform_locking::init_by_last(){
    // can be removed?
    for(int i = 0; i < _n - 1; i++)
        set(i, _data + _length*(_n - 1));
}

double Waveform_locking::get_err_norm_sq_last(double * in){
    double val, res = 0;
    MPI_Win_lock(MPI_LOCK_SHARED, ID_SELF, 0, *_win);
    for(int i = 0; i < _length; i++){
        val = _data[_length*(_n - 1) + i] - in[i];
        res += val * val;
    }
    MPI_Win_unlock(ID_SELF, *_win);
    return res;
}

#endif // WAVEFORM_LOCKING_CPP_
