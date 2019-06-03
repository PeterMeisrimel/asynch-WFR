/*
Authors: Peter Meisrimel
December 2018
*/

#ifndef WAVEFORM_CPP_
#define WAVEFORM_CPP_

#include "waveform.h"

Waveform::Waveform(int n, int length, double * t, double * data){
	_n      = n;
    _length = length;
	_times  = t;
	_data   = data;
}

Waveform::Waveform(int n, int length, double * t, double * data, double * initvec){
	_n      = n;
    _length = length;
	_times  = t;
	_data   = data;
    for(int i = 0; i < _n; i++)
        set(i, initvec);
}

void Waveform::get_time_full(double * out){
    for(int i = 0; i < _n; i++)
        out[i] = _times[i];
}

void Waveform::set_time_full(double * in){
    for(int i = 0; i < _n; i++)
    _times[i] = in[i];
}

void Waveform::time_shift(double shift){
    for(int i = 0; i < _n; i++)
        _times[i] += shift;
}

// linear interpolation for values in the times vector, linear extrapolation for times outside
void Waveform::eval(double t, double * out){
	int idx;
	// find smallest idx such that _times[idx] > t (yes, strict)
	for(idx = 1; idx < _n; idx++)
		if(_times[idx] > t)
			break;
    if(idx == _n)
        idx--;
    idx--;
    double t_fac = (t - _times[idx])/(_times[idx+1] - _times[idx]);
    for(int i = 0; i < _length; i++)
        out[i] = (1 - t_fac) * _data[idx * _length + i] + t_fac * _data[(idx + 1) * _length + i];
}

void Waveform::set(int idx, double * vec){
    for(int i = 0; i < _length; i++)
        _data[idx*_length + i] = vec[i];
}

void Waveform::set_last(double * vec){
    set(_n - 1, vec);
}

void Waveform::get_last(double * out){
    for(int i = 0; i < _length; i++)
        out[i] = _data[(_n - 1)*_length + i];
}

void Waveform::init_by_last(){
    for(int i = 0; i < _n - 1; i++)
        set(i, _data + _length*(_n - 1));
}

double Waveform::get_err_norm_sq_last(double * in){
    double val, res = 0;
    for(int i = 0; i < _length; i++){
        val = _data[_length*(_n - 1) + i] - in[i];
        res += val * val;
    }
    return res;
}

double Waveform::get_norm_sq_last(){
    double val, res = 0;
    for(int i = 0; i < _length; i++){
        val = _data[_length*(_n - 1) + i];
        res += val * val;
    }
    return res;
}

void WF_swap_data_pointers(Waveform * wf1, Waveform * wf2){
    double * tmp;
    tmp = wf2 -> get_data_p();
    wf2 -> set_data_p(wf1 -> get_data_p());
    wf1 -> set_data_p(tmp);
}

#endif // WAVEFORM_CPP_
