energy detection
=====

### debug
ensure the j is sqrt(-1); otherwise you may use j in the iter
simplify the data model to find bugs. (logic backtracing, 逻辑回溯) ;[72]

### random signal detection
* threshold prefix
./ed_thres.m => Pfa = 0.05, thres = 1.16

### ofdm system information, Dist/ed_ofdm

* ofdm signal 
./gen_ofdm_sys_v01.m => signal_ofdm.mat

* threshold prefix
./ed_thres_ofdm.m => threshold.mat   
use 'empty signal + noise' for training [72]

* Pd and Pfa (rand signal, ofdm signal, practical(training) threshold; [71]
