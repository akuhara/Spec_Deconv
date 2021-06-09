# Spec_Deconv

Receiver function estimation by spectral division with freqeuncy-dependent damping

## Requirements

* Fortran compiler

* FFTW library

## Install

Type `make` in the root directory. You may need to edit `Makefile` in accordance with the compiler you use and the path to the FFTW library.

## Synopsis

`bin/spec_devonv Input_SAC_file_list (option flag)`

* Execution without any arguments displays a brief usage. 

### Format of Input SAC file list

* A Filename of Z component records for each line

* The Z component filename must be appended by ".z".

* R and Z component filenames must be appended by ".r" and ".t", respectively (if aaa.z appears in the list, then both aaa.r and aaa.t must exist).


### Option flags
 

|Flag | function | default value |
|:---:|:---:|:---:|
| a | parameter for Gaussian low-pass filter|  1.0 |
| t | Timing of P onset in sec. (relative to SAC header B)| 200.0 |
| l | Length of signal in sec. | 50.0 |
| p | Fraction of taper (0 -- 0.5) | 0.05 | 
| o | Output directory | ./ |
| n | Lengh of negative (acausal) time for output in sec. | 3.0 |  
| s | Output stacked RF (T) or not (F) | F | 
| f | Output frequency-domain RF (T) or not (F) | F | 
| sp | Calculate S-receiver function (T) or not (F) | F |


### Examples
* `bin/spec_deconv z_file_list a=2.5 t=200 l=50 p=0.05`


### Outputs 

Note that all ouputs are SAC file format.

* *.r.rft: Radial component RF (R deconvolved by Z) 

* *.t.rft: Transverse component RF (T deconvolved by Z)

* *.z.rft: Vertical component RF (Z deconvolved by Z)

Belows are optional outputs

* <Input_SAC_file_list>.r.rft:  Stacked R component RF

* <Input_SAC_file_list>.t.rft:  Stacked T component RF

* <Input_SAC_file_list>.z.rft:  Stacked Z component RF

* *.r.rff: Radial component RF in frequency domain

* *.t.rff: Transverse component RF in frequency domain

* *.z.rff: Vertical component RF in frequency domain


## Notes

* The amount of damping factor is determiend from Z component pre-event noise (for S-receiver function, R component is used instead). For this purpose, input SAC file must contain pre-signal time window that has the same length as the signal time window, which is specified by 'l' flag. 

* Output RFs are normalized such that the signal deconvolved by itself has a unit peak amplitude. However, you may not obtain the exact unit amplitude because of additive damping factor. 

* For S-receiver function mode, output time series undergoes time and polarity reversal. Also, output filenames has a different manner: *.r.rft -> R deconvolved by R; *.t.rft -> T deconvolved by R; *.z.rft => Z devoncolbed by R. 

