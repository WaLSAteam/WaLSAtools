; -----------------------------------------------------------------------------------------------------
; WaLSAtools: Wave analysis tools
; Copyright (C) 2025 WaLSA Team - Shahin Jafarzadeh et al.
;
; Licensed under the Apache License, Version 2.0 (the "License");
; you may not use this file except in compliance with the License.
; You may obtain a copy of the License at
; 
; http://www.apache.org/licenses/LICENSE-2.0
; 
; Unless required by applicable law or agreed to in writing, software
; distributed under the License is distributed on an "AS IS" BASIS,
; WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
; See the License for the specific language governing permissions and
; limitations under the License.
; 
; Note: If you use WaLSAtools for research, please consider citing:
; Jafarzadeh, S., Jess, D. B., Stangalini, M. et al. 2025, Nature Reviews Methods Primers, in press.
; -----------------------------------------------------------------------------------------------------


; $Id: //depot/Release/ENVI53_IDL85/idl/idldir/lib/hilbert.pro#1 $
;
;+
; NAME:
;	HILBERT --> walsa_hilbert (slight modification for double precesion)
;
; PURPOSE:
;	Return a series that has all periodic terms shifted by 90 degrees.
;
; CATEGORY:
;	G2 - Correlation and regression analysis
;	A1 - Real arithmetic, number theory.
;
; CALLING SEQUENCE:
;	Result = HILBERT(X [, D])
;
; INPUT:
;	X:	A floating- or complex-valued vector containing any number
;		of elements.
;
; OPTIONAL INPUT:
;	D:	A flag for rotation direction.  Set D to +1 for a
;		positive rotation.  Set D to -1 for a negative rotation.
;		If D is not provided, a positive rotation results.
;
; OUTPUTS:
;	Returns the Hilbert transform of the data vector, X.  The output is
;	a complex-valued vector with the same size as the input vector.
;
; COMMON BLOCKS:
;	None.
;
; SIDE EFFECTS:
;	HILBERT uses FFT() so this procedure exhibits the same side
;	effects with respect to input arguments as that function.
;
; PROCEDURE:
;	A Hilbert transform is a series	that has had all periodic components
;	phase-shifted by 90 degrees.  It has the interesting property that the
;	correlation between a series and its own Hilbert transform is
;	mathematically zero.
;
;	The method consists of generating the fast Fourier transform using
;	the FFT() function and shifting the first half of the transform
;	products by +90 degrees and the second half by -90 degrees.  The
;	constant elements in the transform are not changed.
;
;	Angle shifting is accomplished by multiplying or dividing by the
;	complex number, I=(0.0000, 1.0000).  The shifted vector is then
;	submitted to FFT() for transformation back to the "time" domain and the
;	output is divided by the number elements in the vector to correct for
;	multiplication effect peculiar to the FFT algorithm.
;
; REVISION HISTORY:
;	JUNE, 1985,	Written, Leonard Kramer, IPST (U. of Maryland) on site
;			contractor to NASA(Goddard Sp. Flgt. Cntr.)
;-
FUNCTION walsa_hilbert,X,D   ; performs the Hilbert transform of some data.
	ON_ERROR,2           ; Return to caller if an error occurs
	Y=FFT(X,-1,/DOUBLE)     ; go to freq. domain.
	N=N_ELEMENTS(Y)
	I=COMPLEX(0.0,1.0)
	IF N_PARAMS(X) EQ 2 THEN I=I*D
	N2=N/2-1	     ; effect of odd and even # of elements
			     ; considered here.
	Y[1]=Y[1:N2]*I       ; multiplying by I rotates counter c.w. 90 deg.
	N2=N-N2
	Y[N2]=Y[N2:N-1]/I
	Y=FFT(Y,1,/DOUBLE)  ; go back to time domain
	RETURN,Y
END