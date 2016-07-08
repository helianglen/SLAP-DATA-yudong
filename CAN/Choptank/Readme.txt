
7/8/2016

This directory contains in-situ water data for the Oct. 15 De-Mar-Var flight, to be used as
a calibration target. 

near_surf_temp_cond_sal.txt: extracted from "SALP-Choptank Data Oct 15 2015.xlsx". See below: 


Fromm: "Kim, Edward J. (GSFC-6170)" <edward.j.kim@nasa.gov>
Date: Thursday, June 16, 2016 at 6:14 PM
To: "Tian, Yudong (GSFC-617.0)[UNIV OF MARYLAND]" <yudong.tian-1@nasa.gov>
Subject: FW: Choptank Data: SLAP Water Calibration, Oct 15, 2015

Hi Yudong,
 
Here’s the water cal data for the DelMarVa flight in Oct 2015, just before going to Canada.
The file to use is “SLAP-Choptank Data Oct 15 2015.xlsx”
 
Paul tested some prototype gadgets in addition to the ‘real’ YSI probe.  I shaded the prototype data dark grey.  Just ignore it.
 
The best data for SLAP TB calibration is the ‘YSI near surface’ data, shaded green.  Using the average should be fine; not much variability.
The yellow data is from 50cm below the water surface; I wouldn’t use it for L-band calibration, but the difference with depth is negligible.
Blue data is metadata on the wind strength & wave height in case a roughness correction is needed.  Let’s try without first.
So, let’s use the green data average.
 
The probe measured temperature, salinity, and conductivity; all should be considered calibrated.
Salinity units are ppt.  So the water Paul & Gabrielle measured was ~half as salty as open ocean water, which might be important
to factor into the TB.
 
YSI also measured water electrical conductivity (due to salinity).  That’s what ‘cond’ means.
Units are micro-Siemens per cm (uS/cm).
 
The equations for dielectric constant of salty water in Ulaby (page 2023-25) appear to use salinity and temp, but not conductivity.
How do you calculate TB?
 
Ed
