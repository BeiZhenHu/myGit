Readme file to explain the different files and parameters:

There are 220 ABC database folders, starting from "ABC_database_7" to "ABC_database_226" the numbering is from 7 to 226 is the number of the ABC boards

Each folder ABC_database_XX contains 4 text files

calib_info.txt : contains Charge Calibration parameters for Preamplifier Gain of 20 for HG and LG calibration : when there is 0 for a parameters in LG ( moslty at the 4 last entries, it means that there is no calibration for this channel in LG, DUE TO SMALL NUMBER OF EVENTS < 200 IN MEASURMENT )


charge_info.txt : contains the value of the maximum charge reached in HG ( normally close to 800 ADCu )
Pedestal_info.txt : contains the Pedestal charge paramters, ( mean and sigma of the pedestal distribution for ping and pong and HG anf LG)
time_info.txt : contains the Fine time information for HG/LG threshold at 520 DACu, the fine_time min value and maximum value from the finetime distribution 


Notes on ABC problems:

ABC number 29, 42, 161, 177, 180 : ASIC 3 ( channel_nb: 48 -> 63): No pedestal information, thus no calibration done ( entries in calib_info.txt  = 0 )

ABC number 41, channel 60, not responding, no calibration, no pedestal info

For the three ABC boards : "ABC_159 ABC_179 ABC_201" : not calibrated, no calib file, need to be recalibrated 



