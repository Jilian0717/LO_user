### 2-layer DO budget
#### get_DO_bgc_air_sea_shallow_deep.py
* Calculate integrated DO production & consumptions at shallow and deep layers
* The dividing depth between the shallow and the deep layer was given in line 176, in this case, it is -20 m

#### 2-layer DO.png 
* an example for the Salish Sea domain. 
* The exchange flow DO budget term was calculated using the hourly extracted netcdf files in TEF2 after running extrac_section.py. Here, the Exchange flow DO budget at a cross section == oxygen * vel * area. I then separate it into the shallow layer and the deep layer. Therefore, the exchange flow DO budget term caluclated here is in spatial coordinate, not in salinity coordinate. 
* The vertical exchange term is a residual == storage - all other DO budget terms that can be calculated from history files. Ideally, the vertical exchange flow term for the shallow layer and for the deep layer should be very close (their absolute value)
* For simplicity, I put the river term into the shallow layer. It includes both river and wwtps. More strictly, the wwtps should be added to the deep layer
* DO consumption in the shallow layer == nitrification + remineralization; while in the deep layer ==  nitrification + remineralization + SOD.

#### LO_NPZD.pdf
* NPZD structure in fennel.h with customization in LO


