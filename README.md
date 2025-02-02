This is the release of an updated version of the urban ecohydrological model Urban Tethys-Chloris including a Building Energy Model (UT&C-BEM). Urban Tethys-Chloris is a mechanistic urban ecohydrological model combining principles of ecohydrological land surface modelling with urban canopy modelling. UT&C is a fully coupled energy and water balance model that calculates 2 m air temperature, 2 m humidity, and surface temperatures. It explicitly resolves biophysical and ecophysiological characteristics of ground vegetation, urban trees, and green roofs and models all urban water fluxes including evapotranspiration, canopy interception, infiltration, and soil moisture transport. UT&C accounts for variations in urban densities, building properties, and urban irrigation schemes. The model is able to account for the effects of different plant types on urban climate and hydrology, as well as the effects of the urban environment on plant well-being and performance. 

 
Publications describing model development:

Inclusion of a building energy model in UT&C:

Meili, N., Zheng, X., Takane, Y., Nakajima, K., Yamaguchi, K., Chi, D., Zhu, Y., Wang, J., Qiu, Y., Paschalis, A., Manoli, G., Burlando, P., Tan, P.Y., Fatichi, S.: Modeling the effect of trees on energy demand for indoor cooling and dehumidification across cities and climates, Journal of Advances in Modeling Earth Systems, 2025.

UT&C v1.0 model description article including a detailed technical reference mannual in the supplementary material:

Meili, N., Manoli, G., Burlando, P., Bou-Zeid, E., Chow, W. T. L., Coutts, A. M., Daly, E., Nice, K. A., Roth, M., Tapper, N. J., Velasco, E., Vivoni, E. R., and Fatichi, S.: An urban ecohydrological model to quantify the effect of vegetation on urban climate and hydrology (UT&C v1.0), Geosci. Model Dev., 13, 335–362, https://doi.org/10.5194/gmd-13-335-2020, 2020.

Inclusion of outdoor thermal comfort calculation in UT&C:

Meili, N., Acero, J.A., Peleg, N., Manoli, G., Burlando, P., Fatichi, S. 2021, Vegetation cover and plant-trait effects on outdoor thermal comfort in a tropical city, Building and Environment, 195, 2021, https://doi.org/10.1016/j.buildenv.2021.107733.

 
Model usage tips:

The code on this site includes all scripts and some exemplary forcing data to run the UT&C-BEM model. Start the main computation by running the "ComputationCityName.m" file. Urban geometry and forcing data specifics are specified in the files "Data_UEHM_site_xxx.m" and "UEHMForcingData_xxx.m" stored in the folder "+data_functions". The forcing data (TMY in this case) is stored in the same folder "+data_functions". The most essential simulation outputs are described in the word document "EssentialOutputVariables.doc".

A video tutorial on how to run UT&C (without building energy model) can be found on youtube @UrbanTethysChloris: [How to run UT&C](https://www.youtube.com/watch?v=1dqxOrhI1Dk)

Similarly a video on the radiation partitioning for the meteorological input can also be found on the same channel: [Radiation Partitioning](https://www.youtube.com/watch?v=VOj_lmjdU0s)

 
