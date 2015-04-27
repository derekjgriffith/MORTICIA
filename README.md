# MORTICIA
MORTICIA is Monte-carlo Optical Rendering for Theatre Investigations of Capability under the Influence of the Atmosphere, 
an open-source optical surveillance modelling framework written mainly in Python.

  The theatre of operations is a geographical region defined by a polygon in which an optical surveillance system or 
systems are deployed or to be deployed. Besides the engineering characteristics of the optical surveillance system, 
the atmosphere plays a key role in determining the range at which targets can be Detected, Recognised or Identified (DRI).
  
  Some of the key capabilities of MORTICIA are:
  

 + Physically accurate rendering of targets of interest in atmospheric conditions that are statistically representative 
 of the theatre of operations.
 + Statistically valid Monte-carlo simulations of the Detect, Recognition and Identification (DRI) ranges of such targets
 in the theatre of operations.
 
 ## Theatre Climatology
 The climatology within the theatre must be taken into account in order to obtain statistically valid results. Some of the
 most important elements are:
 
 + Cloud climatology
 + Aerosol climatology
 + Surface reflectance climatology (e.g. vegetation state)
 
 In addition, there are certain fixed geometrical aspects relevant to surveillance that are defined by selection of the
 theatre. These include solar/lunar azimuth and elevation statistics and surface topography.
 
 ## Outputs
 The main output of MORTICIA takes the form of a database of results that define target appearance and DRI, having chosen
 a sensor position and a target position within the theatre at a particular time of day and day of year. The atmospheric 
 conditions can be chosen explicitly or generated randomly on the basis of the theatre climatology.
 
 
 
 
 
 

