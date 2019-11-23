# AQUILA

This is a InGaAs/GaAs version of  AQUILA based on the work of Dr. Martin Rother  
You are also find his work via http://www.rotherland.de/en/aquila.html

The redistributed version include strain effect to model symmetric InGaAs/GaAs quantum well. 
So strain effects are not included in GaAs/AlGaAs system. In addition, assymetric quantum well
could generate wrong result due to our symmetric code.

If you want to use  InGaAs/GaAs  system,
1. use  buildstructure_InGaAs instead of buildstructure
2. specify barrier mole fraction number, such as  aquila_material.barrier=0.8
