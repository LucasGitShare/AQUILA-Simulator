# AQUILA
Copyright (c) 1999, Martin Rother, 2019,Lucas Yang

All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1) Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer. 
2) Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution. 
3) Neither the name of the <ORGANIZATION> nor the names of its contributors
   may be used to endorse or promote products derived from this software
   without specific prior written permission. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
----------------------------------------------------------------------------------

This is a InGaAs/GaAs version of  AQUILA based on the work of Dr. Martin Rother  
Thank Dr. Martin Rother for his permission and feedback.
You are also find his work via http://www.rotherland.de/en/aquila.html

The redistributed version include strain effect to simulate symmetric InGaAs/GaAs quantum well. 
So strain effects are not included in GaAs/AlGaAs system. In addition, assymetric quantum well
could generate wrong result due to our symmetric code.

If you want to use  InGaAs/GaAs  system,
1. use  buildstructure_InGaAs instead of buildstructure
2. specify barrier mole fraction number, such as  aquila_material.barrier=0.8
