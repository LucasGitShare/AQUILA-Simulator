function a=ingaasmaterial(x,prop)

%ingaasmaterial material database for GaAs/InGaAs mertial system
%
%a=ingaasmaterial(x,prop)
%
%takes the x-content in field x and returns a material property
%according to the string prop.
%This database is somewhat rudimentary and incomplete and does not always
%give correct results. It should be ok for computations with Gamma electrons.
%X- and L-electrons are probably not described very well.
%Valid property strings are:
%
%all energies in eV, all masses in M0
%
%E_G6G8    energy between Gamma6 and Gamma8 point = normal gap
%E_G6G7	  Gamma6-Gamma7 = split-off
%E_G8X6    Gamma8-X6 = indirect gap for x>.45
%E_G8L6    Gamma8-L6 = indirect gap f x>.45      
%
%M_eG      electron mass at Gamma point = DOS-mass = cond masse = t = l
%M_eXd     electron dos mass at X point
%M_eXl     electron longitudinal X mass
%M_eXt     electron transversal X mass
%M_eX      electron X mass
%M_eLd     electron dos mass at L point
%M_eLl     electron longitudinal L mass
%M_eLt     electron transversal L mass
%M_eL      electron L mass
%
%M_hh001   heavy hole 001 mass
%M_hh110   heavy hole 110 mass
%M_hh111   heavy hole 111 mass
%M_lh001   light hole 001 mass
%M_lh110   light hole 110 mass
%M_lh111   light hole 111 mass
%M_hhd     heavy hole DOS mass
%M_lhd     light hole DOS mass
%M_so      split-off hole mass
%
%Gamma1    Luttinger parameter
%Gamma2    Luttinger parameter
%Gamma3    Luttinger parameter
%
%The following energies give absolute values of the corresponding quantities.
%These energies are normalized such, that the GaAs CB edge at T=0 has energy 0.
%E_c       conduction band edge = min(E_g,E_x,E_l)
%E_v       valence band edge
%E_x       X band edge
%E_l       L band edge
%E_g       Gamma band edge
%E_so      split-off hole band edge
%E_hh      heavy hole band edge
%E_lh      light hole band edge
%E_d       donator level
%Delta_Ev  valence band offset
%
%eps       relative dielectric constant

%Copyright 1999 Martin Rother
%
%This file is part of AQUILA.
%
%AQUILA is free software; you can redistribute it and/or modify
%it under the terms of the BSD License as published by
%the Open Source Initiative according to the License Policy
%on MATLAB(R)CENTRAL.
%
%AQUILA is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%BSD License for more details.

global aquila_control aquila_material

%check for correct calling sequence
if bitget(aquila_control.progress_check,1)==0
   error('material: INITAQUILA must be run before accessing the materials database !')
end

a=zeros(size(x));
p=find(x>0.69);
if ~isempty(p)
   disp('material: x>0.69 not supported');
   return;
end

%some values are taken from tables
%the interpolation is computed only once and then used on subsequent calls.
if ~exist('aquila_material.alpha') %initialize values for interpolation, only done once
   aquila_material.e0=interp1([0 0.18 0.27 0.53 0.69],[1.517 1.771 1.932 2.251 2.485],[0:0.01:0.69]);
   aquila_material.alpha=interp1([0 0.18 0.27 0.53 0.69],[5.5 6.3 6.58 7.04 7.88]*1e-4,[0:0.01:0.69]);
   aquila_material.beta=interp1([0 0.18 0.27 0.53 0.69],[0.225 0.236 0.248 0.261 0.302],[0:0.01:0.69]);
end
T=aquila_control.T;

%return the requested material property
switch prop      
case {'E_G6G8'} % Gamma6-Gamma8 = normal gap
   index=find(x>0);
   x2=x(index);
   s=size(x2);
   x2=x2(:);
   e0=interp1([0:0.01:0.69],aquila_material.e0,x2,'nearest');
   alpha=interp1([0:0.01:0.69],aquila_material.alpha,x2,'nearest');
   beta=interp1([0:0.01:0.69],aquila_material.beta,x2,'nearest');
   a2= e0-T*T*alpha./(beta+T)  + ingaasmaterial(x2,'dEc')-ingaasmaterial(x2,'dEhh');
   size(a2);
   a(index)=a2;
   index=find(x<=0); %formulas for InGaAs, implemented extremely rudimentary
   x2=x(index);
   s=size(x2);
   x2=-x2(:);
 % a2=1.517-5.5e-4*T^2/(T+200)+(0.42-2.76e-4*T^2/(T+83)-1.517+5.5e-4*T^2/(T+200))*x2-0.475*x2.*(1-x2);
   a2= 0.42 + 0.625*(1-x2) -( 5.8/(T+300)-4.19/(T+271) )*1e-4*T^2*(1-x2) - 4.19e-4*T^2/(T+271) + 0.475*(1-x2).^2  +ingaasmaterial(x2,'dEc')-ingaasmaterial(x2,'dEhh');
   a(index)=a2;
case {'E_G6G7'}		% G6-G7 = split-off
   index=find(x>0);
   x2=x(index);
   a2=1.786+1.028*x2+0.407*x2.*x2 +ingaasmaterial(x2,'dEc')-ingaasmaterial(x2,'dEso');
   a(index)=a2;
   index=find(x<=0); 
   x2=x(index);
   x2=-x2(:);
   a2=1.517+0.39+(0.341-0.39)*(1-x2)+ingaasmaterial(x2,'dEc')-ingaasmaterial(x2,'dEso');
   a(index)=a2;
case {'E_G8X6'}		% G8-X6 = ind Gap for x>.45
   index=find(x>0);
   x2=x(index);
   a2=1.900+1.25*x2+0.143*x2.*x2+ ingaasmaterial(x2,'dEc')-ingaasmaterial(x2,'dEhh');
   a(index)=a2;
   index=find(x<=0); 
   x2=x(index);
   x2=-x2(:);
   a2=1.981+0.39+(0.341-0.39)*(1-x2)+ingaasmaterial(x2,'dEc')-ingaasmaterial(x2,'dEhh'); 
   a(index)=a2;
case {'E_G8L6'}		% G8-L6 = ind Gap for x>.45  
   index=find(x>0);
   x2=x(index);
   a2=1.708+1.642*x2+ingaasmaterial(x2,'dEc')-ingaasmaterial(x2,'dEhh');
   a(index)=a2;
   index=find(x<=0); 
   x2=x(index);
   x2=-x2(:);
   a2=1.815+0.39+(0.341-0.39)*(1-x2)+ingaasmaterial(x2,'dEc')-ingaasmaterial(x2,'dEhh'); 
   a(index)=a2;    

   
case {'M_eG'}		% El-G-Mass (DOS-Mass=cond-Mass=t=l)
   index=find(x>0);
   a(index)=0.067+0.083*x(index);
   index=find(x<=0);
   %a(index)=0.023*ones(size(x(index)));
   x2=x(index);
   x2=-x2(:);
   a2=0.026+(0.067-0.026)*(1-x2);
   a(index)=a2;    
case {'M_eXd'}		% El-X-Mass (DOS)
   index=find(x>0);
   a(index)=(0.85-x(index)*0.14);
   index=find(x<=0);
   %a(index)=0.023*ones(size(x(index)));
   x2=x(index);
   x2=-x2(:);
   a2=0.64+(0.85-0.64)*(1-x2);
   a(index)=a2;  
case {'M_eXl'}		% El-X-Mass long
   index=find(x>0);
   a(index)=(1.3-x(index)*(1.3-1.1));
   index=find(x<=0);
   %a(index)=0.023*ones(size(x(index)));
   x2=x(index);
   x2=-x2(:);
   a2=1.13+(1.3-1.13)*(1-x2);
   a(index)=a2; 
case {'M_eXt'}		% El-X-Mass trans
   index=find(x>0);
   a(index)=(0.23-x(index)*(0.23-0.19));
   index=find(x<=0);
   %a(index)=0.023*ones(size(x(index)));
   x2=x(index);
   x2=-x2(:);
   a2=0.16+(0.23-0.16)*(1-x2);
   a(index)=a2; 
case {'M_eX'}		% El-X-Mass
   ext=ingaasmaterial(x,'M_eXt');
   exl=ingaasmaterial(x,'M_eXl');
   a=3.0*(ext.*exl./(ext+2.0*exl));
   
case {'M_eLd'}		% El-L-Mass (DOS)
   index=find(x>0);
   a(index)=(0.56+x(index)*0.22);
   index=find(x<=0);
   %a(index)=0.023*ones(size(x(index)));
   x2=x(index);
   x2=-x2(:);
   a2=0.29+(0.56-0.29)*(1-x2);
   a(index)=a2; 
case {'M_eLl'}		% El-L-Mass long
   index=find(x>0);
   a(index)=(1.9-x(index)*(1.9-1.32));
   index=find(x<=0);
   %a(index)=0.023*ones(size(x(index)));
   x2=x(index);
   x2=-x2(:);
   a2=0.64+(1.9-0.64)*(1-x2);
   a(index)=a2;
case {'M_eLt'}		% El-L-Mass trans
   index=find(x>0);
   a(index)=(0.0754+x(index)*(0.15-0.0754));
   index=find(x<=0);
   %a(index)=0.023*ones(size(x(index)));
   x2=x(index);
   x2=-x2(:);
   a2=0.05+(0.0754-0.05)*(1-x2);
   a(index)=a2;
case {'M_eL'}		% El-L-Mass
   elt=ingaasmaterial(x,'M_eLt');
   ell=ingaasmaterial(x,'M_eLl');
   a=3.0*(elt.*ell./(elt+2.0*ell));
   
case {'M_hh001'}	% HH-Mass 001
   g1=ingaasmaterial(x,'Gamma1');
   g2=ingaasmaterial(x,'Gamma2');
   a=1./(g1-2*g2);
case {'M_hh110'}	% HH-Mass 110
   g1=ingaasmaterial(x,'Gamma1');
   g2=ingaasmaterial(x,'Gamma2');
   g3=ingaasmaterial(x,'Gamma3');
   a=2./(2*g1-g2-3*g3);
case {'M_hh111'}	% HH-Mass 111
   g1=ingaasmaterial(x,'Gamma1');
   g3=ingaasmaterial(x,'Gamma3');
   a=1./(g1-2*g3);
case {'M_lh001'}	% LH-Mass 001
   g1=ingaasmaterial(x,'Gamma1');
   g2=ingaasmaterial(x,'Gamma2');
   a=1./(g1+2*g2);
case {'M_lh111'}	% LH-Mass 111
   g1=ingaasmaterial(x,'Gamma1');
   g3=ingaasmaterial(x,'Gamma3');
   a=1./(g1+2*g3);
case {'M_lh110'}	% LH-Mass 110
   g1=ingaasmaterial(x,'Gamma1');
   g2=ingaasmaterial(x,'Gamma2');
   g3=ingaasmaterial(x,'Gamma3');
   a=2./(2*g1+g2+3*g3);
case {'M_hhd'}		% HH-DOS-Mass
   index=find(x>0);
   x2=x(index);
   a2=(0.55+0.26*x2);
   a(index)=a2;
   index=find(x<=0); 
   x2=x(index);
   x2=-x2(:);
   a2=0.57+(0.55-0.57)*(1-x2); 
   a(index)=a2;
case {'M_lhd'}		% LH-DOS-Mass
  index=find(x>0);
   x2=x(index);
   a2=(0.08+0.08*x2);
   a(index)=a2;
   index=find(x<=0); 
   x2=x(index);
   x2=-x2(:);
   a2=0.025+(0.08-0.025)*(1-x2); 
   a(index)=a2;
case {'M_so'}		% SO-Mass
  index=find(x>0);
   x2=x(index);
   a2=(0.165+0.135*x2);      
   a(index)=a2;
   index=find(x<=0); 
   x2=x(index);
   x2=-x2(:);
   a2=0.14+(0.165-0.14)*(1-x2); 
   a(index)=a2;

   
case {'Gamma1'}		% Lutt Param
   index=find(x>0);
   x2=x(index);
   a2=1.0./((1.0-x2)./7.10+x2./3.76);
   a(index)=a2;
   index=find(x<=0); 
   x2=x(index);
   x2=-x2(:);
   a2=20.0+(6.98-20)*(1-x2); 
   a(index)=a2;
case {'Gamma2'}		% Lutt Param
   index=find(x>0);
   x2=x(index);
   a2=1.0./((1.0-x2)./2.02+x2./0.90);
   a(index)=a2;
   index=find(x<=0); 
   x2=x(index);
   x2=-x2(:);
   a2=8.5+(2.06-8.5)*(1-x2); 
   a(index)=a2;
case {'Gamma3'}		% Lutt Param
   index=find(x>0);
   x2=x(index);
   a2=1.0./((1.0-x2)./2.91+x2./1.42);
   a(index)=a2;
   index=find(x<=0); 
   x2=x(index);
   x2=-x2(:);
   a2=9.2+(2.93-9.2)*(1-x2);
   a(index)=a2;

case {'dEc'}	
   Ga_frac_barrier=aquila_material.barrier;
   aalattice_GaAs=5.65325+3.88e-5*(T-300);
   aalattice_InAs=6.0583+2.74e-5*(T-300);
   aalattice_barriar=aalattice_InAs+(aalattice_GaAs-aalattice_InAs)*Ga_frac_barrier;	
   index=find(x>0);
   x2=x(index);
   a2=0;
   a(index)=a2;
   index=find(x<=0); 
   x2=x(index);
   x2=-x2(:);
   aalattice=aalattice_InAs+(aalattice_GaAs-aalattice_InAs)*(1-x2);
   cc11=832.9+(1221-832.9)*(1-x2); 
   cc12=452.6+(566-452.6)*(1-x2);
   aav=1.00+(1.16-1.00)*(1-x2); 
   aac=-5.08+(-7.17+5.08)*(1-x2);
   bbbb=-1.8+(-2.0+1.8)*(1-x2);
   sstrain=(aalattice_barriar-aalattice)/aalattice_barriar;
   sstrainz=-2*cc12/cc11*(aalattice_GaAs-aalattice)/aalattice_barriar;
   a2=aac.*(2*sstrain+sstrainz);  %% Ec shift
   a(index)=a2;



case {'dEhh'}	
   Ga_frac_barrier=aquila_material.barrier;
   aalattice_GaAs=5.65325+3.88e-5*(T-300);
   aalattice_InAs=6.0583+2.74e-5*(T-300);
   aalattice_barriar=aalattice_InAs+(aalattice_GaAs-aalattice_InAs)*Ga_frac_barrier;	
   index=find(x>0);
   x2=x(index);
   a2=0;
   a(index)=a2;
   index=find(x<=0); 
   x2=x(index);
   x2=-x2(:);
   aalattice=aalattice_InAs+(aalattice_GaAs-aalattice_InAs)*(1-x2);
   cc11=832.9+(1221-832.9)*(1-x2); 
   cc12=452.6+(566-452.6)*(1-x2);
   aav=1.00+(1.16-1.00)*(1-x2); 
   aac=-5.08+(-7.17+5.08)*(1-x2);
   bbbb=-1.8+(-2.0+1.8)*(1-x2);
   sstrain=(aalattice_barriar-aalattice)./aalattice_barriar;
   sstrainz=-2*cc12./cc11.*(aalattice_GaAs-aalattice)./aalattice_barriar;
   a2=aav.*(2*sstrain+sstrainz)+bbbb./2.*(2*sstrain-2*sstrainz) ; %% Ehh shift
   a(index)=a2;

case {'dElh'}	
   Ga_frac_barrier=aquila_material.barrier;
   aalattice_GaAs=5.65325+3.88e-5*(T-300);
   aalattice_InAs=6.0583+2.74e-5*(T-300);
   aalattice_barriar=aalattice_InAs+(aalattice_GaAs-aalattice_InAs)*Ga_frac_barrier;	
   index=find(x>0);
   x2=x(index);
   a2=0;
   a(index)=a2;
   index=find(x<=0); 
   x2=x(index);
   x2=-x2(:);
   aalattice=aalattice_InAs+(aalattice_GaAs-aalattice_InAs)*(1-x2);
   cc11=832.9+(1221-832.9)*(1-x2); 
   cc12=452.6+(566-452.6)*(1-x2);
   aav=1.00+(1.16-1.00)*(1-x2); 
   aac=-5.08+(-7.17+5.08)*(1-x2);
   bbbb=-1.8+(-2.0+1.8)*(1-x2);
   sstrain=(aalattice_barriar-aalattice)./aalattice_barriar;
   sstrainz=-2*cc12./cc11.*(aalattice_GaAs-aalattice)./aalattice_barriar;
   Dso=0.39+(0.341-0.39)*(1-x2);
   Qe=-bbbb./2.*(2*sstrain-2*sstrainz);
   R=Qe-Dso+sqrt( Dso.*Dso+9*Qe.*Qe +2*Qe.*Dso) ;
 % a2=+aav*(2*sstrain+sstrainz)-bbbb/2*(2*sstrain-2*sstrainz); %% Elh shift
   a2=+aav.*(2*sstrain+sstrainz)+ 0.5*R ; %% Elh shift
   a(index)=a2;

case {'dEso'}	
   Ga_frac_barrier=aquila_material.barrier;
   aalattice_GaAs=5.65325+3.88e-5*(T-300);
   aalattice_InAs=6.0583+2.74e-5*(T-300);
   aalattice_barriar=aalattice_InAs+(aalattice_GaAs-aalattice_InAs)*Ga_frac_barrier;	
   index=find(x>0);
   x2=x(index);
   a2=0;
   a(index)=a2;
   index=find(x<=0); 
   x2=x(index);
   x2=-x2(:);
   aalattice=aalattice_InAs+(aalattice_GaAs-aalattice_InAs)*(1-x2);
   cc11=832.9+(1221-832.9)*(1-x2); 
   cc12=452.6+(566-452.6)*(1-x2);
   aav=1.00+(1.16-1.00)*(1-x2); 
   aac=-5.08+(-7.17+5.08)*(1-x2);
   bbbb=-1.8+(-2.0+1.8)*(1-x2);
   sstrain=(aalattice_barriar-aalattice)./aalattice_barriar;
   sstrainz=-2*cc12./cc11.*(aalattice_GaAs-aalattice)./aalattice_barriar;
   Dso=0.39+(0.341-0.39)*(1-x2);
   Qe=-bbbb./2.*(2*sstrain-2*sstrainz);
   R=Qe-Dso-sqrt( Dso.*Dso+9*Qe.*Qe +2*Qe.*Dso) ;
 % a2=+aav*(2*sstrain+sstrainz)-bbbb/2*(2*sstrain-2*sstrainz); %% Elh shift
   a2=+aav.*(2*sstrain+sstrainz)+ 0.5*R ; %% Elh shift
   a(index)=a2;





case {'E_c'}		% CB-Edge
   p=find(x<0.45);
   eg=ingaasmaterial(x,'E_g');
   a(p)=eg(p);
   p=find(x>=0.45);
   ex=ingaasmaterial(x,'E_x');
   a(p)=ex(p);
case {'E_v'}		% VB-Edge
   a=ingaasmaterial(x,'E_g')-ingaasmaterial(x,'E_G6G8');
case {'E_x'}		% X-Point
   a=ingaasmaterial(x,'E_v')+ingaasmaterial(x,'E_G8X6');
case {'E_l'}	   % L-Point
   a=ingaasmaterial(x,'E_v')+ingaasmaterial(x,'E_G8L6');
case {'E_g'}		% G-Point
   a=ingaasmaterial(x,'E_G6G8')-1.517-ingaasmaterial(x,'Delta_Ev');		% therefore Eg=0 for x=0
case {'E_so'}		% SO-Edge
   a=ingaasmaterial(x,'E_g')-ingaasmaterial(x,'E_G6G7');
   
case {'Delta_Ev'}	% VB-Offset
   index=find(x>0);
   a(index)=0.51*x(index);
   index=find(x<=0);
   p=1.5+(4.88-1.5)/.9*(-x(index)-0.1);
   a(index)=-(ingaasmaterial(zeros(size(x(index))),'E_G6G8')-ingaasmaterial(x(index),'E_G6G8'))./(1+p);      
case {'eps'}		% Epsilon
   index=find(x>0);
   a(index)=13.18-3.12*x(index);
   index=find(x<=0);
   x2=x(index);
   x2=-x2(:);
   a(index)= 15.1-2.87*(1-x2)+0.67*(1-x2).^2;


case {'E_d'}		% Donatorlevel
   ec=ingaasmaterial(x,'E_c');
   p=find(x<0.22);
   a(p)=ec(p)-.006;
   p=find((x>=0.22)&(x<0.35));
   a(p)=ec(p)-(0.006+(x(p)-0.22)/(.35-.22)*(0.019-0.006));
   p=find(x>=0.35);   
   a(p)=ec(p)-0.019;
otherwise
   error('material: property undefined!');
end
