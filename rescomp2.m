function [res] = rescomp2(hi, Gi, hri, hmi, Gmi, hrmi, hpi, Gpi, hrpi, ri, hipre, Gipre, hripre, hppppi, hpppi, hppi, equ)
%%%
%%% Crater Defect Formation Model (Evans et al., J. Colloid Interface Sci., 2000)
%%% Started by nikoscham (2020)
%%%

%%% Residuals computatation using semi-forward differences

  global Ht L co so sinf Akink B Mc eo G1 xo Tscale cdry A d dr dt;

  ci = hri/hi; %%% Resin concentration (i)
  mi = ((1.0+B)*(ci/co)^Mc-B); %%% Paint viscosity (i)
  cmi = hrmi/hmi; %%% Resin concentration (i-1)
  mmi = ((1.0+B)*(cmi/co)^Mc-B); %%% Paint viscosity (i-1)
  cpi = hrpi/hpi; %%% Resin concentration (i+1)
  mpi = ((1.0+B)*(cpi/co)^Mc-B); %%% Paint viscosity (i+1)

  xi = Akink*(Gi-G1); %%% For surface tension computation (i)
  sdimi = (so-sinf)*((1.0+xi^2)^0.5-xi)/((1.0+xo^2)^0.5-xo)+sinf; %%% Surface tension (dimensional) (i-1)
  si = (sdimi-sinf)/(so-sinf); %%% Surface tension (dimensionless) (i)
  xmi = Akink*(Gmi-G1); %%% For surface tension computation (i-1)
  sdimmi = (so-sinf)*((1.0+xmi^2)^0.5-xmi)/((1.0+xo^2)^0.5-xo)+sinf; %%% Surface tension (dimensional) (i-1)
  smi = (sdimmi-sinf)/(so-sinf); %%% Surface tension (dimensionless) (i-1)
  xpi = Akink*(Gpi-G1); %%% For surface tension computation (i+1)
  sdimpi = (so-sinf)*((1.0+xpi^2)^0.5-xpi)/((1.0+xo^2)^0.5-xo)+sinf; %%% Surface tension (dimensional) (i+1)
  spi = (sdimpi-sinf)/(so-sinf); %%% Surface tension (dimensionless) (i+1)

%  if (ci >= cdry)
%    eoi = 0.0;
%  else
%    eoi = eo;
%  end
  eoi = eo*(1.0-erf((ci-0.7)*10))/2.0;
  Ei = Tscale*eoi/Ht; %%% Dimensionless parameter

%%%
%%% Finite difference schemes
%%%
%%% Film thickness
  hder3 = (-3.0*hmi+10.0*hi-12.0*hpi+6.0*hppi-1.0*hpppi)/(2.0*dr^3); 
  hder4 = (2.0*hmi-9.0*hi+16.0*hpi-14.0*hppi+6.0*hpppi-1.0*hppppi)/(dr^4); 
%%% Film thickness
  hder1 = (hpi-hmi)/(2.0*dr); 
  hder2 = (hpi-2.0*hi+hmi)/(dr^2);
%%% Surfactant concentration
  Gder1 = (Gpi-Gmi)/(2.0*dr); 
  Gder2 = (Gpi-2.0*Gi+Gmi)/(dr^2);
%%% Surface tension	  
  sder1 = (spi-smi)/(2.0*dr); 
  sder2 = (spi-2.0*si+smi)/(dr^2);
%%% Resin depth
  hrder1 = (hrpi-hrmi)/(2.0*dr);
%%% Viscosity
  mder1 = (mpi-mmi)/(2.0*dr);

  if (equ == 1) %%% Equation 1 - Film thickness Evolution
    resh = - hi + hipre + ((-1.0/ri)*(-3.0*hi^2*hder1^2*1.0/mi*1.0/ri + 1.0/mi^2*hi^3/ri*hder1*mder1 ...
    + hi^3/mi*1.0/ri^2*hder1 -  hi^3/mi*1.0/ri*hder2 ...
    + 3.0*hi^2*hder1*1.0/mi*hder2 - 1.0/mi^2*hi^3*hder2*mder1 + hi^3/mi*hder3 ...
    + hi^3/mi*hder3 + 3.0*hi^2*hder1*ri/mi*hder3 - 1.0/mi^2*ri*hi^3*hder3*mder1 + ri*hi^3/mi*hder4 ...
    + A*hi^2/mi*sder1 + 2.0*hi*hder1*ri*A*sder1 - 1.0/mi^2*ri*A*hi^2*sder1*mder1 + ri*A*hi^2/mi*sder2)-Ei)*dt;
    res = resh;
  elseif (equ == 2) %%% Equation 2 - Surfactant concentration
    resG = - Gi + Gipre + ((-1.0/ri)*(-Gder1*3.0/2.0*hi^2/mi*1.0/ri*hder1 - Gi*3.0/2.0*2.0*hi*hder1^2*1.0/mi*1.0/ri ...
    + Gi*3.0/2.0*hi^2*1.0/mi^2*1.0/ri*hder1*mder1 + Gi*3.0/2.0*hi^2/mi*1.0/ri^2*hder1 - Gi*3.0/2.0*hi^2/mi*1.0/ri*hder2 ...
    + Gder1*3.0/2.0*hi^2/mi*hder2 + Gi*3.0/2.0*2.0*hi*hder1*1.0/mi*hder2 - 1.0/mi^2*Gi*3.0/2.0*hi^2*hder2*mder1 + Gi*3.0/2.0*hi^2/mi*hder3 ...
    + Gder1*ri*3.0/2.0*hi^2/mi*hder3 + Gi*3.0/2.0*hi^2/mi*hder3 + Gi*ri*3.0/2.0*2.0*hi*hder1*1.0/mi*hder3 - 1.0/mi^2*mder1*Gi*ri*3.0/2.0*hi^2*hder3 + Gi*ri*3.0/2.0*hi^2/mi*hder4 ...
    + Gder1*ri*2.0*A*hi/mi*sder1 + Gi*2.0*A*hi/mi*sder1 + Gi*ri*2.0*A*hder1*1.0/mi*sder1 - 1.0/mi^2*Gi*ri*2.0*A*hi*sder1*mder1 +  Gi*ri*2.0*A*hi/mi*sder2)+d*(1.0/ri)*(Gder1+ri*Gder2))*dt;
    res = resG;
  elseif (equ == 3) %%% Equation 3 - Resin depth
    reshr = - hri + hripre + (-1.0/ri)*(-2.0*hi*hder1^2*1.0/mi*1.0/ri*hri + 1.0/mi^2*hi^2*1.0/ri*hder1*hri*mder1 ...
    + 1.0/ri^2*hi^2/mi*hder1*hri - hi^2/mi*1.0/ri*hder2*hri - hi^2/mi*1.0/ri*hder1*hrder1 ...
    + 2.0*hi*hder1*1.0/mi*hder2*hri - 1.0/mi^2*hi^2*hder2*hri*mder1 + hi^2/mi*hder3*hri + hi^2/mi*hder2*hrder1 ...
    + hi^2/mi*hder3*hri + 2.0*hi*hder1*ri/mi*hder3*hri - 1.0/mi^2*ri*hi^2*hder3*hri*mder1 + ri*hi^2/mi*hder4*hri ...
    + ri*hi^2/mi*hder3*hrder1 + A*hi/mi*sder1*hri + ri*A*1.0/mi*sder1*hri*hder1 - 1.0/mi^2*ri*A*hi*sder1*hri*mder1 ...
    + ri*A*hi/mi*sder2*hri + ri*A*hi/mi*sder1*hrder1)*dt;
    res = reshr;
  end

end
