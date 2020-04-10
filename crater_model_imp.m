%%%
%%% Crater Defect Formation Model (Evans et al., J. Colloid Interface Sci., 2000)
%%% Started by nikoscham (2020)
%%%

tic
clear all

global Ht L co so sinf Akink B Mc eo G1 xo Tscale cdry A d dr dt;

%%%
%%% Model Parameters
%%%
Ht = 0.003; %%% Initial film thickness (cm)
L = 0.025; %%% Horizontal length scale (cm)
co = 0.6; %%% Initial resin fraction
so = 30.0; %%% Reference surface tension (dyn/cm)
sinf = 22.0; %%% Minimum surface tension (dyn/cm)
Akink = 10.0; %%% Surface tension law 'kink'
%Ds = 8.6e-6; %%% Surfactant diffusivity (cm2/s)
Ds = 1.0e-3; %%% Surfactant diffusivity (cm2/s)
mo = 5.0; %%% Initial paint viscosity (P)
B = 4.0; %%% Viscosity law constant
Mc = 25; %%% Viscosity law constant
%eo = 2e-6; %%% Drying rate (cm/s)
eo = 6.7e-4; %%% Drying rate (cm/s) - Valspar EzDex 4000W56R/11HB
Tscale = (3.0*mo*L^4)/(so*Ht^3); %%% Time scale
d = Tscale*Ds/L^2; %%% Dimensionless diffusion parameter
G0 = 1.0 ; %%% Bolus concentration
G1 = G0/2.0; %%% Critical micelle concentration
xo = -Akink*G1; %%% For surface tension computation
A = (3.0/2.0)*(so-sinf)/so*(L/Ht)^2; %%% Dimensionless aspect ratio
E = Tscale*eo/Ht;
dryt = (1.0-co)/E %%% Drying time (s)
cdry = 0.995; %%% Dry resin concentration
dryt = 15;

%%%
%%% Set up the mesh
%%%
N = 200;
rend = 12.0;
r = linspace(0,rend,N+1)';
%r = rend*pow2(0:-0.1:-10)';
dr = rend/real(N);
hjac = 1.0e-4*dr; %%% Used for jacobian calc

%%%
%%% Set up time step
%%%
dt = 1.0e-4;

%%%
%%% Initialize variables
%%%
k = 1;
t(k)=0;
disp ('Time step:'), disp(k)
disp ('Time (s):'), disp(t(k)*Tscale)
Amat = sparse(3*N+3,3*N+3);
for i = 1:N+1
  h(i,1) = 1.0; %%% Film thickness
  c(i,1) = co; %%% Resin concentration
  hr(i,1) = c(i,1)*h(i,1); %%% Resin depth
%  G(i,1) = 1.0/(0.1+(r(i)/L)); %%% Surfactant concentration
  G(i,1) = (1.0-erf((r(i)-1.0)*2.5))/2.0; %%% Surfactant concentration
%  G(i,1) =  G0*exp(-(r(i)/L)^2); %%% Surfactant concentration
end
G(1,1) = G(2,1); %%% Surfactant concentration

while ((t(k)*Tscale) <= dryt) %%% Loop over time
  k = k+1;
  disp ('Time step:'), disp(k)
  t(k) = t(k-1)+dt;
  disp ('Time (s):'), disp(t(k)*Tscale)
  
  for i = 1:N+1  %%% Solution guessing
    h(i,k) = h(i,k-1);
	G(i,k) = G(i,k-1);
	hr(i,k) = hr(i,k-1);
    c(i,k) = hr(i,k)./h(i,k);
  end

%%%
%%% Start Newton iterations
%%%  
  newtonerror = 1.0;
  while (newtonerror >= 1e-6)

	Amat(:,:) = 0.0;
	equ = 1; %%% Film thickness Evolution
    for i = 2:N %%% Loop over space

      if (i == 2) %%% Semi-forward differences
        [res(i,k)] = rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ);
        Amat(i,i) = (rescomp2(h(i,k)+hjac, G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ) - rescomp2(h(i,k)-hjac, G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ))/(2.0*hjac);
        Amat(i,i+1) = (rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k)+hjac, G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ) - rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k)-hjac, G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ))/(2.0*hjac);
        Amat(i,i-1) = (rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k)+hjac, G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ) - rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k)-hjac, G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ))/(2.0*hjac);
        Amat(i,i+2) = (rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k)+hjac, equ) - rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k)-hjac, equ))/(2.0*hjac);
        Amat(i,i+3) = (rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k)+hjac, h(i+2,k), equ) - rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k)-hjac, h(i+2,k), equ))/(2.0*hjac);
        Amat(i,i+4) = (rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k)+hjac, h(i+3,k), h(i+2,k), equ) - rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k)-hjac, h(i+3,k), h(i+2,k), equ))/(2.0*hjac);
		
        Amat(i,N+1+i) = (rescomp2(h(i,k), G(i,k)+hjac, hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ) - rescomp2(h(i,k), G(i,k)-hjac, hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ))/(2.0*hjac);
        Amat(i,N+1+i+1) = (rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k)+hjac, hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ) - rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k)-hjac, hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ))/(2.0*hjac);
        Amat(i,N+1+i-1) = (rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k)+hjac, hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ) - rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k)-hjac, hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ))/(2.0*hjac);
		
        Amat(i,2*N+2+i) = (rescomp2(h(i,k), G(i,k), hr(i,k)+hjac, h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ) - rescomp2(h(i,k), G(i,k), hr(i,k)-hjac, h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ))/(2.0*hjac);
        Amat(i,2*N+2+i+1) = (rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k)+hjac, r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ) - rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k)-hjac, r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ))/(2.0*hjac);
        Amat(i,2*N+2+i-1) = (rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k)+hjac, h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ) - rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k)-hjac, h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ))/(2.0*hjac);

	  elseif (i == N) %%% Semi-backward differences
        [res(i,k)] = rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ);
        Amat(i,i) = (rescompn(h(i,k)+hjac, G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ)-rescompn(h(i,k)-hjac, G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ))/(2.0*hjac);
        Amat(i,i-1) = (rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k)+hjac, G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ)-rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k)-hjac, G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ))/(2.0*hjac);
        Amat(i,i+1) = (rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k)+hjac, G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ)-rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k)-hjac, G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ))/(2.0*hjac);
        Amat(i,i-2) = (rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k)+hjac, equ)-rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k)-hjac, equ))/(2.0*hjac);
        Amat(i,i-3) = (rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k)+hjac, h(i-2,k), equ)-rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k)-hjac, h(i-2,k), equ))/(2.0*hjac);
        Amat(i,i-4) = (rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k)+hjac, h(i-3,k), h(i-2,k), equ)-rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k)-hjac, h(i-3,k), h(i-2,k), equ))/(2.0*hjac);
		
        Amat(i,N+1+i) = (rescompn(h(i,k), G(i,k)+hjac, hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ)-rescompn(h(i,k), G(i,k)-hjac, hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ))/(2.0*hjac);
        Amat(i,N+1+i-1) = (rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k)+hjac, hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ)-rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k)-hjac, hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ))/(2.0*hjac);
        Amat(i,N+1+i+1) = (rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k)+hjac, hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ)-rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k)-hjac, hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ))/(2.0*hjac);
		
        Amat(i,2*N+2+i) = (rescompn(h(i,k), G(i,k), hr(i,k)+hjac, h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ)-rescompn(h(i,k), G(i,k), hr(i,k)-hjac, h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ))/(2.0*hjac);
        Amat(i,2*N+2+i-1) = (rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k)+hjac, h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ)-rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k)-hjac, h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ))/(2.0*hjac);
        Amat(i,2*N+2+i+1) = (rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k)+hjac, r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ)-rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k)-hjac, r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ))/(2.0*hjac);
		
	  else %%% Central differences
        [res(i,k)] = rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ);
        Amat(i,i) = (rescomp(h(i,k)+hjac, G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ) - rescomp(h(i,k)-hjac, G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ))/(2.0*hjac);
		Amat(i,i+1) = (rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k)+hjac, G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ) - rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k)-hjac, G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ))/(2.0*hjac);
		Amat(i,i+2) = (rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k)+hjac, equ) - rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k)-hjac, equ))/(2.0*hjac);
		Amat(i,i-1) = (rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k)+hjac, G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ) - rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k)-hjac, G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ))/(2.0*hjac);
		Amat(i,i-2) = (rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k)+hjac, h(i+2,k), equ) - rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k)-hjac, h(i+2,k), equ))/(2.0*hjac);
		
        Amat(i,N+1+i) = (rescomp(h(i,k), G(i,k)+hjac, hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ) - rescomp(h(i,k), G(i,k)-hjac, hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ))/(2.0*hjac);
		Amat(i,N+1+i+1) = (rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k)+hjac, hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ) - rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k)-hjac, hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ))/(2.0*hjac);
		Amat(i,N+1+i-1) = (rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k)+hjac, hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ) - rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k)-hjac, hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ))/(2.0*hjac);
		
        Amat(i,2*N+2+i) = (rescomp(h(i,k), G(i,k), hr(i,k)+hjac, h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ) - rescomp(h(i,k), G(i,k), hr(i,k)-hjac, h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ))/(2.0*hjac);
		Amat(i,2*N+2+i+1) = (rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k)+hjac, r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ) - rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k)-hjac, r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ))/(2.0*hjac);
		Amat(i,2*N+2+i-1) = (rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k)+hjac, h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ) - rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k)-hjac, h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ))/(2.0*hjac);
	  end
	
    end
	
	equ = 2; %%% Surfactant concentration
	i = 2;
    for jj = (N+3):(2*N+1) %%% Loop over space
	
      if (i == 2) %%% Semi-forward differences
        [res(jj,k)] = rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ);
        Amat(jj,i) = (rescomp2(h(i,k)+hjac, G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ) - rescomp2(h(i,k)-hjac, G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ))/(2.0*hjac);
        Amat(jj,i+1) = (rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k)+hjac, G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ) - rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k)-hjac, G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ))/(2.0*hjac);
        Amat(jj,i-1) = (rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k)+hjac, G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ) - rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k)-hjac, G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ))/(2.0*hjac);
        Amat(jj,i+2) = (rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k)+hjac, equ) - rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k)-hjac, equ))/(2.0*hjac);
        Amat(jj,i+3) = (rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k)+hjac, h(i+2,k), equ) - rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k)-hjac, h(i+2,k), equ))/(2.0*hjac);
        Amat(jj,i+4) = (rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k)+hjac, h(i+3,k), h(i+2,k), equ) - rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k)-hjac, h(i+3,k), h(i+2,k), equ))/(2.0*hjac);
		
        Amat(jj,N+1+i) = (rescomp2(h(i,k), G(i,k)+hjac, hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ) - rescomp2(h(i,k), G(i,k)-hjac, hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ))/(2.0*hjac);
        Amat(jj,N+1+i+1) = (rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k)+hjac, hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ) - rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k)-hjac, hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ))/(2.0*hjac);
        Amat(jj,N+1+i-1) = (rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k)+hjac, hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ) - rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k)-hjac, hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ))/(2.0*hjac);
		
        Amat(jj,2*N+2+i) = (rescomp2(h(i,k), G(i,k), hr(i,k)+hjac, h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ) - rescomp2(h(i,k), G(i,k), hr(i,k)-hjac, h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ))/(2.0*hjac);
        Amat(jj,2*N+2+i+1) = (rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k)+hjac, r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ) - rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k)-hjac, r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ))/(2.0*hjac);
        Amat(jj,2*N+2+i-1) = (rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k)+hjac, h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ) - rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k)-hjac, h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ))/(2.0*hjac);

	  elseif (i == N) %%% Semi-backward differences
        [res(jj,k)] = rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ);
        Amat(jj,i) = (rescompn(h(i,k)+hjac, G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ)-rescompn(h(i,k)-hjac, G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ))/(2.0*hjac);
        Amat(jj,i-1) = (rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k)+hjac, G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ)-rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k)-hjac, G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ))/(2.0*hjac);
        Amat(jj,i+1) = (rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k)+hjac, G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ)-rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k)-hjac, G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ))/(2.0*hjac);
        Amat(jj,i-2) = (rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k)+hjac, equ)-rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k)-hjac, equ))/(2.0*hjac);
        Amat(jj,i-3) = (rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k)+hjac, h(i-2,k), equ)-rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k)-hjac, h(i-2,k), equ))/(2.0*hjac);
        Amat(jj,i-4) = (rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k)+hjac, h(i-3,k), h(i-2,k), equ)-rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k)-hjac, h(i-3,k), h(i-2,k), equ))/(2.0*hjac);
		
        Amat(jj,N+1+i) = (rescompn(h(i,k), G(i,k)+hjac, hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ)-rescompn(h(i,k), G(i,k)-hjac, hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ))/(2.0*hjac);
        Amat(jj,N+1+i-1) = (rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k)+hjac, hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ)-rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k)-hjac, hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ))/(2.0*hjac);
        Amat(jj,N+1+i+1) = (rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k)+hjac, hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ)-rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k)-hjac, hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ))/(2.0*hjac);
		
        Amat(jj,2*N+2+i) = (rescompn(h(i,k), G(i,k), hr(i,k)+hjac, h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ)-rescompn(h(i,k), G(i,k), hr(i,k)-hjac, h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ))/(2.0*hjac);
        Amat(jj,2*N+2+i-1) = (rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k)+hjac, h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ)-rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k)-hjac, h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ))/(2.0*hjac);
        Amat(jj,2*N+2+i+1) = (rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k)+hjac, r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ)-rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k)-hjac, r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ))/(2.0*hjac);
		
	  else %%% Central differences
        [res(jj,k)] = rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ);
        Amat(jj,i) = (rescomp(h(i,k)+hjac, G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ) - rescomp(h(i,k)-hjac, G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ))/(2.0*hjac);
		Amat(jj,i+1) = (rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k)+hjac, G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ) - rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k)-hjac, G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ))/(2.0*hjac);
		Amat(jj,i+2) = (rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k)+hjac, equ) - rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k)-hjac, equ))/(2.0*hjac);
		Amat(jj,i-1) = (rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k)+hjac, G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ) - rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k)-hjac, G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ))/(2.0*hjac);
		Amat(jj,i-2) = (rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k)+hjac, h(i+2,k), equ) - rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k)-hjac, h(i+2,k), equ))/(2.0*hjac);
		
        Amat(jj,N+1+i) = (rescomp(h(i,k), G(i,k)+hjac, hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ) - rescomp(h(i,k), G(i,k)-hjac, hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ))/(2.0*hjac);
		Amat(jj,N+1+i+1) = (rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k)+hjac, hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ) - rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k)-hjac, hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ))/(2.0*hjac);
		Amat(jj,N+1+i-1) = (rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k)+hjac, hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ) - rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k)-hjac, hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ))/(2.0*hjac);
		
        Amat(jj,2*N+2+i) = (rescomp(h(i,k), G(i,k), hr(i,k)+hjac, h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ) - rescomp(h(i,k), G(i,k), hr(i,k)-hjac, h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ))/(2.0*hjac);
		Amat(jj,2*N+2+i+1) = (rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k)+hjac, r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ) - rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k)-hjac, r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ))/(2.0*hjac);
		Amat(jj,2*N+2+i-1) = (rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k)+hjac, h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ) - rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k)-hjac, h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ))/(2.0*hjac);
	  end
	  i = i + 1;
	end

	equ = 3; %%% Resin depth
	i = 2;
	for kk = (2*N+4):(3*N+2) %%% Loop over space
      if (i == 2) %%% Semi-forward differences
        [res(kk,k)] = rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ);
        Amat(kk,i) = (rescomp2(h(i,k)+hjac, G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ) - rescomp2(h(i,k)-hjac, G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ))/(2.0*hjac);
        Amat(kk,i+1) = (rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k)+hjac, G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ) - rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k)-hjac, G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ))/(2.0*hjac);
        Amat(kk,i-1) = (rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k)+hjac, G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ) - rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k)-hjac, G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ))/(2.0*hjac);
        Amat(kk,i+2) = (rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k)+hjac, equ) - rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k)-hjac, equ))/(2.0*hjac);
        Amat(kk,i+3) = (rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k)+hjac, h(i+2,k), equ) - rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k)-hjac, h(i+2,k), equ))/(2.0*hjac);
        Amat(kk,i+4) = (rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k)+hjac, h(i+3,k), h(i+2,k), equ) - rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k)-hjac, h(i+3,k), h(i+2,k), equ))/(2.0*hjac);
		
        Amat(kk,N+1+i) = (rescomp2(h(i,k), G(i,k)+hjac, hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ) - rescomp2(h(i,k), G(i,k)-hjac, hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ))/(2.0*hjac);
        Amat(kk,N+1+i+1) = (rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k)+hjac, hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ) - rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k)-hjac, hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ))/(2.0*hjac);
        Amat(kk,N+1+i-1) = (rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k)+hjac, hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ) - rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k)-hjac, hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ))/(2.0*hjac);
		
        Amat(kk,2*N+2+i) = (rescomp2(h(i,k), G(i,k), hr(i,k)+hjac, h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ) - rescomp2(h(i,k), G(i,k), hr(i,k)-hjac, h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ))/(2.0*hjac);
        Amat(kk,2*N+2+i+1) = (rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k)+hjac, r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ) - rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k)-hjac, r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ))/(2.0*hjac);
        Amat(kk,2*N+2+i-1) = (rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k)+hjac, h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ) - rescomp2(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k)-hjac, h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i+4,k), h(i+3,k), h(i+2,k), equ))/(2.0*hjac);

	  elseif (i == N) %%% Semi-backward differences
        [res(kk,k)] = rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ);
        Amat(kk,i) = (rescompn(h(i,k)+hjac, G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ)-rescompn(h(i,k)-hjac, G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ))/(2.0*hjac);
        Amat(kk,i-1) = (rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k)+hjac, G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ)-rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k)-hjac, G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ))/(2.0*hjac);
        Amat(kk,i+1) = (rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k)+hjac, G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ)-rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k)-hjac, G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ))/(2.0*hjac);
        Amat(kk,i-2) = (rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k)+hjac, equ)-rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k)-hjac, equ))/(2.0*hjac);
        Amat(kk,i-3) = (rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k)+hjac, h(i-2,k), equ)-rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k)-hjac, h(i-2,k), equ))/(2.0*hjac);
        Amat(kk,i-4) = (rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k)+hjac, h(i-3,k), h(i-2,k), equ)-rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k)-hjac, h(i-3,k), h(i-2,k), equ))/(2.0*hjac);
		
        Amat(kk,N+1+i) = (rescompn(h(i,k), G(i,k)+hjac, hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ)-rescompn(h(i,k), G(i,k)-hjac, hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ))/(2.0*hjac);
        Amat(kk,N+1+i-1) = (rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k)+hjac, hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ)-rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k)-hjac, hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ))/(2.0*hjac);
        Amat(kk,N+1+i+1) = (rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k)+hjac, hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ)-rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k)-hjac, hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ))/(2.0*hjac);
		
        Amat(kk,2*N+2+i) = (rescompn(h(i,k), G(i,k), hr(i,k)+hjac, h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ)-rescompn(h(i,k), G(i,k), hr(i,k)-hjac, h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ))/(2.0*hjac);
        Amat(kk,2*N+2+i-1) = (rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k)+hjac, h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ)-rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k)-hjac, h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ))/(2.0*hjac);
        Amat(kk,2*N+2+i+1) = (rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k)+hjac, r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ)-rescompn(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k)-hjac, r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-4,k), h(i-3,k), h(i-2,k), equ))/(2.0*hjac);
		
	  else %%% Central differences
        [res(kk,k)] = rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ);
        Amat(kk,i) = (rescomp(h(i,k)+hjac, G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ) - rescomp(h(i,k)-hjac, G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ))/(2.0*hjac);
		Amat(kk,i+1) = (rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k)+hjac, G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ) - rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k)-hjac, G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ))/(2.0*hjac);
		Amat(kk,i+2) = (rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k)+hjac, equ) - rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k)-hjac, equ))/(2.0*hjac);
		Amat(kk,i-1) = (rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k)+hjac, G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ) - rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k)-hjac, G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ))/(2.0*hjac);
		Amat(kk,i-2) = (rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k)+hjac, h(i+2,k), equ) - rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k)-hjac, h(i+2,k), equ))/(2.0*hjac);
		
        Amat(kk,N+1+i) = (rescomp(h(i,k), G(i,k)+hjac, hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ) - rescomp(h(i,k), G(i,k)-hjac, hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ))/(2.0*hjac);
		Amat(kk,N+1+i+1) = (rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k)+hjac, hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ) - rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k)-hjac, hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ))/(2.0*hjac);
		Amat(kk,N+1+i-1) = (rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k)+hjac, hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ) - rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k)-hjac, hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ))/(2.0*hjac);
		
        Amat(kk,2*N+2+i) = (rescomp(h(i,k), G(i,k), hr(i,k)+hjac, h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ) - rescomp(h(i,k), G(i,k), hr(i,k)-hjac, h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ))/(2.0*hjac);
		Amat(kk,2*N+2+i+1) = (rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k)+hjac, r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ) - rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k), h(i+1,k), G(i+1,k), hr(i+1,k)-hjac, r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ))/(2.0*hjac);
		Amat(kk,2*N+2+i-1) = (rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k)+hjac, h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ) - rescomp(h(i,k), G(i,k), hr(i,k), h(i-1,k), G(i-1,k), hr(i-1,k)-hjac, h(i+1,k), G(i+1,k), hr(i+1,k), r(i), h(i,k-1), G(i,k-1), hr(i,k-1), h(i-2,k), h(i+2,k), equ))/(2.0*hjac);
	  end	
	  i = i + 1;    
    end
    
    res(1,k) = - 3.0*h(1,k) + 4.0*h(2,k) - 1.0*h(3,k); %%% No-flux boundary condition (N=0)
    res(N+1,k) = 3.0*h(N+1,k) - 4.0*h(N,k) + 1.0*h(N-1,k); %%% No-flux boundary condition (N=Lmax)
    Amat(1,1) = -3.0;
    Amat(1,2) = 4.0;
    Amat(1,3) = -1.0;
    Amat(N+1,N+1) = 3.0;
    Amat(N+1,N) = -4.0;
    Amat(N+1,N-1) = 1.0;
	
    res(N+2,k) = - 3.0*G(1,k) + 4.0*G(2,k) - 1.0*G(3,k); %%% No-flux boundary condition (N=0)
    res(2*N+2,k) = 3.0*G(N+1,k) - 4.0*G(N,k) + 1.0*G(N-1,k); %%% No-flux boundary condition (N=Lmax)
    Amat(N+2,N+2) = -3.0;
    Amat(N+2,N+3) = 4.0;
    Amat(N+2,N+4) = -1.0;
    Amat(2*N+2,2*N+2) = 3.0;
    Amat(2*N+2,2*N+1) = -4.0;
    Amat(2*N+2,2*N) = 1.0;
	
    res(2*N+3,k) = - 3.0*hr(1,k) + 4.0*hr(2,k) - 1.0*hr(3,k); %%% No-flux boundary condition (N=0)
    res(3*N+3,k) = 3.0*hr(N+1,k) - 4.0*hr(N,k) + 1.0*hr(N-1,k); %%% No-flux boundary condition (N=Lmax)
    Amat(2*N+3,2*N+3) = -3.0;
    Amat(2*N+3,2*N+4) = 4.0;
    Amat(2*N+3,2*N+5) = -1.0;
    Amat(3*N+3,3*N+3) = 3.0;
    Amat(3*N+3,3*N+2) = -4.0;
    Amat(3*N+3,3*N+1) = 1.0;		

	dsol = Amat\res(:,k); % Solve linear system
	
    for i = 1:N+1  %%% Update solution
	  h(i,k) = h(i,k) - dsol(i);
	  G(i,k) = G(i,k) - dsol(N+1+i);
	  hr(i,k) = hr(i,k) - dsol(2*N+2+i);
	  c(i,k) = hr(i,k)./h(i,k);
	  eoi(i,k) = eo*(1.0-erf((hr(i,k)./h(i,k)-0.7)*10))/2.0;
    end 
    
	newtonerror = 0.0;
    for i = 1:3*N+3
      newtonerror = newtonerror + res(i,k)^2;
    end
	newtonerror
    
  end
%%%
%%% Stop Newton iterations
%%%

%%%
%%% Plot results
%%%
  if (k == 2)
    subplot (1, 2, 1)
    p = plot (r*L*10.0, h(:,k)*Ht*10.0, 'linewidth', 1);
	hold on
%    p2 = plot (r*L*10.0, hr(:,k)*Ht*10.0, 'linewidth', 1);
	grid on;
    axis ([0, 10*L*10.0, 0, 0.10001], 'square');
    xlabel ('r (mm)');
    ylabel ('z (mm)');
    title ('Film thickness');
	
    subplot (1, 2, 2)
    pp = plot (r*L*10.0, G(:,k)/G0, 'linewidth', 1);
	grid on;
    axis ([0, 10*L*10.0, 0, 1.0], 'square');
    xlabel ('r (mm)');
    ylabel ('G/G_{init}');
    title ('Surfactant concentration');
	
    hhtext = text (1.3, 0.9, 't =');	
    htext = text (1.8, 0.9, num2str (t(k)*Tscale,2));
  else
    set (p, 'xdata', r*L*10.0);
    set (p, 'ydata', h(:,k)*Ht*10.0);
	
%    set (p2, 'xdata', r*L*10.0);
%    set (p2, 'ydata', hr(:,k)*Ht*10.0);
	
    set (pp, 'xdata', r*L*10.0);
    set (pp, 'ydata', G(:,k)/G0);
	
	delete(htext);
    htext = text (1.8, 0.9, num2str (t(k)*Tscale,2));
  end
  fname = sprintf ('img%03i.png', k);
  print ('-dpng', '-r100', fname);
 
%%%
%%% Update time step
%%% 
dt = dt*1.1;

end

toc