%%The following code is used to calculate the C-matrix i.e. damping ratio

function [C,fs] = Damping(K_r,M_r,zeta)
OPTS.issym=1;
  OPTS.isreal=1;
  %To increase accuracy: (Agnes)
  OPTS.tol=eps/10;
  [minvals,minvallocs]=sort(diag(K_r)./diag(M_r));
  shift=minvals(min(2,length(minvals)));
  [fms,f]=eigs((K_r+K_r')/2+shift*(M_r+M_r')/2,(M_r+M_r')/2,min([ size(K_r,1) ...
		    max([floor(sqrt(size(K_r,1))) 2])]),0,OPTS);
  fs=sqrt(diag(f)-shift)/2/pi;
  
  zeta_omega1 = 2*zeta*fs*2*pi;
  zeta_omega = diag(zeta_omega1);
  C = fms'\zeta_omega/fms;
  

