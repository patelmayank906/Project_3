%%The following code is used to calculate the C-matrix i.e. damping ratio


load('reduced_K_M')

OPTS.issym=1;
  OPTS.isreal=1;
  %To increase accuracy: (Agnes)
  OPTS.tol=eps/10;
  [minvals,minvallocs]=sort(diag(K_r)./diag(M_r));
  shift=minvals(min(7,length(minvals)));
  [fms,f]=eigs((K_r+K_r')/2+shift*(M_r+M_r')/2,(M_r+M_r')/2,min([ size(K_r,1) ...
		    max([floor(sqrt(size(K_r,1))) 150])]),0,OPTS);
  fs=sqrt(diag(f)-shift)/2/pi;
  
  zeta = 0.02;
  zeta_omega1 = 2*zeta*fs*2*pi;
  zeta_omega = diag(zeta_omega1);
  C = fms'\zeta_omega/fms;
  I4 = fms'*M_r*fms      % To check the C_matrix