function c = transport(a,n,dt,nsteps,mflag,nplot)

% Compute the solution of the pure advection equation u_t - a u_x = 0
% with homogeneous Dirichlet boundary conditions, starting from the
% initial solution  u0(x). 
% a   is the advection velocity (a > 0 is assumed).
% n   is the number of unknowns at x1=h to xn=1 on the interval [0,1].
% dt  is the time step.
% nsteps  is the number of time steps.
% nplot   plot the solution every nplot time steps.
% mflag = 0  Forward-time backward-space FD scheme
% mflag = 1  Explicite centered-space FD
% mflag = 2  Implicite centered-space FD
% mflag = 3 Lax-wendroff
h = 1/n;    %space step
% initial solution
  x = h:h:1;
  u = u0(x)';
  xx = 0:h:1;
  yy = u0(xx);
% computation of the transfer matrix
  bet= a*dt/h;
  switch mflag
      
  case 0 
  %Forward-time backward-space FD scheme
  D = sparse(1:n,1:n,(1-bet)*ones(n,1),n,n);
  L = sparse(2:n,1:n-1,bet*ones(n-1,1),n,n);
  A = L+D;
  case 1
  % Explicit centered-space FD scheme
  D = sparse(1:n,1:n,ones(n,1),n,n);
  L = sparse(2:n,1:n-1,bet/2*ones(n-1,1),n,n);
  U = sparse(1:n-1,2:n,-bet/2*ones(n-1,1),n,n);
  A = L+D+U;
  case 2
  % Implicit centered-space FD
  D = sparse(1:n,1:n,ones(n,1),n,n);
  L = sparse(2:n,1:n-1,-1*(bet/2)*ones(n-1,1),n,n);
  U = sparse(1:n-1,2:n,(bet/2)*ones(n-1,1),n,n);
  A = L+D+U;
  case 3
  %Lax-wendroff scheme
  D = sparse(1:n,1:n,(1- bet*bet)*ones(n,1),n,n);
  L = sparse(2:n,1:n-1,0.5*(bet+ bet*bet)*ones(n-1,1),n,n);
  U = sparse(1:n-1,2:n,0.5*(-bet+ bet*bet)*ones(n-1,1),n,n);
  A = L+D+U;
  end; 
for iter = 1:nsteps
    if (mod(iter,nplot)==0)
      yy(2:n+1) = u(1:n,1)';
      plot(xx,yy,'b')
      axis([0,1,-1,1.5]);
      title(sprintf('Iteration number %u',iter));
      FF = getframe;
    end;
    switch mflag
      
      case 0
        % Forward-time backward-space FD scheme
        u =  A*u + dt*f(x)';

      case 1
        % Explicit centered-space FD
        u =  A*u + dt*f(x)';
      
      case 2
        % Implicit centered-space FD
        u = A\(u + dt*f(x)');
      
      case 3
        %Lax-wendroff scheme
        u = A*u +dt*f(x)';
    end;   
    
  end;

  c = u;

%The functions defining the source term and intial conditions

 function f=f(x)
 f=0.0;
 %f=exp(-200*(x-0.25).^2);
 
 function u0=u0(x)
 u0= 0*(x<1/8)+ 1*(x>=1/8).*(x<=1/4)+0*(x>1/4);
 %u0=sin(20*x);
 %u0=exp(-100*(x-0.25).^2);