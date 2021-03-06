function [V] =  SSpotential(Nx,Ny,Lx,Ly,R,Rs,eps, a)

x = zeros(1,Nx);
y = zeros(1,Ny);

dx = Lx / Nx;
dy = Ly / Ny;

%% Shift x and y
for i = 0:Nx-1 
  if( i <= Nx/2)
    x(i + 1) = dx * i;
  else
    x(i + 1) = dx * ( i - Nx); 
  end
end

for i = 0: Ny-1
  if( i <= Ny/2)
    y(i + 1) = dy * i;
  else
    y(i + 1) = dy * ( i - Ny); 
  end
end

[y2, x2] =  meshgrid(y,x);

V = eps * ( exp( -  ( (x2.^2 + y2 .^2 ) / R .^2 ) .^ 4 ) + ...
a * exp( -  ( (x2.^2 + y2 .^2 ) / Rs .^2 ) .^ 4 )  );

end
