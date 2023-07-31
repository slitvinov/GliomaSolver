%================================
%
%  get exact solution for Hypre tests
%
%================================



function [uexact, fAna] = getExactSolution(X,Y,testCase)

[nx,ny] = size(X);
uexact = zeros(nx,ny);
fAna = zeros(ny,ny);

switch testCase
    case {0}
        for ix = 1:nx
            for iy = 1:ny
                c = 4;
                uexact(ix,iy) = cos(c*pi*X(ix,iy))*cos(c*pi*Y(ix,iy));
                fAna(ix,iy) = -1*( -2*c*c*pi*pi + 1) * uexact(ix,iy);
                
                
            end;
        end;
        
        
        
    case {1}
        for ix = 1:nx
            for iy = 1:ny
                r = sqrt((X(ix,iy)-0.5)^2 + (Y(ix,iy)-0.5)^2);
                c = 4;
                alpha = c*pi*r;
                ifactor = 1./(c*c*pi*pi);
                uexact(ix,iy) =((sin(alpha)  - alpha * cos(alpha))* ifactor);
                fAna(ix,iy) =  -1*(sin(alpha)*(2 - ifactor) + alpha * cos(alpha)*(1+ ifactor));
            end;
        end;
        
    case {2}
        for ix = 1:nx
            for iy = 1:ny
                r = sqrt((X(ix,iy)-0.5)^2 + (Y(ix,iy)-0.5)^2);
                c = 4;
                alpha = c*pi*r;
                ifactor = 1./(c*c*pi*pi);
                uexact(ix,iy) =-1*((sin(alpha)  - alpha * cos(alpha))* ifactor) + 1;
                fAna(ix,iy) =  -1*(sin(alpha)*(2 - ifactor) + alpha * cos(alpha)*(1+ ifactor));
            end;
        end;
        
end;



