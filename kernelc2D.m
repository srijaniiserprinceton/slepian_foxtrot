function D=kernelc2D(XY,XYP,K)
% D=KERNELC2D(XY,XYP,K)
%
% Calculates the 2D Cartesian Slepian concentration (Bessel-1) kernel.
% SIMONS & WANG, doi: 10.1007/s13137-011-0016-z, eq. (54c)
%
% INPUT:
%
% XY      [X(:) Y(:)] A two-column matrix with one set of coordinate pairs 
% XYP     [XP(:) YP(:)] A two-column matrix with another set of points
% K       The wavenumber of the circular bandlimitation
%
% OUTPUT:
%
% D       The desired spatiospectral localization kernel, a matrix of
%         [length(XY(:))-by-length(XYP(:))]
%
% SEE ALSO:
%
% LOCALIZATION2D
%
% Last modified by dongwang-at-princeton.edu, 02/22/2008
% Last modified by fjsimons-at-alum.mit.edu, 07/27/2022

t0=clock;
% Make all the required combinations, i.e. the pairs of the unwrapped
% pairs of points; switch order so the dimensions are right
[XX,XXP]=meshgrid(XYP(:,1),XY(:,1));
[YY,YYP]=meshgrid(XYP(:,2),XY(:,2));

% Calculate the distance norm
md=sqrt((XX-XXP).^2+(YY-YYP).^2);

% Calculate the actual kernel, but watch out for zero division
% So - not actually a singular kernel, right?
% Remember in full 2D this is at fixed order-1, as opposed to when the
% radial functions are calculated at single orders m, see SWDISK
warning off MATLAB:divideByZero
D=K*besselj(1,K*md)/2/pi./md;
warning on MATLAB:divideByZero

% Supply the correct form of the kernel where the argument was zero. Note
% that this is not necessarily on the diagonal - the matrix D might not
% even be square
D(find(md<eps))=K^2/4/pi;

% Report on calculation time
disp(sprintf('KERNELC2D took %8.4f s',etime(clock,t0)))
