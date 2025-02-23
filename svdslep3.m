function varargout=svdslep3(XY,KXY,J,tol,ngro,xver)
% [E,V,c11cmnR,c11cmnK,SE,KXY]=SVDSLEP3(XY,KXY,J,tol,ngro,xver)
%
% Two-dimensional Slepian functions with arbitrary concentration/limiting
% regions in the Cartesian spatial and (half-)spectral domains.
%
% INPUT:
%
% XY       [X(:) Y(:)] coordinates of a SPATIAL-domain curve
% KXY      [X(:) Y(:)] coordinates of a SPECTRAL-domain half-curve 
%          i.e. in the positive (including zero) halfplane!
%          For various reasons the "curve" is not the boundary but rather
%          the last set of "elements" on the "grid".
% J        Number of eigentapers [default: 10]
% tol      abs(log10(tolerance)) for EIGS [default: 12]
% ngro     The computational "growth factor" [default: 8]
% xver     Performs excessive verification [default: 0]

%
% OUTPUT:
%
% E        The eigenfunctions of the concentration problem
% V        The eigenvalues of the concentration problem
% c11cmnR  The spatial coordinates of the top left corner
% c11cmnK  The spectral coordinates of the top left corner
% SE       The periodogram of the eigenfunctions
% KXY      The symmetrized spectral-space domain
%
% EXAMPLE:
%
% svdslep3('demo1')
% 
% SEE ALSO:
%
% LOCALIZATION2D
%
% Last modified by fjsimons-at-alum.mit.edu, 07/28/2022

% Default values
defval('J', 10);
defval('ngro',3);
defval('xver',1);
defval('tol',12);

% Default curve is a CIRCLE in pixel space, of some radius and pixelization
defval('cR',30)
defval('cN',41)
defval('XY',...
       cR*[cos(linspace(0,2*pi,cN)) ; sin(linspace(0,2*pi,cN))]')
% And some (half-)square in the spectral (half-)space
% Here you use the notion of the Shannon ratio as in SVDSLEP2, which is
% relative to the grown area and has unit values in the CORNERS
defval('R',0.1)
defval('KXY',...
       R*[-1 1  1 -1 -1; 1 1  0  0  1]')

if ~isstr(XY)
  % Check the curves and return the range on the inside 
  % For the SPATIAL part
  [QinR,xylimt,QX,QY]=ccheck(XY);
  % For the SPECTRAL part - including the mirrored symmetry
  [QinK,kxylimt,QXK,QYK,KXY]=ccheck(KXY,[],1);
  
  % Check if you must
  if xver==1
    clf
    % Plot the SPATIAL region of interest exactly as input
    subplot(221)
    plot(XY(:,1),XY(:,2),'b'); hold on
    % Compute the centroid - my function is from SAGA, but R2022a has its own
    [X0,Y0]=centroid(XY(:,1),XY(:,2));
    plot(X0,Y0,'b+')
    axis equal
    xlim(minmax(XY(:))*1.5)
    ylim(minmax(XY(:))*1.5)

    % Plot the SPECTRAL region of interest exactly as input
    subplot(222)
    plot(KXY(:,1),KXY(:,2),'r'); hold on
    % Compute the centroid - my function is from SAGA, but R2022a has its own
    [KX0,KY0]=centroid(KXY(:,1),KXY(:,2));
    plot(KX0,KY0,'r+')
    axis equal
    axis([-1 1 -1 1]/sqrt(2))

keyboard

    subplot(223)
    dom=zeros(size(QX)); dom(QinR)=1; spy(dom)
    axis ij
    subplot(224)
    dom=zeros(size(QXK));
    [Kn,a,b,dci,dcn]=knum2(size(dom));
    imagesc(Kn); axis image ; hold on
    dom(QinK)=1; spy(dom,'kx')
    % Check that this is Hermitian
    difer(isreal(ifft2(ifftshift(dom)))-1,[],[],NaN)
    pause(2)
  end

  % Now embed these in a larger size index matrix to get rid of the edge
  % effects 
  newsize=ngro*max(size(QX),size(QXK));

  [QinR,c11cmnR]=cexpand(QinR,QX,QY,newsize);
  
  [QinK,c11cmnK]=cexpand(QinK,QXK,QYK,newsize);
    
  % But now this needs to be turned into a FFTSHIFT
  QinK=indeks(fftshift(v2s(1:prod(newsize))),QinK);

  if xver==1
    clf
    subplot(121)
    dom=zeros(newsize); dom(QinR)=1; spy(dom)
    axis ij
    aa=axis;
    subplot(122)
    dom=zeros(newsize); dom(QinK)=1; spy(fftshift(dom))
  end

  % Now make the operators that we are trying to diagonalize
  P=@(x) proj(x,QinR);
  % We're finding VECTORS that are going to be 2-D images!
  Q= @(x) fft2vec(x);
  Qi=@(y) ifft2vec(y);
  L=@(y) proj(y,QinK);
  H=@(x) P(Qi(L(Q(P(x)))));

  % And then find the eigenvectors and eigenvalues
  OPTS.isreal=false;
  OPTS.disp=0;
  defval('tolerance',10^-tol);
  OPTS.tol=tolerance;
  OPTS.maxit=500;

  % Remember to specify the output size
  [E,V]=eigs(H,prod(newsize),J,'LR',OPTS);
  
  [V,i]=sort(diag(V),'descend');
  E=E(:,i); V=V(1:J); E=E(:,1:J);

  % Define some kind of tolerance level
  tol=sqrt(eps); 

  % Make them real as we know they should be
  if any(imag(V)>tol)
    error('Complex eigenvalues');
  else
    V=real(V);
    % Check imaginary part of the "good" eigenfunctions
    disp(sprintf('mean(abs(imag(E))) = %8.3e out of %8.3e\n',...
		 mean(mean(abs(imag(E(:,V>tol))))),...
		 mean(mean(abs(E(:,V>tol))))))
    % Note that they were normalized in the complex plane
    E=real(E); E=E./repmat(diag(sqrt(E'*E))',size(E,1),1);
  end

  if nargout>4
    % Get the power spectrum
    SE=zeros(prod(newsize),size(E,2));
    for i=1:size(E,2),
      SE(:,i)=indeks((abs(fft2(v2s(E(:,i)))).^2),':');
    end
  else
    SE=NaN;
  end

  % Output
  varns={E,V,c11cmnR,c11cmnK,SE,KXY};
  varargout=varns(1:nargout);
elseif strcmp(XY,'demo1')
  defval('R',300);
  % A circle in space
  XY=R*[cos(linspace(0,2*pi,circn)) ; sin(linspace(0,2*pi,circn))]';
  % A half triangle in spectral space
  KXY=sqrt(pi*(R)^2)/2*[0  1/2 -1/2 0;...
    		          0  1  1 0]';
  % Clockwise rotation i the Fourier domain!
  rr=rotz(pi/5);
  KXY=[rr(1:2,1:2)*KXY']';
    
  % How many eigenfunctions?
  J=40;
  % Compute the eigenfunctions
  [E,V,c11cmnR,c11cmnK,SE,KXY]=svdslep3(XY,KXY,J);
  
  % Quick fix
  c11=c11cmnR(1:2);
  cmn=c11cmnR(3:4);

  % Make the figures
  offs=0;

  figure(1)
  clf
  [ah,ha]=krijetem(subnum(2,3));
  for ind=1:3
    axes(ah(ind))
    % This is in actual units
    imagefnan(c11,cmn,v2s(E(:,ind+offs)))
    hold on
    plot(XY(:,1),XY(:,2),'k'); hold off
    title(sprintf('%s = %i','\alpha',ind+offs))
    
    axes(ha(2*ind))
    psdens=fftshift(decibel(v2s(SE(:,ind+offs))));
    psdens(psdens<-80)=NaN;
    imagefnan(c11,cmn,psdens);
    hold on
    plot(KXY(:,1),KXY(:,2),'k'); hold off
  end
  fig2print(gcf,'landscape')

  % Also try this one here
  figure(2)
  clf
  subplot(121)
  EE=sum(repmat(V(:)',length(E),1).*E.^2,2);
  SEE=sum(repmat(V(:)',length(E),1).*SE.^2,2);
  imagefnan(c11,cmn,v2s(EE)); axis image 
  hold on
  plot(XY(:,1),XY(:,2),'k'); hold off
  subplot(122)
  psdens=fftshift(decibel(v2s(SEE)));
  psdens(psdens<-80)=NaN;
  imagefnan(c11,cmn,psdens); axis image
  axis([c11(1) cmn(1) cmn(2) c11(2)]/4)
  hold on
  plot(KXY(:,1),KXY(:,2),'k'); hold off

  % Also try this one here
  figure(3)
  clf
  plot(V,'o-')
  title(sprintf('sum of the eigenvalues %8.3f',sum(V)))
  longticks(gca,2)
  ylim([-0.1 1.1])
  grid on
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Qin,xylimt,QX,QY,XY]=ccheck(XY,dxdy,isk)
% Given coordinates of a closed curve XY, computes

% The following was stripped from LOCALIZATION2D 

% In the spectral domain the input curve needs to be doubled
defval('isk',0)

% Make sure the XY of the curve has two columns
if ~isempty(find(size(XY,2)==2))
  if size(XY,1)<size(XY,2)
    XY=XY';
  end
else
  error('Coordinates of the curve not as expected')
end
keyboard
% Find limits in x and y so as to contain the curves
if isk==0
  % We're in x-space
  xlimt=minmax(XY(:,1)); 
  ylimt=minmax(XY(:,2)); 
elseif isk==1
  % We're in k-space, symmetrize the half curve
  xlimt=[-1 1]*max(abs(XY(:,1)));
  % Note that Y must be >= 0
  ylimt=[-1 1]*max(XY(:,2));
  % Mirror the curve itself
  XY=[XY ; flipud(-XY)];
end
xrange=xlimt(2)-xlimt(1);
pyrange=ylimt(2)-ylimt(1);

% Do this for the benefit of the output
xylimt=[xlimt ylimt];

% Construct a grid with the region inscribed
qdx=dxdy(1); qdy=dxdy(2);

% Make it square for now
Nx=round(xrange/qdx);
Ny=round(yrange/qdy);
[Nx,Ny]=deal(max(Nx,Ny));

% Don't be a fanatic but strive for symmetry - compare LOCALIZATION2D
qx=linspace(xlimt(1),xlimt(2),Nx);
qy=linspace(ylimt(1),ylimt(2),Ny);
  
% Remember that this may not be square
[QX,QY]=meshgrid(qx,qy);

% Should look into this later 
QY=flipud(QY);

% The midpoint indices of this subset that fall inside of the region...
Qin=find(inpolygon(QX,QY,XY(:,1),XY(:,2)));

xver=0;
if xver==1
  clf
  plot(QX,QY,'k.'); hold on 
  plot(QX(Qin),QY(Qin),'o'); 
  twoplot(XY)
  hold off
  axis image
  pause(2)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [newQin,c11cmn]=cexpand(Qin,QX,QY,newsize)
oldsize=size(QX); 
addon=round([newsize-oldsize]/2);
[i,j]=ind2sub(oldsize,Qin);
newQin=sub2ind(newsize,i+addon(1),j+addon(2));
addx=range(QX(1,:))/size(QX,2)*addon(2);
addy=range(QY(:,1))/size(QY,1)*addon(1);
c11=[min(QX(1,:)) max(QY(:,1))]+[-addx addy];
cmn=[max(QX(1,:)) min(QY(:,1))]+[addx -addy];
c11cmn=[c11 cmn];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Pv=proj(v,p)
% Projects the vector v on the indices p
Pv=zeros(size(v));
Pv(p)=v(p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Fv=fft2vec(v)
% Returns the two-dimensional FFT of a vector
Fv=fft2(v2s(v));
Fv=Fv(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function iFv=ifft2vec(Fv)
% Returns the two-dimensional IFFT of a vector
iFv=ifft2(v2s(Fv));
iFv=iFv(:);

