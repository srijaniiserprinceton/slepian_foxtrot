function varargout=svdslep4(yes1,no1,yes2,NW,offs)
% [E,V,W,ngro]=SVDSLEP4(yes1,no1,yes2,NW,offs)
%
% Given a target time window, compute a Slepian window, similar to
% SVDSLEP1 but for arbitrary time and frequency supports
%
% INPUT:
%
% yes1       A number of samples to include, followed by:
% no1        A number of samples to exclude, followed by:
% yes2       A number of samples to include
% NW         Half the Shannon number (time-bandwidth product) [default: 5]
% offs       Frequency offset in case you are building a sweep [default: 15]
% 
% OUPUT:
%
% E          The obtained Slepian time window
% V          The obtained Slepian eigenvalue
% W          The target time window
% NW         Half the Shannon number (time-bandwidth product) [default: 5]
% offs       Frequency offset in case you are building a sweep [default: 15]
% ngro       The growth factor used
%
% EXAMPLE:
%
% svdslep4('demo1')
%
% SEE ALSO:
%
% SURFACEWIN, PARTITA, PARTITF, PARTITS
%
% Last modified by fjsimons-at-alum.mit.edu, 08/17/2020

% I shouldn't forget to look into the Dirichlet distribution... 
% Change-making problems, partition problems, string cutting, etc....

% Default values
defval('Nmax',100)
defval('yes1',randi(Nmax))
defval('no1',randi(Nmax))
defval('yes2',randi(Nmax))
% Implicit method... see SVDSLEP1... need a growth factor for stability
defval('ngro',8);
% How about a subset instead so that you can sweep through frequencies
defval('NW',5);
defval('offs',25);
  
% The number of eigenvalues sought
defval('K',1)

if ~strcmp(yes1,'demo1')
  % Not a demo - defaults don't add up to anything special
  % Make the target time window
  W=[ones(1,yes1) zeros(1,no1) ones(1,yes2)]';
  % Now work with all the frequencies and find the Slepian window
  N=length(W);
  Nd=N*ngro;  
  nonz=(Nd-N)/2+[1:N];
  desir=zeros(1,Nd);
  desir(nonz)=W;
  
  % The time-projection
  P =@(x) proj(x,find(desir));
  % The Fourier operator
  Q =@(x) fft(x);
  % The inverse Fourier operator
  Qi=@(x) ifft(x);
  % The frequency projection
  nonz=[ngro*offs:ngro*offs+ngro*NW+1 Nd-ngro*NW+1-ngro*offs:Nd-ngro*offs];
  % This would show what we are doing
  % stem(proj(ones(1,Nd),nonz))
  % The below would be ALL the frequencies and defeat the point
  % nonz=[1:Nd];

  % The frequency projection
  L =@(y) proj(y,nonz);
  % The composite operator
  H =@(x) P(Qi(L(Q(P(x)))));
  % Acknowledge that H is complex (though it is symmetric)
  OPTS.isreal=false;
  OPTS.disp=0;
  
  % How many do you compute?
  K=100;
  % Remember to specify the output size
  [E,V]=eigs(H,Nd,K,'LR',OPTS);
  % Sorting
  if K>1
    [V,i]=sort(diag(V),'descend');
    E=E(:,i); V=V(1:K); E=E(:,1:K);
  end
  % Take out only the central part
  E=E((Nd-N)/2+1:(Nd-N)/2+N,:);
  % Note that they were normalized in the complex plane
  E=real(E); E=E./repmat(diag(sqrt(E'*E))',size(E,1),1);
else
  % The first input was 'demo1' so the second can be reused for pix

  % Default values
  defval('Nmax',100)
  % h=randfixedsumint(1,3,Nmax);
  h=randis(1,3,Nmax)';
  h=[400 200 100]';
  % Calculate the Slepian window
  [E,V,W,NW,offs]=svdslep4(h(1),h(2),h(3));

  % Just pick one out already for illustration
  defval('no1',[]); pix=no1;
  defval('pix',randi(size(E,2)));
  E=E(:,pix);
  V=V(pix);

  % Calculate the power spectral densities
  SW=abs(fft(W-mean(W)).^2);
  SE=abs(fft(E-mean(E)).^2);

  % Make the frequency axis for the positive frequencies only
  [fax,selekt]=fftaxis1D(E,length(SE),1);
  SE=SE(selekt);
  SW=SW(selekt);

  % Make the plot
  clf
  [ah,ha,H]=krijetem(subnum(2,2));

  axes(ah(1))
  p(1)=plot(W);
  ylim([-0.1 1.1])
  grid on
  xl(1)=xlabel('sample');
  yl(1)=ylabel('ideal time window');
  xlim([0 length(W)-1])

  axes(ah(2))
  p(2)=plot(fax,decibel(SW));
  grid on
  xl(2)=xlabel('frequency');
  yl(2)=ylabel('ideal spectral window');
  xlim([1 length(SW)-1])
  ylim([-50 1])
  
  axes(ah(3))
  % Norm adjustment is purely for visual pleasure
  p(3)=plot(E*sqrt(norm(W)));
  ylim([-1.1 1.1])
  grid on
  xl(3)=xlabel('sample');
  yl(3)=ylabel(sprintf('Slepian time window %i',pix));
  xlim([0 length(E)-1])
    
  axes(ah(4))
  p(4)=plot(fax,decibel(SE));
  grid on
  xl(4)=xlabel('frequency');
  yl(4)=ylabel(sprintf('Slepian spectral window %i',pix));
  xlim([1 length(SE)-1])
  ylim([-50 1])
  hold on
  pp=plot([offs offs+NW; offs offs+NW],[ylim; ylim]');
  set(pp,'Color','k')
  hold off
  bottom(p,ah(4))
  l=legend(sprintf('%s = %8.6f','\lambda',V),'Location','NorthEast');
  
  % Final cosmetics
  set(p(1:2),'Color','b')
  if V>0.5
    set(p(3:4),'Color','r')
  else
    set(p(3),'Color','r')
    set(p(4),'Color',grey)
  end
  
  set(ah,'FontSize',9)
  set(yl,'FontSize',9)
  set(xl,'FontSize',9)
  delete(xl(1:2))
  nolabels(ah(1:2),1)
  serre(H',1/2,'down')
  longticks(ah)
  
  % Print command
  figdisp([],sprintf('%2.2i',pix),[],2)
end

% Optional output
varns={E,V,W,NW,offs,ngro};
varargout=varns(1:nargout);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Pv=proj(v,p)
% Projects the vector v on the indices p
Pv=zeros(size(v));
Pv(p)=v(p);
