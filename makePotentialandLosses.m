clear
close all
clc

Ly = 256;
My = 256;
Lx = 256;
Mx = 256;
dx = Lx/Mx;
dy = Ly/My;
x = -Lx/2:dx:Lx/2-dx;
y = -Ly/2:dy:Ly/2-dy;
[X,Y] = meshgrid(x,y);

dkx = 2*pi/Lx; dky = 2*pi/Ly;
kxmax = 2*pi/dx; kymax = 2*pi/dy;
kx = -kxmax/2:dkx:kxmax/2-dkx;
ky = -kymax/2:dky:kymax/2-dky;
[Kx,Ky] = meshgrid(kx,ky);


R0 = 80;
njet = 2;
angles = 0:2*pi/njet:2*pi*(1-1/njet)
V =  X.^2 + Y.^2 > R0.^2 ;

test = zeros(size(V));
gamma = zeros(size(V));
loss = zeros(size(V));
for jj = 1:length(angles)
    
    test = test + (X > 0.9*R0 & X < R0+40 & abs(Y) < 5);
    Z = (X+1i*Y).*exp(-1i*angles(2));
    X = real(Z); Y = imag(Z);
    if floor(jj/2) == jj/2
        gamma =  gamma + 2*(X.^2 + Y.^2 > (R0+30).^2).*(abs(angle(Z)) < angles(2)/2);
        loss =  loss + 2.5*(X.^2 + Y.^2 > (R0+30).^2).*(abs(angle(Z)) < angles(2)/2);
    else
        gamma = gamma + 1*(X.^2 + Y.^2 > (R0+30).^2).*(abs(angle(Z)) < angles(2)/2);
        loss = loss +5*(X.^2 + Y.^2 > (R0+30).^2).*(abs(angle(Z)) < angles(2)/2);
    end
end
V = V & ~ test;

% [X,Y] = meshgrid(x,y);
% V = V - 0.5*(X.^2 + Y.^2 < R0.^2 & X.^2 + Y.^2 > (R0-5).^2);
% for jj = 1:length(angles)
%     Z = (X+1i*Y).*exp(-1i*angles(2));
%     X = real(Z); Y = imag(Z);
%     if floor(jj/2) == jj/2
%        disp(num2str(jj))
%         V = V + 0.5*(X > R0-5 & X < R0 & abs(Y) <5);
%     end
%     
% end


figure
subplot(121)
imagesc(x,y,V)
set(gca,'PlotBoxAspectRatio',[Lx/Ly 1 1])


V = fftshift(fft2(V));

V = V.*exp(-(X.^2+Y.^2)/40^2);
V = real(ifft2(ifftshift(V)));

subplot(122)
imagesc(x,y,V)
set(gca,'PlotBoxAspectRatio',[Lx/Ly 1 1])
hold on
plot(R0*exp(1i*angles),'sqr')

V = 3*(V);
hdf5write('Potential.h5','/1/V',V,'/1/x',x,'/1/y',y);
hdf5write('Losses.h5','/1/V2',gamma,'/1/loss',loss,'/1/x',x,'/1/y',y);

Vx = real(ifft2(ifftshift(1i*Kx.*fftshift(fft2(V)))));
Vy = real(ifft2(ifftshift(1i*Ky.*fftshift(fft2(V)))));


figure
imagesc(gamma)
%% Look at the Ground State

phiR = hdf5read('Groundstate.h5','/1/phiR');
phiI = hdf5read('Groundstate.h5','/1/phiI');
wavefunction = phiR + 1i*phiI;
figure
imagesc(x,x,abs(wavefunction).^2)
set(gca,'PlotBoxAspectRatio',[Lx/Ly 1 1])


%%
addpath ~/Documents/2DQT_Analysis_v1.0
addpath ~/Documents/2DQT_Analysis_v1.0/RCA_v3/
v = VideoWriter('peaks.avi');
open(v);
B = bwboundaries(V<0.25);
B = B{1};

RP = [];
RM = []
for ii =175:220
    
    loadfile = [num2str(ii) '.h5'];
    phiR = hdf5read(loadfile,'/1/phiR');
    phiI = hdf5read(loadfile,'/1/phiI');
    wavefunction = phiR + 1i*phiI;
    
    figure(2)
    clf
    subplot(131)
    imagesc(x,x,abs(wavefunction).^2)
    set(gca,'PlotBoxAspectRatio',[Lx/Ly 1 1],'ydir','normal')
    Fx(ii) = sum(abs(wavefunction(:)).^2.*Vx(:))*dx*dy;
    Fy(ii) = sum(abs(wavefunction(:)).^2.*Vy(:))*dx*dy;
    colorbar
    subplot(132)
    imagesc(x,x,angle(wavefunction))
    set(gca,'PlotBoxAspectRatio',[Lx/Ly 1 1],'ydir','normal')
    [n,nplus,nminus,xi,yi,kappai] = GetVortexPositions3(wavefunction,Lx,Ly,V < 0.25);
    subplot(133)
    plot(xi(1:nplus),yi(1:nplus),'.r')
    rp = [xi(1:nplus) yi(1:nplus)];
    rm = [xi(nplus+1:end) yi(nplus+1:end)];
    
    RP = [RP; rp];
    RM = [RM; rm];
    NP(ii) = nplus;
    NM(ii) = nminus;
%           if ~isempty(rp) && ~isempty(rm)
%           [PCLUSTERS,NCLUSTERS,DIPOLES] = RCA(rp,rm,Lx,Ly,1);
%          end
    xlim([-Lx/2 Lx/2])
    ylim([-Ly/2 Ly/2])
    set(gca,'PlotBoxAspectRatio',[Lx/Ly 1 1],'ydir','normal')
    box on
    hold on
    plot(xi(nplus+1:end),yi(nplus+1:end),'.b')
    plot(B(:,2)-Lx/2-1,B(:,1)-Ly/2-1,'-k')
    frame = getframe(gcf);
    writeVideo(v,frame);
    
    
end
close(v);


phase = angle(wavefunction);
[Ux,~] = gradient(unwrap(phase,[],2));
[~,Uy] = gradient(unwrap(phase,[],1));


