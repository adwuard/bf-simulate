% BeamForming 3D simulation at fixed frequency f
% circular wavefront from point source on circumference with radius R
% plot surface to show intensity at given angle of arrival



% This script performs a 3D simulation: for every angle of arrival in the 3D space, for the
% specified frequency, it computes the intensity at a reference distance R. Set f as desired
% before calling the script.

if ~exist('config','var'), BFconfig; end;

N = 100; % subdivisions in the range, should be even
R = 8; % reference radius [m], intensity computed at this distance
[Xm,Ym,Zm]=sphere(N); Xm=Xm*R; Ym=Ym*R; Zm=Zm*R;
intM=zeros(N+1);
wT = 1/f; wl = wT*v; % wavecycle time interval [s] and wavelength [m]
mv = zeros(Nmic,1); % aux vect
for ix = 1:N+1,
    for iy = 1:N+1,
        x = Xm(ix,iy); y = Ym(ix,iy); z = Zm(ix,iy);
        % circular wavefront (accurate everywhere, including near field)
        for im = 1:Nmic,
            xd = m(im,1)-x; yd = m(im,2)-y; zd = m(im,3)-z;
            ph1 = (rem(mt(im),wT)/wT*2*pi); % phase delay from time delay
            ph2 = (rem(sqrt(xd*xd+yd*yd+zd*zd),wl)/wl*2*pi); % phase from distance
            mv(im) = ph1+ph2;
        end
        intM(ix,iy) = 20*log10(abs(sum(mw.*exp(1i.*mv)))); % max intensity in dB
    end
end
i = isfinite(intM(:)); intM=intM-min(intM(i)); intM(~i)=0;
dBmax=ceil(max(intM(:))/10)*10;
intMX = Xm.*intM/R; intMY = Ym.*intM/R; intMZ = Zm.*intM/R;
axv=linspace(-dBmax,+dBmax,2); [axM,~]=meshgrid(axv,axv); axM1=ones(length(axv));
i=round(linspace(0,1,256).^2*255)+1; % increase visibility of bluish tones (attenuation)
figure; hold on; cmap=parula(256); cmap=cmap(i,:); colormap(cmap);
h = surf(intMX,intMY,intMZ,intM+dBmax); axis equal; axis off
h.FaceAlpha = 0.6; h.FaceColor = 'interp'; h.EdgeColor = 'none';
h = patch([-dBmax -dBmax +dBmax +dBmax],[-dBmax +dBmax +dBmax -dBmax],[0 0 0 0],0);
h.FaceAlpha = 0.4; h.FaceColor = [0.85 0.85 0.85]; h.EdgeColor = 'none';
mygray = [0 0 0];
surf(axM,axM',-dBmax*axM1,0,'EdgeColor',mygray,'FaceColor','w');
surf(axM,+dBmax*axM1,axM',0,'EdgeColor',mygray,'FaceColor','w');
surf(+dBmax*axM1,axM,axM',0,'EdgeColor',mygray,'FaceColor','w');
N2 = N/2+1; mygray = [0.5 0.5 0.5];
plot3(intMX(N2,:),intMY(N2,:),intMZ(N2,:),'k-','LineWidth',2); % XY plane
plot3(intMX(N2,:),intMY(N2,:),intMZ(N2,:)-dBmax,'Color',mygray); % XY plane, projection
plot3(intMX(:,1),intMY(:,1),intMZ(:,1),'k:','LineWidth',2); % XZ plane, X<0
plot3(intMX(:,N2),intMY(:,N2),intMZ(:,N2),'k:','LineWidth',2); % XZ plane, X>0
plot3(intMX(:,1),intMY(:,1)+dBmax,intMZ(:,1),':','Color',mygray); % XZ plane, X<0, projection
plot3(intMX(:,N2),intMY(:,N2)+dBmax,intMZ(:,N2),':','Color',mygray); % XZ plane, X>0, projection
camlight left; lighting gouraud; view(3);
title(sprintf('%s, f=%.1fkHz',BFstr,f/1e3));