% This scripts is the same as BFsimAF.m except that the plot is a ring instead of a plane. The
% radius corresponds to a given frequency, the angle corresponds to the angle of arrival.


% BeamForming 2D simulation
% circular wavefront from point source on circumference with radius R
% plot ring to show intensity at given angle of arrival and frequency
if ~exist('config','var'), BFconfig; end;
N = 5*100; % subdivisions in the range
S = max(abs(m(:,1)))*10; % reference to define x,y range [m]
xrange = [-S,+S]; % [m]
yrange = [-S,+S]; % [m]
xv = linspace(min(xrange),max(xrange),N);
yv = linspace(min(yrange),max(yrange),N);
mb=mean(m,1); % barycenter of mic position
d = sqrt(m(:,1).^2+m(:,2).^2); % distance origin-mic
XR = [xrange; 0 0]; YR = [0 0; yrange];
DR = sqrt(XR.^2+YR.^2); % distance origin-middlepoints
R1 = 1.05*max(d); % min radius correspond to f1
R2 = 0.95*min(DR(:)); % max radius correspond to f2
R = 0.5; % reference radius [m], intensity computed at this distance
intM = NaN(N); intA = zeros(N); % intensity and alpha matrix
mv = zeros(Nmic,1); % aux vect
for ix = 1:N,
    for iy = 1:N,
        x = xv(ix); y = yv(iy); % x,y -> freq,angle
        d = sqrt(x*x+y*y);
        if d<R1, continue; end;
        if d>R2, continue; end;
        ft = f1 + (f2-f1)*(d-R1)/(R2-R1); wTt = 1/ft; wlt = wTt*v; % wavelength [m]
        a = atan2(y,x); % angle -> x,y
        x = R*cos(a); y = R*sin(a);
        
        % circular wavefront (accurate everywhere, including near field)
        for im = 1:Nmic,
            xd = m(im,1)-x; yd = m(im,2)-y;
            ph1 = (rem(mt(im),wTt)/wTt*2*pi); % phase delay from time delay
            ph2 = (rem(sqrt(xd*xd+yd*yd),wlt)/wlt*2*pi); % phase from distance
            mv(im) = ph1+ph2; % here +pi is equivalent to weight*(-1)
        end
        intM(iy,ix) = abs(sum(mw.*exp(1i.*mv))); % max intensity
        intA(iy,ix) = 1;
    end
end
%sc = max(intM(:)); % true max value
sc = maxgain; % largest possible value
%intM=intM/sc; % normalize
figure; hold on;
h = imagesc(xv,yv,intM); set(h,'AlphaData',intA); axis xy; axis equal; colorbar;
for r = linspace(R1,R2,6),
    h = polar(linspace(0,2*pi,50),r*ones(1,50),'b-'); set(h,'Color',[0.25 0.25 0.25]);
end
contour(xv,yv,intM,6,'Color',[1 1 1]);
for r = linspace(R1,R2,6),
    rf = f1 + (f2-f1)*(r-R1)/(R2-R1); if r<R2, sunit=''; else sunit='kHz'; end
    text(r*cos(pi/4),r*sin(pi/4),sprintf('\n%.1f%s',rf/1e3,sunit));
end
h = gca; set(h,'XTick',[]); set(h,'YTick',[]); set(h,'YColor',[1 1 1]); set(h,'XColor',[1 1 1])
for im = 1:Nmic, plot(m(im,1),m(im,2),'ks','LineWidth',2); end;
for a = +pi : -pi/6 : -pi+0.001,
    if cos(a)>0, hstr='left'; else hstr='right'; end; if abs(cos(a))<0.01, hstr='center'; end;
    if sin(a)>0, vstr='bottom'; else vstr='top'; end; if abs(sin(a))<0.01, vstr='middle'; end;
    h=text('Position',[R2*cos(a),R2*sin(a),0],'String',sprintf('%.0f',a*180/pi), ...
        'HorizontalAlignment',hstr,'VerticalAlignment',vstr);
end
title(sprintf('%s, R=%.1fcm',BFstr,R*100));