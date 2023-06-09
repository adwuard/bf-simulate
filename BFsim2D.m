% This script performs a 2D simulation: for every point in the 2D plane of the board, it
% simulates a circular and planar wavefront and computes the intensity of the beamforming
% output. Set the frequency f as desired before calling the script.

% BeamForming 2D simulation at fixed frequency f
% . circular wavefront from near source in 2D space
% . planar wavefront from infinitely far source in 2D space
if ~exist('config','var'), BFconfig; end;
N = 100; % subdivisions in the range
S = 8; % reference to define x,y range [m]
xrange = [-S,+S]; % [m]
yrange = [-S/2,+S/2]; % [m]
yrange = [-S,+S]; % [m]
xv = linspace(min(xrange),max(xrange),N);
yv = linspace(min(yrange),max(yrange),N);
intMc = zeros(N); % NxN intensity matrix for circular wavefront
intMp = zeros(N); % NxN intensity matrix for planar wavefront
wT = 1/f; wl = wT*v; % wavecycle time interval [s] and wavelength [m]
mv = zeros(Nmic,1); % aux vect
mb = mean(m,1); % barycenter of mic position
for ix = 1:N,
    for iy = 1:N,
        x = xv(ix); y = yv(iy);
        
        % circular wavefront (accurate everywhere, including near field)
        for im = 1:Nmic,
            xd = m(im,1)-x; yd = m(im,2)-y;
            ph1 = (rem(mt(im),wT)/wT*2*pi); % phase delay from time delay
            ph2 = (rem(sqrt(xd*xd+yd*yd),wl)/wl*2*pi); % phase from distance
            mv(im) = ph1+ph2; % here +pi is equivalent to weight*(-1)
        end
        intMc(iy,ix) = abs(sum(mw.*exp(1i.*mv))); % max intensity
        
        % planar wavefront (accurate only in far field, NOT in near field)
        x0 = m(1,1); y0 = m(1,2);
        mv(1) = (rem(mt(1),wT)/wT*2*pi); % 1st mic is reference
        for im = 2:Nmic,
            % line through x1,y1=barycenter and x2,y2=x,y has given orientation
            % line through x1,y1=n-th mic and x2,y2=x+xd,y+yd has same orientation
            xd = m(im,1)-mb(1); yd = m(im,2)-mb(2);
            x1 = m(im,1); y1 = m(im,2);
            x2 = x+xd; y2 = y+yd;
            A = y2-y1; B = x1-x2; %C = x2*y1-y2*x1; % coefficients for line equation, Ax+By+C=0
            Ap = -B; Bp = A; % coefficients for perpendicular line, dot=A*Ap+B*Bp=0
            C = -Ap*m(im,1)-Bp*m(im,2); % must pass through 2nd mic
            d = (Ap*x0 + Bp*y0 + C) / sqrt(Ap*Ap+Bp*Bp); % perpendicular distance, sign here matters!
            ph1 = (rem(mt(im),wT)/wT*2*pi); % phase from time delay
            ph2 = (rem(d,wl)/wl*2*pi); % phase from distance
            mv(im) = ph1+ph2; % here +pi is equivalent to weight*(-1)
        end
        intMp(iy,ix) = abs(sum(mw.*exp(1i.*mv))); % max intensity
    end
end
sc = max(max(intMc(:)),max(intMp(:))); % true max value
%sc = maxgain; % largest possible value
intMcp=intMc-intMp;
dsc=max(abs(intMcp(:))); % [-dsc +dsc] range so that 0 (no diff) is in the middle
figure; hold on;
imagesc(xv,yv,intMcp,[-dsc +dsc]); axis xy; axis equal; colorbar;
xlabel('x [m]'); ylabel('y [m]');
for im = 1:Nmic, plot(m(im,1),m(im,2),'ks','LineWidth',2); end;
plot(mb(1),mb(2),'kx','LineWidth',2);
title(sprintf('circ.-planar, %s, f=%.1fkHz',BFstr,f/1e3));
figure; hold on;
imagesc(xv,yv,intMp,[0 sc]); axis xy; axis equal; colorbar;
xlabel('x [m]'); ylabel('y [m]');
for im = 1:Nmic, plot(m(im,1),m(im,2),'ks','LineWidth',2); end;
plot(mb(1),mb(2),'kx','LineWidth',2);
title(sprintf('planar wave, %s, f=%.1fkHz',BFstr,f/1e3));
figure; hold on;
imagesc(xv,yv,intMc,[0 sc]); axis xy; axis equal; colorbar;
xlabel('x [m]'); ylabel('y [m]');
for im = 1:Nmic, plot(m(im,1),m(im,2),'ks','LineWidth',2); end;
title(sprintf('circ. wave, %s, f=%.1fkHz',BFstr,f/1e3));