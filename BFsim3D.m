
% This script performs a 3D simulation: for every point in the 3D space, it simulates a circular
% wavefront and computes the intensity of the beamforming output. Set the frequency f as
% desired before calling the script.

% BeamForming 3D simulation at fixed frequency f
% circular wavefront from point source in 3D space
if ~exist('config','var'), BFconfig; end;
N = 100; % subdivisions in the range
S = 8; % reference to define x,y range [m]
xrange = [-S,+S]; % [m]
yrange = [-S,+S]; % [m]
zrange = [-S,+S]; % [m]
xv = linspace(min(xrange),max(xrange),N);
yv = linspace(min(yrange),max(yrange),N);
zv = linspace(min(zrange),max(zrange),N);
intM = zeros(N,N,N); % NxNxN intensity matrix for circular wavefront
wT = 1/f; wl = wT*v; % wavecycle time interval [s] and wavelength [m]
mv = zeros(Nmic,1); % aux vect
for ix = 1:N,
    for iy = 1:N,
        for iz = 1:N,
            x = xv(ix); y = yv(iy); z = zv(iz);
            % circular wavefront (accurate everywhere, including near field)
            for im = 1:Nmic,
                xd = m(im,1)-x; yd = m(im,2)-y; zd = m(im,3)-z;
                ph1 = (rem(mt(im),wT)/wT*2*pi); % phase delay from time delay
                ph2 = (rem(sqrt(xd*xd+yd*yd+zd*zd),wl)/wl*2*pi); % phase from distance
                mv(im) = ph1+ph2;
            end
            intM(iy,ix,iz) = abs(sum(mw.*exp(1i.*mv))); % max intensity
        end
    end
end
alpha=0.9;
[XV,YV,ZV]=meshgrid(xv,yv,zv);
cmax = max(intM(:)); cmin = min(intM(:)); % true range
%cmax = maxgain; cmin = 0; % largest possible range
c1 = cmin+0.9*(cmax-cmin); % isosurface to show captured directions
c2 = cmin+0.1*(cmax-cmin); % isosurface to show rejected directions
figure; hold on;
%hy = slice(XV,YV,ZV,intM,[],0,[]);
%hy.FaceAlpha = 0.5; hy.EdgeColor = 'none';
hz = slice(XV,YV,ZV,intM,[],[],0);
hz.EdgeColor = 'none';
set(gca,'CLim',[cmin cmax]);
p1 = patch(isosurface(XV,YV,ZV,intM,c1));
isonormals(XV,YV,ZV,intM,p1)
p1.FaceAlpha = alpha; p1.FaceColor = 'yellow'; p1.EdgeColor = 'none';
p1c = patch(isocaps(XV,YV,ZV,intM,c1,'above'));
p1c.FaceAlpha = alpha; p1c.FaceColor = 'yellow'; p1c.EdgeColor = 'none';
if abs(c2-c1)>0.01,
    p2 = patch(isosurface(XV,YV,ZV,intM,c2));
    isonormals(XV,YV,ZV,intM,p2)
    p2.FaceAlpha = alpha; p2.FaceColor = 'blue'; p2.EdgeColor = 'none';
    p2c = patch(isocaps(XV,YV,ZV,intM,c2,'below'));
    p2c.FaceAlpha = alpha; p2c.FaceColor = 'blue'; p2c.EdgeColor = 'none';
end
camlight left
lighting gouraud
view(3); axis equal; axis([-S +S -S +S -S +S]); grid on;
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
for im = 1:Nmic, plot3(m(im,1),m(im,2),m(im,3),'ks','LineWidth',2); end;
title(sprintf('%s, f=%.1fkHz',BFstr,f/1e3));