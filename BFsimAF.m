% This script performs a 2D simulation: for every angle of arrival in the 2D plane of the board,
% and for all frequencies in the specified range, it computes the intensity at a reference
% distance R. This script uses two auxiliary scripts: mypolarsetup.m and mypolar.m. Set f1
% and f2 as desired before calling the script.

% BeamForming 2D simulation
% circular wavefront from point source on circumference with radius R
% plot surface to show intensity at given angle of arrival and frequency
% . extract frequency response H(f) at given angle of arrival
% . extract response as function of angle H(a) at given frequency
if ~exist('config','var'), BFconfig; end;
R = .5; % reference radius [m], intensity computed at this distance
avect = linspace(-pi,+pi,360); % angle vector, direction of arrival [rad]
fvect = linspace(f1,f2,200); % frequency vector [Hz]
NA = length(avect);
NF = length(fvect);
intM = zeros(NF,NA); % intensity matrix
mv = zeros(Nmic,1); % aux vect
for i = 1:NF,
    for ia = 1:NA,
        ft = fvect(i); wTt = 1/ft; wlt = wTt*v; % wavelength [m]
        a = avect(ia); x = R*cos(a); y = R*sin(a);
        
        % circular wavefront (accurate everywhere, including near field)
        for im = 1:Nmic,
            xd = m(im,1)-x; 
            yd = m(im,2)-y;
            ph1 = (rem(mt(im),wTt)/wTt*2*pi); % phase delay from time delay
            ph2 = (rem(sqrt(xd*xd+yd*yd),wlt)/wlt*2*pi); % phase from distance
            mv(im) = ph1+ph2; % here +pi is equivalent to weight*(-1)
        end
        intM(i,ia) = abs(sum(mw.*exp(1i.*mv))); % max intensity
    end
end
%sc = max(intM(:)); % true max value
sc = maxgain; % largest possible value
%intM=intM/sc; % normalize
intMdB = 20*log10(intM);
dBmax = ceil(max(intMdB(:))/10)*10;
dBmin = floor(min(intMdB(:))/10)*10;
dBnr = round((dBmax-dBmin)/20)+1;
% 3D plot (copy figure bitmap)
figure; hold on;
h = surf(avect*180/pi,fvect/1e3,intM);
h.FaceAlpha = 0.9; h.FaceColor = 'interp'; h.EdgeColor = 'none';
camlight left; lighting gouraud
contour3(avect*180/pi,fvect/1e3,intM);
grid on; zoom on; view(3);
xlabel('angle [deg]'); ylabel('freq [kHz]'); zlabel('intensity');
title(sprintf('%s, R=%.1fcm',BFstr,R*100));
% 2D plot
figure; hold on;
imagesc(avect*180/pi,fvect/1e3,intM); axis xy; axis tight; colorbar;
contour(avect*180/pi,fvect/1e3,intM,'Color',[1 1 1]);
xlabel('angle [deg]'); ylabel('freq [kHz]'); zlabel('intensity');
title(sprintf('%s, R=%.1fcm',BFstr,R*100));
% log plot of frequency response at given angle
a2vect = [0, 45, 90, 135, 180]*pi/180;
L = length(a2vect); cmap = copper(L); labels = {};
figure;
for i = 1:L,
    j = find(avect>=a2vect(i));
    if ~isempty(j),
        j=j(1);
        semilogx(fvect/1e3,20*log10(intM(:,j)),'Color',cmap(i,:),'LineWidth',2); hold on;
        labels{length(labels)+1} = sprintf('%.0fdeg',a2vect(i)*180/pi);
    end
end
legend(labels,'Location','southeast');
xlabel('freq [kHz]'); ylabel('intensity [dB]');
title(sprintf('%s, R=%.1fcm',BFstr,R*100));
axis([f1/1e3 f2/1e3 dBmin dBmax]); grid on; zoom on;
% polar plot of angle response at given frequency
% dBmin = -90; dBmax = +10; dBnr = 6; % override dB range
f2vect = [0.5 1 3 5 10 20]*1e3;
L = length(f2vect); cmap = copper(L); labels = {}; handles = [];
mypolarsetup([],linspace(dBmin,dBmax,dBnr),'%.0fdB');
for i = 1:L,
    j = find(fvect>=f2vect(i));
    if ~isempty(j),
        j=j(1);
        h = mypolar(avect,intMdB(j,:),[dBmin dBmax]);
        h.Color = cmap(i,:); h.LineWidth = 2;
        handles(length(handles)+1) = h;
        labels{length(labels)+1} = sprintf('%.1fkHz',f2vect(i)/1e3);
    end
end
legend(handles,labels,'Location','southeast');
title(sprintf('%s, R=%.1fcm',BFstr,R*100));