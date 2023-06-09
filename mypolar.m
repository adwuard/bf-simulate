function h = mypolar(a,r,sc) % angle vect, radius vect, scaling [rmin rmax]
rmin = min(sc); rmax = max(sc); rn = (r-rmin)/(rmax-rmin);
rn = max(0,rn); % avoid mirroring around center
h = plot(cos(a).*rn,sin(a).*rn); hold on;
end