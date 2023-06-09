function mypolarsetup(at,rt,fmt) % angle and radius ticks
if isempty(at), at = linspace(0,330,12)*pi/180; end
if isempty(rt), rt = linspace(0,1,6); end
N = 100; a = linspace(0,2*pi,N); ax = pi*75/180;
figure; axis off; axis equal; axis([-1 1 -1 1]); hold on;
rtmax = max(rt); rtmin = min(rt); rtn = (rt-rtmin)/(rtmax-rtmin);
p = patch(cos(a),sin(a),[1 1 1]); p.EdgeColor = 'none';
for i = 1:length(rt),
    line(cos(a)*rtn(i),sin(a)*rtn(i),'Color',[0.85 0.85 0.85]);
    text(cos(ax)*rtn(i),sin(ax)*rtn(i),sprintf(fmt,rt(i)));
end
for i = 1:length(at),
    line([0 cos(at(i))],[0 sin(at(i))],'Color',[0.85 0.85 0.85]);
    t = text(1.1*cos(at(i)),1.1*sin(at(i)),sprintf('%3.0f',at(i)*180/pi));
    t.HorizontalAlignment = 'center';
    t.VerticalAlignment = 'middle';
end
end