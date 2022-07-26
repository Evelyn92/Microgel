function hydro=createsc_rh(circ_x,circ_y,Rh)
hold on
th = 0:pi/50:2*pi;
xunit = Rh * cos(th) + circ_x;
yunit = Rh * sin (th) + circ_y;
zunit = 0;
for i=1:length(xunit)-1
    zunit=[zunit 0];
end

hydro = plot3(xunit,yunit,zunit,'-','LineWidth',5,'Color',[0 0 0]);
end
