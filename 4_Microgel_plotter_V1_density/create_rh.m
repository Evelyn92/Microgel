function hydro=create_rh(circ_x,circ_y,Rh)
hold on
th = 0:pi/50:2*pi;
xunit = Rh * cos(th) + circ_x;
yunit = Rh * sin (th) + circ_y;
hydro = plot(xunit,yunit,'-','LineWidth',5,'Color',[0 0 0]);
end
