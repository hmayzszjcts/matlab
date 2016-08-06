
function draw(x,y,r,points)

ang=0:0.01:2*pi; 
xp=r*cos(ang);
yp=r*sin(ang);
plot(x+xp,y+yp);
axis equal
points_mat=points;
hold on
text(x,y,'*') 
text(points_mat(:,1),points_mat(:,2),'o') 



end