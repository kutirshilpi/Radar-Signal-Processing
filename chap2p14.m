
D=4;
theta = 1:.01:2*pi;
c=3e8;
F=0.08e9;
lamda=c/F;
RCS=4*(abs(cos(pi*D*sin(theta)/lamda))).^2;

plot(theta,db(RCS/max(RCS)),'*-')

