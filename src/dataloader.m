clear all;
close all;

u=load("u.dat");
us=load("us.dat");
u0=load("u0.dat");

v=load("v.dat");
vs=load("vs.dat");
v0=load("v0.dat");

RHSu=load("RHSu.dat");
RHSv=load("RHSv.dat");

uD1=load("Du.dat");
vD1=load("Dv.dat");

ci=load("ci.dat");
cj=load("cj.dat");

cis=load("cis.dat");
cjs=load("cjs.dat");

f=load("f.dat");

massd=load("massd.dat");
massnew=load("massdnew.dat");

x=load("x.dat");
y=load("y.dat");
xc=load("xc.dat");
yc=load("yc.dat");

dx=load("dx.dat");
dy=load("dy.dat");
delx=load("delx.dat");
dely=load("dely.dat");

pgx=load("xgradp.dat");
pgy=load("ygradp.dat");

cleft=cis(1,2:end);
cright=cis(end,2:end);