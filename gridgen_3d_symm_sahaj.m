clear all;
close all;
%% get BL thickness
nu=1/1000;    %% KINEMATIC VISCOSITY
lref=1.0;
uref=1.0;
Re=uref*lref/nu;
blt=5*lref/sqrt(Re);
dx=blt/10;
%% start!
box_x1=0;
box_x2=2;
box_y1=-0.6;
box_y2=-box_y1;
upstream_x=-30;
downstream_x=17;
side_y=30;
div_fact=8; %%(ensures number of grid points is divisible by this number)

us_str_p=10;
ds_str_p=5;
sd_str_p=9;

Re=lref*uref/nu;
local_length=5/sqrt(Re); % (or refer to pointwise y+ calculator)
delta_box=local_length/10;
delta_box=1/32;
%% construct upstream
x_us=[];
us_rat=(100+us_str_p)/100;
delta_current=delta_box*us_rat;
x_us(1)=box_x1-delta_current;
pres_p=x_us(1);
i=2;
while pres_p>upstream_x
    delta_current=delta_current*us_rat;
    x_us(i)=x_us(i-1)-delta_current;
    pres_p=x_us(i);
    i=i+1;
end
delta_us=delta_current;
%% x_box
x_box=[];
x_box(1)=box_x1;
pres_p=x_box(1);
i=2;
while pres_p<box_x2
    delta_current=delta_box;
    x_box(i)=x_box(i-1)+delta_box;
    pres_p=x_box(i);
    i=i+1;
end
delta_box_x=delta_box;
i_box_max=i-1;

%% construct downstream
x_ds=[];
ds_rat=(100+ds_str_p)/100;
delta_current=delta_box*ds_rat;
x_ds(1)=x_box(i_box_max)+delta_current;
pres_p=x_ds(1);
i=2;
while pres_p<downstream_x
    delta_current=delta_current*ds_rat;
    x_ds(i)=x_ds(i-1)+delta_current;
    pres_p=x_ds(i);
    i=i+1;
end
delta_ds=delta_current;

%% x assembly
x_grid=[flip(x_us) x_box x_ds];
%check divisibility by 32
N_x=length(x_grid);
check =0;
usds=0;
iter=0;
dsch=0;
pointsleft=div_fact-mod(N_x-1,div_fact);
pa_us=pointsleft/2;
%% add two points downstream for 1 point upstream
div_fact_x=div_fact;
if mod(N_x-1,div_fact_x)>0.5
    while check<0.5
        iter=iter+1;
        if mod(N_x,div_fact_x)<0.5
            check=1.0;
        end
        if usds<0.5 %% add point to upstream
            delta_ds=delta_ds*ds_rat;
            x_usp=x_grid(N_x)+delta_ds;
            x_grid=[x_grid x_usp];
            x_ds=[x_ds x_usp];
            N_x=length(x_grid);
            dsch=dsch+1;
            if dsch>2.5
                dsch=0.0;
                usds=1.0;
            end
        elseif usds>0.5 %% add point to upstream
            delta_us=delta_us*us_rat;
            x_dsp=x_grid(1)-delta_us;
            x_grid=[x_dsp x_grid];
            N_x=length(x_grid);
            x_us=[x_us x_dsp];
            usds=0.0;
        end
    end
end
%% upstream compress
us_length=abs(upstream_x-box_x1);
ain=delta_box;
n=length(x_us);
fun=@(r) ain*(r.^n - 1)/(r-1)-us_length;
r0=us_rat;
r_us=fzero(fun,r0);
x_us(1)=box_x1-delta_box;
h_us=delta_box*r_us;
for i=2:n
    x_us(i)=x_us(i-1)-h_us;
    h_us=h_us*r_us;
end
x_us(end)=upstream_x;
x_grid=[flip(x_us) x_box x_ds];
%% downstream compress
ds_length=abs(downstream_x-box_x2);
ain=delta_box;
n=length(x_ds);
fun=@(r) ain*(r.^n - 1)/(r-1)-ds_length;
r0=ds_rat;
r_ds=fzero(fun,r0);
x_ds(1)=x_box(end)+delta_box;
h_ds=delta_box*r_ds;
for i=2:n
    x_ds(i)=x_ds(i-1)+h_ds;
    h_ds=h_ds*r_ds;
end
x_ds(end)=downstream_x;
x_grid=[flip(x_us) x_box x_ds];
%% Construct y box: y_box should have odd number of elements
% y_box=[];
y_mid=(box_y1+box_y2)*0.5;
y_box=box_y1:delta_box:box_y2;
N_y_box=length(y_box);
if mod(N_y_box,2)<0.5
    y_box=linspace(box_y1,box_y2,N_y_box+1);
else
    y_box=linspace(box_y1,box_y2,N_y_box);
end
delta_box_y=abs(y_box(2)-y_box(1));
j_box_max=length(y_box);
%% construct sides
y_sd=[];
sd_rat=(100+sd_str_p)/100;
delta_current=delta_box_y*sd_rat;
y_sd(1)=y_box(j_box_max)+delta_current;
pres_p=y_sd(1);
j=2;
while pres_p<side_y
    delta_current=delta_current*sd_rat;
    y_sd(j)=y_sd(j-1)+delta_current;
    pres_p=y_sd(j);
    j=j+1;
end
delta_sd=delta_current;
%% y assembly
y_grid=[flip(-y_sd) y_box y_sd];
N_y=length(y_grid);
check =0;
usds=0;
iter=0;
pointsleft=div_fact-mod(N_x-1,div_fact);
pa_us=pointsleft/2;
while mod(N_y-1,div_fact)>0.5
        iter=iter+1;
        delta_sd=delta_sd*sd_rat;
        y_sd(length(y_sd)+1)=y_sd(length(y_sd))+delta_sd;
        y_grid=[flip(-y_sd) y_box y_sd];
        N_y=length(y_grid);
end
%% side compress
sd_length=abs(side_y-box_y2);
ain=delta_box;
n=length(y_sd);
fun=@(r) ain*(r.^n - 1)/(r-1)-sd_length;
r0=sd_rat;
r_sd=fzero(fun,r0);
y_sd(1)=box_y2+delta_box;
h_sd=delta_box*r_sd;
for i=2:n
    y_sd(i)=y_sd(i-1)+h_sd;
    h_sd=h_sd*r_sd;
end
y_sd(end)=side_y;
y_grid=[flip(-y_sd) y_box y_sd];
%% set min
xmin=-min(x_grid);
ymin=-min(y_grid);
x_grid=x_grid+xmin;
y_grid=y_grid+ymin;

%% plot grids
x_box_length=abs(x_box(end)-x_box(1));
y_box_length=abs(y_box(end)-y_box(1));
x_zi_l=xmin+x_box(1)-0.5*x_box_length;
x_zi_r=xmin+x_box(end)+0.5*x_box_length;
y_zi_up=y_box(end)+0.5*y_box_length;
[xg,yg]=meshgrid(x_grid,y_grid);
zg=xg.*0+yg.*0;
s=mesh(xg,yg,zg);
s.LineStyle='-';
s.EdgeColor='[0 0 0]';
s.FaceColor='[1.0000    0.9922    0.8157]';
view(2);
xlim([x_grid(1) x_grid(N_x)]);
ylim([y_grid(1) y_grid(N_y)]);
pbaspect([x_grid(N_x)-x_grid(1), y_grid(N_y)-y_grid(1), 1.0]);
xrect=x_box(1)+xmin;
yrect=y_box(1)+ymin;
xrect_center=x_box(1)+xmin +x_box_length/2;
yrect_center=y_box(1)+ymin+y_box_length/2;
rectangle('Position',[x_box(1)+xmin y_box(1)+ymin x_box_length y_box_length],'EdgeColor','r','LineWidth',2)
xlabel("x");
ylabel("y");
box on;
ax = gca;

exportgraphics(ax,"zoomedout_2dgrid_side.png",'Resolution',900);

%% zoomed in plot

xlim([x_zi_l x_zi_r]);
ylim([-y_zi_up+ymin y_zi_up+ymin]);
box on;
radius=0.5;
xcirc=xrect+0.6-radius;
ycirc=yrect_center-radius;
xcirc_center=xcirc+radius;
ycirc_center=ycirc+radius;
pos = [xcirc ycirc radius*2 radius*2]; 
rectangle('Position',pos,'Curvature',[1 1],'EdgeColor','b','LineWidth',2)
pbaspect([-x_zi_l+x_zi_r, 2*y_zi_up, 1.0]);
ax = gca;

exportgraphics(ax,"zoomedin_2dgrid_side.png",'Resolution',900);

%% write grid information to file

fileID = fopen('3d_gridinfo.txt','w');
fprintf(fileID,'3D Sim Info \r\n');
fprintf(fileID,'This is a grid generated for the following case: \r\n');
fprintf(fileID,'Lref= %8.3f \r\n',lref);
fprintf(fileID,'Uref= %8.3f  \r\n',uref);
fprintf(fileID,'Re=1/nu= %8.3f \r\n',1/nu);
fprintf(fileID,'Relevant Physical Length Scale chosen as BL thickness is, l_p= 5/sqrt(Re)= %8.3f Lref \r\n',local_length);
Ngrid=N_x*N_y*N_y;
str=CommaFormat(Ngrid);
fprintf(fileID,"Total number of grid points, N_x *N_y *N_y="+str);
fprintf(fileID,'\r\n \r\n');


fprintf(fileID,'3D Grid Info \r\n');
fprintf(fileID,'Direction: x \r\n');
fprintf(fileID,'N_x= %d \r\n', N_x);
fprintf(fileID,'Uniform Box Region: Centered at  %4.2f, Starts at %4.2f , ends at %4.2f  and has Delta x = l_p / %8.3f = %8.3f  \r\n', xmin,x_box(1)+xmin,x_box(end)+xmin,local_length/delta_box_x,delta_box_x );
fprintf(fileID,'Upstream Region: Starts at %4.2f , ends at %4.2f  and has stretch factor = %4.2f  \r\n', x_grid(1),x_box(1)+xmin,r_us );
fprintf(fileID,'Downstream Region: Starts at %4.2f , ends at %4.2f  and has stretch factor = %4.2f  \r\n \r\n', x_box(end)+xmin,x_grid(end),r_ds );

fprintf(fileID,'3D Grid Info \r\n');
fprintf(fileID,'Direction: y and z \r\n');
fprintf(fileID,'N_y= %d \r\n', N_y);

fprintf(fileID,'Uniform Box Region: Centered at %4.2f ,Starts at %4.2f , ends at %4.2f  and has Delta y = l_p / %8.3f = %8.3f  \r\n', ymin,ymin+y_box(1),ymin+y_box(end),local_length/delta_box_y,delta_box_y );
fprintf(fileID,'Down Region: Starts at %4.2f , ends at %4.2f  and has stretch factor = %4.2f  \r\n', y_grid(1),ymin+y_box(1),r_sd );
fprintf(fileID,'Up Region: Starts at %4.2f , ends at %4.2f  and has stretch factor = %4.2f  \r\n \r\n', ymin+y_box(end),y_grid(end),r_sd );

fprintf(fileID,'TimeStep size selection \r\n');
fprintf(fileID,'Smallest grid spacing is Delta= x = %8.3f \r\n', min(delta_box_y,delta_box_x) );
fprintf(fileID,'For Convective CFL = 1 (C-CFL=uref*dt/dx+2*uref*dt/dy), Delta t= %8.3f \r\n', 1.0/(uref/delta_box_x+2*uref/delta_box_y));
fprintf(fileID,'For Diffusive CFL = 1 (D-CFL=2*nu*dt/dx^2+2*2*nu*dt/dy^2), Delta t= %8.3f \r\n', 1.0/(2*nu/(delta_box_x*delta_box_x)+2*2*nu/(delta_box_y*delta_box_y)));
fprintf(fileID,'Cell Peclet Number (Pe=uref*dx/nu+2*uref*dy/nu) is Pe= %8.3f \r\n', uref*delta_box_x/nu+2*uref*delta_box_y/nu);
fclose(fileID);
fclose('all')
%% Write xgrid
fxg = fopen('xgrid.dat','w');
for i=1:length(x_grid)
    fprintf(fxg,'%d %15.10f \n',i,x_grid(i));
end
%% Write ygrid
fxg = fopen('ygrid.dat','w');
for i=1:length(y_grid)
    fprintf(fxg,'%d %15.10f \n',i,y_grid(i));
end
%% Write zgrid
fxg = fopen('zgrid.dat','w');
for i=1:length(y_grid)
    fprintf(fxg,'%d %15.10f \n',i,y_grid(i));
end
fclose('all')


%% pretty num
%=====================================================================
% Takes a function and inserts commas for the thousands separators.
function [commaFormattedString] = CommaFormat(value)
	try
		% Split into integer part and fractional part.
		[integerPart, decimalPart]=strtok(num2str(value),'.'); 
		% Reverse the integer-part string.
		integerPart=integerPart(end:-1:1); 
		% Insert commas every third entry.
		integerPart=[sscanf(integerPart,'%c',[3,inf])' ... 
				repmat(',',ceil(length(integerPart)/3),1)]'; 
		integerPart=integerPart(:)'; 
		% Strip off any trailing commas.
		integerPart=deblank(integerPart(1:(end-1)));
		% Piece the integer part and fractional part back together again.
		commaFormattedString = [integerPart(end:-1:1) decimalPart];
	catch ME
		errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
			ME.stack(1).name, ME.stack(1).line, ME.message);
		set(handles.txtInfo, 'String', errorMessage);
		WarnUser(errorMessage);
	end
	return;
end
