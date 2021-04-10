function [xx,yy,zz,tt,RMS]=eqloc(indir, x0, y0, z0, t0, v, n, varargin)

% This function is written to estimate a location of an earthquake based on
% an initial guess and information of stations that have recorded the
% earthquake. It works for a one homogeneous velocity layer model.
%
% The file should be in a certain format
%
% Written by Huda Al Alawi - April 3rd, 2021
% Last modified by Huda Al Alawi - April 7th, 2021
%
% INPUT:
%
% indir         The directory of the file at which the input data is saved (including file name)
% x0, y0, z0    The location of the initial guess in regards to a reference point (included in the file)
% t0            Initial guess of origin time
% v             Velocity of the medium
% n             Number of iterations
% varargin      Replace it with 'n' or 'N' if no output files are needed
%               Replace it with 'y' or 'Y' if output files are needed and indicate the directory at which we want to save the output files (only the directory NOT including any names)
%
%
% EXAMPLE
% [x,y,z,t,RMS]=eqloc('/Users/hma/Desktop/HW8/Homework8/loc_dat_1', 0, 0, -10, 0, 6, 10, 'n')
% OR
% [x,y,z,t,RMS]=eqloc('/Users/hma/Desktop/HW8/Homework8/loc_dat_1', 0, 0, -10, 0, 6, 10, 'y', '/Users/hma/Desktop/')
%
%
% OUTPUT:
%
% It will only return the final estimated location. The rest will be plots and files.
% xx, yy, zz            The final estimated location
% tt                    The final estimated time
% RMS                   Last root-mean-square error of residuals
% modelsolutions.txt    Contains x0, y0, z0, and th0 and their uncertainties for all iterations 
% stationsinfo.txt      Contains the stations with their observed and calculated arrival time after obtaining the last solution besides the difference and RMS

% Save the first guess in a vector
guess1=[x0,y0,z0,t0];
% Open the file and load the data
fid=fopen(indir,'r');
% Skip the first line
fgetl(fid);
% Reda the reference now from the second line
% Will read the latitude and longitude which come after the name of the
% station (normally 3-4 characters)
ref=cell2mat(textscan(extractAfter(fgetl(fid),4),'%f%f'));
reflat=ref(1);
reflon=ref(2);

% Read the rest of the stations
table=textscan(fid,'%s%f%f%f%f%f%f');
sta=string(table{1});
lat=[table{2}];
lon=[table{3}];
elev=[table{4}]*0.001;
x=[table{5}];
y=[table{6}];
t_obs=[table{7}];
% Close the file
fclose(fid);

% Open file for iterations updates IF user choose 'y'
if varargin{1}=='y' || varargin{1}=='Y'
    iterf=strcat(varargin{2},'modelsolutions.txt');
    fid1=fopen(iterf,'a');
    % Print the reference info
    fprintf(fid1,'Reference Lat. & Long.      %.3f   %.3f\n\n',reflat,reflon);
end

% Plot the reference station, stations & the initial guess
subplot(1,2,1)
plot3(0,0,0,'o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',6)
hold on
grid on
plot3(x,y,elev,'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5)
plot3(x0,y0,z0,'o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',5)
title(sprintf('Visualization of Earthquake Location Changes \n Initial estimation (%.1f, %.1f, %.1f, %.1f)',x0,y0,z0,t0))
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')

% Repeate for a specific number of iterations
for ii=1:n
    % Find the distance between the stations and the current solution of m
    %(initially, it is the location determined by the user as initial guess)
    dd=sqrt((x(:)-x0).^2 + (y(:)-y0).^2 + (elev(:)-z0).^2);
    % The predicted arrival times is the refrence time + time based on the
    % calculated distances & given velocity
    t_pred=t0+(dd/v);
    
    % Save the first distance and time into matrix:
    if ii==1
        d1=dd(:);
    end
    
    % Calculating distance and time residuals
    dx=x(:)-x0;
    dy=y(:)-y0;
    dz=elev(:)-z0;
    dt=t_pred-t0;
    
    % For the model Am=d, A contains the partial derivatives, then
    dtdx= -dx./(v^2 * dt);
    dtdy= -dy./(v^2 * dt);
    dtdz= -dz./(v^2 * dt);
    dtdt=ones(length(sta),1);
    
    A=[dtdx,dtdy,dtdz,dtdt];
    d=t_obs-t_pred;
    % Now the solution m will be...
    cii=inv(A' * A);
    m=cii*A'*d;
    % Update the relative location
    x0=x0+m(1);
    y0=y0+m(2);
    z0=z0+m(3);
    t0=t0+m(4);
    % Find delta_m
    delta_m=sqrt(diag(cii)).*(sqrt((sum(d.^2)) ./ (length(sta)-4)));
    
    % Find RMS error for each iteration
    rmse(ii)=sum(d.^2)/length(sta);
    
    % Print the results into the file IF specified
    if varargin{1}=='y' || varargin{1}=='Y'
        fprintf(fid1,'Iteration#%d\n',ii);
        fprintf(fid1,'x0=      %.3f   %.3f\n',x0,delta_m(1));
        fprintf(fid1,'y0=      %.3f   %.3f\n',y0,delta_m(2));
        fprintf(fid1,'z0=      %.3f   %.3f\n',z0,delta_m(3));
        fprintf(fid1,'th0=      %.3f   %.3f\n',t0,delta_m(4));
        fprintf(fid1,'-----------------\n');
    end
    
    % Plot the second estimation
    if ii==1
        plot3(x0,y0,z0,'o','MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',5)
    end
end

if varargin{1}=='y' || varargin{1}=='Y'
    fclose(fid1);
end

% Plot the final estimation
plot3(x0,y0,z0,'o','MarkerEdgeColor','c','MarkerFaceColor','c','MarkerSize',5)
legend('Reference station','Stations','Initial estimation','Second estimation','Final estimation')
hold off

% Plot RMS error of residuals
subplot(1,2,2)
plot(rmse,'r','LineWidth',1)
title('RMS Error of Residuals')
xlabel('Iteration')
ylabel('RMS')
grid on

% Plot the fit of the data for the first and last solution
figure
% First
subplot(1,2,1)
plot(d1,t_obs,d1,d1/v,'LineWidth',1)
title('Observed Data and Theoretical Fit after the First Iteration')
xlabel('Distance [km]')
ylabel('Time [s]')
legend('Observed data','Theoretical fit')
grid on
% Last
subplot(1,2,2)
plot(dd,t_obs,dd,dd/v,'LineWidth',1)
title('Observed Data and Theoretical Fit after the Final Iteration')
xlabel('Distance [km]')
ylabel('Time [s]')
legend('Observed data','Theoretical fit')
grid on

% Visualization of how each of x,y,z,t has changed from first to last estimation
figure
X={'x';'y';'z';'t'};
vals=[guess1(1),x0 ; guess1(2),y0 ; guess1(3),z0 ; guess1(4),t0];
b=bar(vals);
set(gca, 'XTickLabel',X, 'XTick',1:numel(X))
grid on
title('Visualization of the Final Estimated x, y, z, t Values with Respect to the Initial Guess')
legend('Initial guess','Final solution')
yticks([round(min(min(vals))):1:round(max(max(vals)))]);



% Open a file to print stations info. IF specified
if varargin{1}=='y' || varargin{1}=='Y'
    stainfo=strcat(varargin{2},'stationsinfo.txt');
    fid2=fopen(stainfo,'a');
    % Print the header
    fprintf(fid2,'Station    x(km)    y(km)    z(km)    t_observed    t_calculated    Residual \n');
    % Print data
    for jj=1:length(sta)
        fprintf(fid2,'%-s %12.3f %9.3f %7.3f %9.3f %14.3f %15.3f \n',sta(jj),x(jj),y(jj),elev(jj),t_obs(jj),t_pred(jj),d(jj));
    end
    % Print last RMS error
    fprintf(fid2,'------ \n RMS of residuals = %.3f',rmse(end));
    fclose(fid2);
end


% Return the values
xx=x0; yy=y0; zz=z0; tt=t0; RMS=rmse(end);

end
