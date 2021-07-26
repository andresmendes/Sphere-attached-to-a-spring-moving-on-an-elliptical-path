%% Sphere attached to a spring moving on an elliptical path 
% Simulation and animation of a sphere of mass m attached to a spring of
% constant k sliding without friction on a horizontal surface. 
%
% SAMPLE PROBLEM 13.8 - Beer, F.P., Johnston Jr, E.R., Mazurek, D.F.,
% Cornwell, P.J., 2013. Vector mechanics for engineers: Statics and
% Dynamics. McGraw-Hill.
%
%%

clear ; close all ; clc

% Parameters 
k = 100;                        % Spring constant               [N/m]
m = 0.6;                        % Sphere mass                   [kg]

playback_speed = 0.05;          % Speed of playback

% Parameters
tf      = 1.5;                  % Final time                    [s]
fR      = 30/playback_speed;    % Frame rate                    [fps]
dt      = 1/fR;                 % Time resolution               [s]
time    = linspace(0,tf,tf*fR); % Time                          [s]

vA      = 20;                   % Initial speed                 [m/s]
i_a     = 60*pi/180;            % Initial angle of speed vector [rad]

% Initial conditions
r0      = 0.5;                  % Initial radial position       [m]
th0     = 0;                    % Initial angular position      [rad]
dr0     = vA*cos(i_a);          % Initial dr/dt                 [m/s]
dth0    = vA*sin(i_a)/r0;       % Initial dth/dt                [rad/s]
Z0      = [ r0 th0 dr0 dth0];   % Initial conditions

% Simulation
options = odeset('RelTol',1e-5);
[TOUT,ZOUT] = ode45(@(t,z) dynamics(t,z,m,k),time,Z0,options);

% Retrieving states
r   = ZOUT(:,1);                % Radial position               [m]
th  = ZOUT(:,2);                % Angular position              [rad]
dr  = ZOUT(:,3);                % dr/dt                         [m/s]
dth = ZOUT(:,4);                % dth/dt                        [rad/s]

% Position
x = r.*cos(th);                 % x position [m]
y = r.*sin(th);                 % y position [m]

% Speed
v = sqrt(dr.^2 + (r.*dth).^2);  % [m/s]

% Energy
T = 1/2*m*v.^2;                 % Kinetic energy                [J]
V = 1/2*k*r.^2;                 % Potential energy              [J]

%% Results

figure
set(gcf,'Position',[50 50 1280 720]) % YouTube: 720p
% set(gcf,'Position',[50 50 854 480]) % YouTube: 480p
% set(gcf,'Position',[50 50 640 640]) % Social

% Create and open video writer object
v = VideoWriter('sphere_spring_elliptical.mp4','MPEG-4');
v.Quality = 100;
open(v);

for i=1:length(time)

    subplot(1,4,1:3)
        cla
        hold on ; grid on ; axis equal
        set(gca,'xlim',[1.2*min(x) 1.2*max(x)],'ylim',[1.1*min(y) 1.1*max(y)])
        plot(x,y,'k')
        plot([0 x(i)],[0 y(i)],'k:')
        plot(x(i),y(i),'ro','MarkerFaceColor','r')
        legend('Path','Spring','Sphere','Location','SouthEastOutside')
        xlabel('x position [m]')
        ylabel('y position [m]')
        title(strcat('m=',num2str(m),' kg ; k=',num2str(k),' N/m ; Time=',num2str(time(i),'%.3f'),' s (Playback speed=',num2str(playback_speed),')'))
    subplot(1,4,4)
        hold on ; grid on
        set(gca,'XTick',[],'ylim',[0 max(T+V)])
        cla
        Y = [T(i) ; V(i)];
        bar(1,Y,'stacked')
        ylabel('Energy [J]')
        legend('Kinetic','Potential','Location','NorthOutside')

    frame = getframe(gcf);
    writeVideo(v,frame);
    
end

close(v);

%% Auxiliary functions

function dz = dynamics(~,z,m,k)

    % States
    r   = z(1);
    % th  = z(2);
    dr  = z(3);
    dth = z(4);

    % Dynamics
    dz(1,1) = dr;
    dz(2,1) = dth;
    dz(3,1) = -(- m*r*dth^2 + k*r)/m;
    dz(4,1) = -(2*dr*dth)/r;

end
