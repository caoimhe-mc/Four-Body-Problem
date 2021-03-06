 function FourBodyProb_Visual(Storage,k)
% Creates a figure showing the location of the Earth, the Moon, Jupiter 
% and the Sun from the [x y z] coordinates in rows 1 to 12 in the kth column
% of the matrix Storage.  Also plots a trace for Earth, the Moon and
% Jupiter showing the 5000 coordinates of the bodies to theleft of the kth
% column
% 
% Inputs: 
% Storage (Matrix , real) At least 12 rows which represent four
% sets of [x y z] coordinates for the Earth, the Moon, Jupiter and the Sun.
% With the bodies coordinates stores as folows
% Earth - Storage(1:3,k)
% Moon - Storage(4:6,k)
% Jupiter - Storage(7:9,k)
% Sun - Storage(10:12,k)
% k (real scalar, positive) the column of Storage matrix that holds the
% coordinates of the four bodies that will be plotted.

% Version 1: Created 05/05/2021 Authors:    Harry Watchorn, 
%                                           Caoimhe McCann,             
%                                           Rachel Naughton
% Error Checking 
if (~ismatrix(Storage)) || size(Storage,1) < 12 || (~isreal(Storage))
    error('Input matrix Storage must contain only real values and must have a minimum of 12 rows.')
end

if (~isscalar(k)) || (~isreal(k)) || k <= 0 
    error('Input argument k must be a positive real scalar')
end

AU = 149597870.7; % astronomical unit in km
% Solar, Jovian, Earth and moon radius in AU
RS = 695700/AU; 
RJ = 69911/AU; 
RE = 6371/AU;
RM = 1737.1/AU;

scale = 100; 
RSs = scale*RS; % amplification of Solar size
scaleJ = 1000; 
RJs = scaleJ*RJ; % amplification of Jovian size
scaleE = 4000;
REs = scaleE*RE; % amplification of Earth size
scaleM = 4500;
RMs = scaleM*RM; % amplification of Moon size

S_im = imread('Sun.jpg'); % read in Sun texture
J_im = imread('Jupiter.jpg'); % read in Jupiter texture
M_im = imread('Moon.jpg'); % read in Moon texture
E_im = imread('Earth.jpg'); % read in Earth texture

[x,y,z] = sphere(50);
%Scaling the Moons coordinates so that tye earth and moon are not on top of
%eachother
newmoon = 60.*(Storage(4:6,k)-Storage(1:3,k))+Storage(1:3,k);
%Getting kth position of all four bodies
v_E = Storage(1:3,k);
v_M = newmoon(1:3);
v_J = Storage(7:9,k);
v_S = Storage(10:12,k);

axis([-12 12 -12 12 -12 12]) ;
set(gca,'color',[0 0 0 ]);
set(gcf,'color',[0 0 0 ]);
axis off
hold on
% adding traces of the orbits of Earth, Moon and Jupiter 
    if k <= 5000
        MoonTrace = 60.*(Storage(4:6,1:k)-Storage(1:3,1:k))+Storage(1:3,1:k);
        plot3(Storage(1,1:k)*2,Storage(2,1:k)*2,Storage(3,1:k)*2,'b')
        plot3(MoonTrace(1,1:k)*2,MoonTrace(2,1:k)*2,MoonTrace(3,1:k)*2,'w')
        plot3(Storage(7,1:k)*2,Storage(8,1:k)*2,Storage(9,1:k)*2,'y')
    else
        MoonTrace = 60.*(Storage(4:6,k-5000:k)-Storage(1:3,k-5000:k))+Storage(1:3,k-5000:k);
        plot3(Storage(1,k-5000:k)*2,Storage(2,k-5000:k)*2,Storage(3,k-5000:k)*2,'b')
        plot3(MoonTrace(1,1:5000)*2,MoonTrace(2,1:5000)*2,MoonTrace(3,1:5000)*2,'w')
        plot3(Storage(7,k-5000:k)*2,Storage(8,k-5000:k)*2,Storage(9,k-5000:k)*2,'y')
    end

%Plotting the bodies of the Sun, Jupiter, Moon and Earth
Sun = surf(RSs*x+v_S(1,1)*2,RSs*y+v_S(2,1)*2,RSs*z+v_S(3,1)*2,'facecolor','texturemap','cdata' ...
    ,S_im,'edgecolor','none');
Jupiter = surf(RJs*x+v_J(1,1)*2,RJs*y+v_J(2,1)*2,RJs*z+v_J(3,1)*2,'facecolor','texturemap','cdata'...
    ,J_im,'edgecolor','none');
Moon = surf(RMs*x+v_M(1,1)*2,RMs*y+v_M(2,1)*2,RMs*z+v_M(3,1)*2,'facecolor','texturemap','cdata'...
    ,M_im,'edgecolor','none');
Earth = surf(REs*x+v_E(1,1)*2,REs*y+v_E(2,1)*2,REs*z+v_E(3,1)*2,'facecolor','texturemap','cdata'...
    ,E_im,'edgecolor','none');

set(Sun,'AmbientStrength',1);
set(Jupiter,'AmbientStrength',0.4);
set(Earth,'AmbientStrength',0.4);
set(Moon,'AmbientStrength',0.6);
light('Style','local','Position',[v_S(1,1) v_S(2,1) v_S(3,1)]);
hold off
end
