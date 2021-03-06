function F = FourBodyProblem_Equations(vec)
% F = FourBodyProblem_Equations(vec)
% This vector function is a Four-body problem involving Earth, Moon, Jupiter and
% the Sun. 
% Newton's Law of gravitation: 
%     F= -G(mass_i*mass_j/|rij^2|)r
% is used to calculate the Acceleration of Four bodies: 
%     A1= -G*{(mass2*(r1-r2)/|r12|^3)+(mass3*(r1-r3)/|r13|^3)+(mass4*(r1-r4)/|r14|^3)}
%
% Input:
% vec (real vector) is a 24 x 1 column vector, containing the initial positions of the
% planets in x, y and z coordinates followed by the initial velocities in x, y and z. 
% The order of input must be as follows:
% Earth, Moon, Jupiter, Sun positions in x,y then z 
% Earth, Moon, Jupiter, Sun velocities in x,y then z
% distance units is AU, time units is days
%
% Output: 
% F (real vector) is a 24 x 1 column vector that contains the velocity and acceleration of the
% bodies in x,y,z coordinates. The order of output data is given as:
% Earth, Moon, Jupiter, Sun velocities
% followed by Earth, Moon, Jupiter, Sun accelerations in x,y and then z.
 
% Version 1: Created 31/04/2021 Authors:  Caoimhe McCann, 
%                                         Harry Watchorn, 
%                                         Rachel Naughton 

% Error checking   
if (~isreal(vec) || ~iscolumn(vec))
    error('Input vector must be real values in a column vector')
end
    
G = 1.4879e-34;     %define the gravitational constant
MJ = 1.898e27;      %define the mass of jupiter
MS = 1.98855e30;    %define the mass of the sun
ME = 5.97219e24;    %define the mass of earth
MM = 7.34767309e22; %define the mass of the moon
    
%calculate the distance between each pair of planets
rerm = sqrt((vec(1,1)-vec(4,1))^2 + (vec(2,1)-vec(5,1))^2 +(vec(3,1)-vec(6,1))^2);
rerj = sqrt((vec(1,1)-vec(7,1))^2 + (vec(2,1)-vec(8,1))^2 + (vec(3,1)-vec(9,1))^2);
rers = sqrt((vec(1,1)-vec(10,1))^2 + (vec(2,1)-vec(11,1))^2 + (vec(3,1)-vec(12,1))^2);
rmrj = sqrt((vec(4,1)-vec(7,1))^2 + (vec(5,1)-vec(8,1))^2 + (vec(6,1)-vec(9,1))^2);
rmrs = sqrt((vec(4,1)-vec(10,1))^2 + (vec(5,1)-vec(11,1))^2 + (vec(6,1)-vec(12,1))^2);
rjrs = sqrt((vec(7,1)-vec(10,1))^2 + (vec(8,1)-vec(11,1))^2 + (vec(9,1)-vec(12,1))^2);

F = zeros(24,1); % set aside storage space

F(1:3,1) = vec(13:15,1);    %earth's x,y and z velocity
F(4:6,1) = vec(16:18,1);    %moon's x,y and z velocity
F(7:9,1) = vec(19:21,1);    %jupiter's x,y and z velocity
F(10:12,1) = vec(22:24,1);  %sun's x,y and z velocity

%calculate the acceleration for the x,y and z component for every body
F(13:15,1) = -G*(MM*((vec(1:3,1)-vec(4:6,1))/(abs(rerm)^3)) ...
+ MJ*((vec(1:3,1)-vec(7:9,1))/(abs(rerj)^3)) + MS*(vec(1:3,1)-vec(10:12,1))/(abs(rers)^3));
F(16:18,1) = -G*(ME*((vec(4:6,1)-vec(1:3,1))/(abs(rerm)^3)) ...
+ MJ*((vec(4:6,1)-vec(7:9,1))/(abs(rmrj)^3)) + MS*(vec(4:6,1)-vec(10:12,1))/(abs(rmrs)^3));
F(19:21,1) = -G*(MM*((vec(7:9,1)-vec(4:6,1))/(abs(rmrj)^3)) ...
+ ME*((vec(7:9,1)-vec(1:3,1))/(abs(rerj)^3)) + MS*(vec(7:9,1)-vec(10:12,1))/(abs(rjrs)^3));
F(22:24,1) = -G*(ME*((vec(10:12,1)-vec(1:3,1))/(abs(rers)^3)) ...
+ MM*((vec(10:12,1)-vec(4:6,1))/(abs(rmrs)^3)) + MJ*(vec(10:12,1)-vec(7:9,1))/(abs(rjrs)^3));    
end
