%%  sliding_animation1
%   Description:
%       This script uses the sliding dynamics and tries to visualize their
%       animation.

clear all; close all; clc;

%% Add Functions to Path
addpath(genpath('../../src/'))

%% Create System
random_seed = 31;
rng(random_seed)
sys1 = FrictionfulSliding();

%% Create figure for two different positions

figure;
subplot(2,1,1)
sys1.show(zeros(2,1));

subplot(2,1,2)
sys1.show([0.5,1.0])

%% Create Simulation of Simple Version
sys1.theta = sampleFromPolytope(sys1.Theta);
g = 10;

x0 = sampleFromPolytope(sys1.X0);
x0(2) = abs(x0(2));
u0 = [8.6;0.0];

N_sim = 250;
dt = 0.01;

N_turnoff = 60;

% Simulate
x_history = [x0];
for k = [0:N_sim-1]
    %Update before turnoff
    x_k = x_history(:,k+1);

    if k < N_turnoff
        x_kp1 = x_k + dt * sys1.dynamics(x_k,u0);
    else
        x_kp1 = x_k + dt * sys1.dynamics(x_k,zeros(2,1));
    end

    x_history = [x_history, x_kp1];


end

figure;
hold on;
plot(x_history(1,:),x_history(2,:))
scatter(x_history(1,1),x_history(2,1),'filled','o','LineWidth',2.0)

%% Create Animation

if ~exist('FileName')
	FileName = 'sliding_test1.mp4';
end
vidObj = VideoWriter(FileName,'MPEG-4');

open(vidObj);

% Plotting
figure;
for k = 0:N_sim
    % Plot
    x_t = x_history(:,k+1);
    if k < N_turnoff
        hs = sys1.show(x_t, 'magenta');
    else
        hs = sys1.show(x_t);
    end

    % Get current frame and write it.
    currFrame = getframe;
    writeVideo(vidObj,currFrame);
    delete(hs);

    % Prepare for plot to be overwritten
    hold off;
end

close(vidObj);