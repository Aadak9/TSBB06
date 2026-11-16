clear;
clc;
close all;

%% Load images and camera matrices

im1 = imread('dinosaur0.png');
im2 = imread('dinosaur5.png');

load('dino_Ps','P');

C1 = P{1};
C2 = P{6};

%% Calculate camera centres

c1 = null(C1);
c2 = null(C2);

%% Calculate epipoles from the camera matrices

e12 = C1*c2;   % Camera 2 projected into image 1
e21 = C2*c1;   % Camera 1 projected into image 2

%% Calculate the fundamental matrix

F21 = liu_crossop(e21)*C2*pinv(C1);
F21 = F21/norm(F21,'fro');

F12 = F21';

%% Select corresponding points

[y1,y2] = correspondences_select(im1,im2);

% Convert to homogeneous coordinates
y1(end+1,:) = 1;
y2(end+1,:) = 1;

%% Calculate epipolar constraint values

num_of_points = size(y1,2);
epipolar_error = zeros(1,num_of_points);

for k = 1:num_of_points
    epipolar_error(k) = y2(:,k)'*F21*y1(:,k);
end

disp('Epipolar constraint values:');
disp(epipolar_error);

%% Calculate epipolar lines

l2 = F21*y1;
l1 = F21'*y2;

for k = 1:num_of_points
    l1(:,k) = -l1(:,k)/norm(l1(1:2,k))*sign(l1(3,k));
    l2(:,k) = -l2(:,k)/norm(l2(1:2,k))*sign(l2(3,k));
end

%% Plot points and epipolar lines

figure(2);
clf;

subplot(1,2,1);
image(im1);
axis image;
hold on;
title('Image 1');

subplot(1,2,2);
image(im2);
axis image;
hold on;
title('Image 2');

for k = 1:num_of_points

    subplot(1,2,1);
    plot(y1(1,k),y1(2,k),'ro','MarkerSize',8);
    drawline(l1(:,k));

    subplot(1,2,2);
    plot(y2(1,k),y2(2,k),'ro','MarkerSize',8);
    drawline(l2(:,k));
end

%% Calculate epipoles from the fundamental matrix

e12 = null(F21);
e21 = null(F21');

% Convert to Cartesian coordinates
e12 = e12/e12(3);
e21 = e21/e21(3);

fprintf('\nEpipole in image 1: (%.2f, %.2f)\n', ...
    e12(1),e12(2));

fprintf('Epipole in image 2: (%.2f, %.2f)\n', ...
    e21(1),e21(2));

%% Determine where the epipoles are located

[h1,w1,~] = size(im1);
[h2,w2,~] = size(im2);

inside1 = e12(1) >= 1 && e12(1) <= w1 && ...
          e12(2) >= 1 && e12(2) <= h1;

inside2 = e21(1) >= 1 && e21(1) <= w2 && ...
          e21(2) >= 1 && e21(2) <= h2;

if inside1
    fprintf('The first epipole is inside image 1.\n');
else
    fprintf('The first epipole is outside image 1.\n');
end

if inside2
    fprintf('The second epipole is inside image 2.\n');
else
    fprintf('The second epipole is outside image 2.\n');
end

% More detailed location for image 1
if e12(1) > w1
    fprintf('The first epipole is to the right of image 1.\n');
elseif e12(1) < 1
    fprintf('The first epipole is to the left of image 1.\n');
end

if e12(2) < 1
    fprintf('The first epipole is above image 1.\n');
elseif e12(2) > h1
    fprintf('The first epipole is below image 1.\n');
end

% More detailed location for image 2
if e21(1) > w2
    fprintf('The second epipole is to the right of image 2.\n');
elseif e21(1) < 1
    fprintf('The second epipole is to the left of image 2.\n');
end

if e21(2) < 1
    fprintf('The second epipole is above image 2.\n');
elseif e21(2) > h2
    fprintf('The second epipole is below image 2.\n');
end

%% Verify the epipoles

fprintf('Error norm(F21*e12): %e\n',norm(F21*e12));
fprintf('Error norm(F21''*e21): %e\n',norm(F21'*e21));

abs(sum(y1.*l1))
abs(sum(y2.*l2))
