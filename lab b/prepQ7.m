clear;
clc;
close all;

%% Load images and camera matrices

im1 = imread('dinosaur0.png');
im2 = imread('dinosaur5.png');

load('dino_Ps', 'P');

C1 = P{1};
C2 = P{6};

[h, w, ~] = size(im1);

%% Calculate camera centres

c1 = null(C1);
c2 = null(C2);

%% Calculate epipoles

e12 = C1*c2;
e21 = C2*c1;

e12 = e12/e12(3);
e21 = e21/e21(3);

%% Calculate the fundamental matrices

F21 = liu_crossop(e21)*C2*pinv(C1);
F21 = F21/norm(F21, 'fro');

F12 = F21';

%% Verify the epipoles

disp('F21*e12 should be approximately zero:');
disp(F21*e12);

disp('F21''*e21 should be approximately zero:');
disp(F21'*e21);

%% Calculate Loop-Zhang distortion

lambda_values = -5:0.01:5;

dist = LoopZhangDistortion( ...
    e12, e21, w, h, lambda_values);

[minimum_distortion, index] = min(dist);
lambda = lambda_values(index);

fprintf('Optimal lambda: %f\n', lambda);
fprintf('Minimum distortion: %f\n', minimum_distortion);

%% Approximate derivative

dist_derivative = diff(dist)./diff(lambda_values);

lambda_derivative = ...
    (lambda_values(1:end-1) + lambda_values(2:end))/2;

%% Plot distortion and derivative

figure;

subplot(2,1,1);
plot(lambda_values, dist);
set(gca, 'YScale', 'log');
grid on;
xlabel('\lambda');
ylabel('Distortion');
title('Loop-Zhang distortion');

subplot(2,1,2);
plot(lambda_derivative, dist_derivative);
grid on;
xlabel('\lambda');
ylabel('Approximate derivative');
title('Derivative of Loop-Zhang distortion');

%% Construct the projective homographies

z = [lambda; 1; 0];

w1 = liu_crossop(e12)*z;
w1 = w1/w1(3);

w2 = F21*z;
w2 = w2/w2(3);

Hp1 = [1 0 0;
       0 1 0;
       w1'];

Hp2 = [1 0 0;
       0 1 0;
       w2'];

%% Construct the similarity homographies

vcp = 0;

Hr1 = [F21(3,2)-w1(2)*F21(3,3), ...
       w1(1)*F21(3,3)-F21(3,1), ...
       0;

       F21(3,1)-w1(1)*F21(3,3), ...
       F21(3,2)-w1(2)*F21(3,3), ...
       F21(3,3)+vcp;

       0 0 1];

Hr2 = [w2(2)*F21(3,3)-F21(2,3), ...
       F21(1,3)-w2(1)*F21(3,3), ...
       0;

       w2(1)*F21(3,3)-F21(1,3), ...
       w2(2)*F21(3,3)-F21(2,3), ...
       vcp;

       0 0 1];

%% Complete rectifying homographies

H1 = Hr1*Hp1;
H2 = Hr2*Hp2;

disp('Hp1 =');
disp(Hp1);

disp('Hp2 =');
disp(Hp2);

disp('Hr1 =');
disp(Hr1);

disp('Hr2 =');
disp(Hr2);

disp('H1 =');
disp(H1);

disp('H2 =');
disp(H2);

%% Verify that the epipoles map to horizontal infinity

e1_rectified = H1*e12;
e2_rectified = H2*e21;

e1_rectified = e1_rectified/e1_rectified(1);
e2_rectified = e2_rectified/e2_rectified(1);

disp('Rectified epipole 1:');
disp(e1_rectified);

disp('Rectified epipole 2:');
disp(e2_rectified);

%% Verify the transformed fundamental matrix

Frect = (H2') \ F21 / H1;

Fcanonical = [0 0  0;
              0 0 -1;
              0 1  0];

Frect_normalized = Frect/norm(Frect, 'fro');
Fcanonical_normalized = ...
    Fcanonical/norm(Fcanonical, 'fro');

verification_error = min( ...
    norm(Frect_normalized-Fcanonical_normalized, 'fro'), ...
    norm(Frect_normalized+Fcanonical_normalized, 'fro'));

disp('Normalized rectified fundamental matrix:');
disp(Frect_normalized);

fprintf('Rectification verification error: %e\n', ...
    verification_error);