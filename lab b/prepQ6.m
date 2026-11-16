
%%Prep Q6
% Load images and camera matrices
im1 = imread('dinosaur0.png');
im2 = imread('dinosaur5.png');

load('dino_Ps','P');

C1 = P{1};
C2 = P{6};

% Calculate camera centres
c1 = null(C1);
c2 = null(C2);

% Calculate epipoles
e12 = C1*c2;
e21 = C2*c1;

% Calculate fundamental matrices
F21 = liu_crossop(e21)*C2*pinv(C1);
F12 = F21';

% Normalize
e12 = e12/e12(3);
e21 = e21/e21(3);
F21 = F21/norm(F21,'fro');
F12 = F12/norm(F12,'fro');

% Calculate Loop-Zhang distortion
w = 720;
h = 576;
lambda = -5:0.01:5;

dist = LoopZhangDistortion(e12,e21,w,h,lambda);

% Plot distortion
figure;
plot(lambda,dist);
set(gca,'YScale','log');
grid on;
xlabel('\lambda');
ylabel('Distortion');
title('Loop-Zhang distortion');

% Display lambda giving the minimum distortion
[minimum_distortion,index] = min(dist);
lambda_min = lambda(index);
minimum_distortion;

