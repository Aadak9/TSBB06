%% Task 3: Rectifying homographies
clear; clc; close all;
labfolder = 'C:\Users\Oneit\Desktop\TSBB06\lab b';
addpath(labfolder);
im1 = imread(fullfile(labfolder,'dinosaur0.png'));
im2 = imread(fullfile(labfolder,'dinosaur5.png'));
load(fullfile(labfolder,'dino_Ps.mat'),'P');
C1 = P{1};
C2 = P{6};
[h,w,~] = size(im1);

%% 3.1 Introduction
figure(1); image([im1 im2]); axis image;
c1 = null(C1);
e21 = C2*c1;
F21 = liu_crossop(e21)*C2*pinv(C1);

disp('Select 4-5 matching pairs, one click in each image per pair. Right-click to finish.');
figure(1); set(gcf,'WindowState','maximized');
[y1,y2] = correspondences_select(im1,im2);
assert(size(y1,2)>=4,'Select at least four matching pairs.');
y1(end+1,:) = 1;
y2(end+1,:) = 1;
l2 = F21*y1;
l1 = F21'*y2;
for k = 1:size(y1,2)
    l1(:,k) = l1(:,k)/norm(l1(1:2,k));
    l2(:,k) = l2(:,k)/norm(l2(1:2,k));
end
figure(2); clf;
subplot(1,2,1); image(im1); axis image; hold on;
subplot(1,2,2); image(im2); axis image; hold on;
for k = 1:size(y1,2)
    subplot(1,2,1); plot(y1(1,k),y1(2,k),'ro'); drawline(l1(:,k));
    subplot(1,2,2); plot(y2(1,k),y2(2,k),'ro'); drawline(l2(:,k));
end
distance1_pixels = abs(sum(y1.*l1,1))
distance2_pixels = abs(sum(y2.*l2,1))
[U,~,V] = svd(F21);
e12 = V(:,3); e12 = e12/e12(3);
e21 = U(:,3); e21 = e21/e21(3);
epipole_image1 = e12(1:2)
epipole_image2 = e21(1:2)

%% Answers: Section 3.1
% Q: With what accuracy can you determine corresponding coordinates?
% ANSWER: A reasonable estimate for manual clicking is one to a few pixels.
% This is your visual estimate, not something the script measures directly.
fprintf('\nANSWER - Clicking accuracy: estimate one to a few pixels; adjust this estimate to your actual experience.\n');

% Q: What are the Cartesian epipoles, and where are they?
fprintf('\nANSWER - Epipoles: image 1 (%.2f, %.2f); image 2 (%.2f, %.2f).\n', ...
    e12(1),e12(2),e21(1),e21(2));
% For these dinosaur cameras, e12 is above/right and e21 is above/left,
% both far outside the 720-by-576 images.
disp('For these cameras, the first epipole is above/right and the second is above/left, far outside the images.');

% Q: Do lines pass through matching points? Are they parallel? Consistent?
% ANSWER: Check Figure 2. Small point-to-line distances mean approximate
% intersection. The lines converge at finite epipoles, so they are not
% exactly parallel, but appear nearly parallel because the epipoles are far away.
fprintf('\nANSWER - Line agreement: mean distances %.3f and %.3f pixels. Inspect Figure 2 for intersection.\n', ...
    mean(distance1_pixels),mean(distance2_pixels));
disp('The lines should look nearly parallel, consistent with the distant finite epipoles.');

% Q: Are distances of the same order as estimated clicking accuracy?
% ANSWER: Compare the printed distances with your estimate above. Values
% around 1-3 pixels are comparable; much larger values deserve inspection.
fprintf('\nANSWER - Accuracy comparison: largest distances %.3f and %.3f pixels.\n', ...
    max(distance1_pixels),max(distance2_pixels));
disp('Distances of one to a few pixels agree with that clicking estimate; assess any larger outliers separately.');

% Q: Are the computed epipoles consistent with the observed lines?
% ANSWER: Each epipole should lie on every epipolar line in its image.
% Distant finite epipoles explain the nearly parallel visible lines.
fprintf('\nANSWER - Epipole consistency: maximum line-incidence residuals %.3g and %.3g (ideally zero).\n', ...
    max(abs(e12'*l1)),max(abs(e21'*l2)));

%% 3.2 Loop and Zhang
lambda = -5:0.01:5;
dist = LoopZhangDistortion(e12,e21,w,h,lambda);
[~,index] = min(dist);
lambda_coarse = lambda(index);
figure(3); clf;
subplot(2,1,1); plot(lambda,dist); set(gca,'YScale','log'); grid on;
xlabel('\lambda'); ylabel('Distortion');
subplot(2,1,2); plot((lambda(1:end-1)+lambda(2:end))/2,diff(dist)./diff(lambda));
grid on; xlabel('\lambda'); ylabel('Approximate derivative');

lambda = lambda_coarse-0.02:0.0001:lambda_coarse+0.02;
dist = LoopZhangDistortion(e12,e21,w,h,lambda);
[minimum_distortion,index] = min(dist);
lmin = lambda(index)
minimum_distortion
figure(4); clf;
subplot(2,1,1); plot(lambda,dist); grid on;
xlabel('\lambda'); ylabel('Distortion');
subplot(2,1,2); plot((lambda(1:end-1)+lambda(2:end))/2,diff(dist)./diff(lambda));
grid on; xlabel('\lambda'); ylabel('Approximate derivative');
% The supplied distortion function is used unchanged, as requested by the PM.
% Q: How did you verify the minimum, and what lambda minimizes the error?
% ANSWER: Search [-5,5], refine around the lowest sample, and inspect the
% fine plot: the derivative should change from negative to positive.
fprintf('\nANSWER - Minimum: refined sampled lambda %.4f, distortion %.6f.\n',lmin,minimum_distortion);
if index>1 && index<numel(dist)
    derivative_left = (dist(index)-dist(index-1))/(lambda(index)-lambda(index-1));
    derivative_right = (dist(index+1)-dist(index))/(lambda(index+1)-lambda(index));
    fprintf('Approximate derivatives immediately left/right: %.6g, %.6g.\n',derivative_left,derivative_right);
    disp('A negative left derivative and positive right derivative support a local minimum.');
end
disp('This is the minimum of the supplied function; its implementation does not exactly match the paper criterion.');
lambda = lmin;
w1 = liu_crossop(e12)*[lambda;1;0]; w1 = w1/w1(3);
w2 = F21*[lambda;1;0]; w2 = w2/w2(3);
Hp1 = [1 0 0;0 1 0;w1'];
Hp2 = [1 0 0;0 1 0;w2'];
vcp = 0;
Hr1 = [F21(3,2)-w1(2)*F21(3,3), w1(1)*F21(3,3)-F21(3,1), 0;
       F21(3,1)-w1(1)*F21(3,3), F21(3,2)-w1(2)*F21(3,3), F21(3,3)+vcp;
       0 0 1];
Hr2 = [w2(2)*F21(3,3)-F21(2,3), F21(1,3)-w2(1)*F21(3,3), 0;
       w2(1)*F21(3,3)-F21(1,3), w2(2)*F21(3,3)-F21(2,3), vcp;
       0 0 1];
H1 = Hr1*Hp1;
H2 = Hr2*Hp2;
Frect = H2'\F21/H1;
Frect = Frect/norm(Frect,'fro')
% Expected up to sign: [0 0 0;0 0 -1;0 1 0]/sqrt(2).
% Q: How can you verify H1 and H2 are rectifying homographies? Does it work?
% ANSWER: H2^(-T)*F21*H1^(-1) must be proportional to the canonical
% rectified fundamental matrix, which requires corresponding rows to agree.
canonical = [0 0 0;0 0 -1;0 1 0]/sqrt(2);
rectification_error = min(norm(Frect-canonical,'fro'),norm(Frect+canonical,'fro'));
fprintf('\nANSWER - Homography verification: canonical-matrix error %.3g.\n',rectification_error);
if rectification_error<1e-8
    disp('Yes: the transformed fundamental matrix has the required rectified form, within numerical precision.');
else
    disp('The canonical-matrix check failed; inspect the homographies and camera geometry.');
end

%% 3.3 Rectifying the images
% Flip both vertical axes only when the transformed image is upside down.
centre = [(w+1)/2;(h+1)/2];
direction = map_points(H1,[centre centre+[0;1]]);
if direction(2,2)<direction(2,1)
    Hflip = diag([1 -1 1]);
    H1 = Hflip*H1;
    H2 = Hflip*H2;
end
oldcorners = [1 w w 1;1 1 h h];
newcorners1 = map_points(H1,oldcorners);
newcorners2 = map_points(H2,oldcorners);
mincol = min([newcorners1(1,:) newcorners2(1,:)]);
minrow = min([newcorners1(2,:) newcorners2(2,:)]);
T = [1 0 -mincol+1;0 1 -minrow+1;0 0 1];
newcorners1 = map_points(T*H1,oldcorners);
newcorners2 = map_points(T*H2,oldcorners);
maxcol = max([newcorners1(1,:) newcorners2(1,:)]);
maxrow = max([newcorners1(2,:) newcorners2(2,:)]);
T = inv(diag([maxcol/w maxrow/h 1]))*T
imr1 = image_resample(double(im1),T*H1,h,w);
imr2 = image_resample(double(im2),T*H2,h,w);
figure(5); image(uint8([imr1 imr2])); axis image;
points1_rectified = map_points(T*H1,y1(1:2,:))
points2_rectified = map_points(T*H2,y2(1:2,:))
row_difference_pixels = points1_rectified(2,:)-points2_rectified(2,:)

%% Answers: Section 3.3
% Q: How can you see T is a scaling and translation transformation?
% ANSWER: T has the form [sx 0 tx;0 sy ty;0 0 1]. Its diagonal entries
% scale x/y, and its last column translates them. The corner coordinates
% transformed by both H1 and H2 determine these scales and translations.
fprintf('\nANSWER - T: x_new = %.6g*x + %.6g; y_new = %.6g*y + %.6g.\n', ...
    T(1,1),T(1,3),T(2,2),T(2,3));
disp('Here x/y are already rectified coordinates; the common T fits both images into the output area.');

% Q: Do the transformed images and points appear rectified?
% ANSWER: Inspect Figure 5. Matching points should lie on approximately the
% same row; remaining discrepancies reflect the selected correspondences.
fprintf('\nANSWER - Rectified points: mean absolute row difference %.3f pixels; largest %.3f pixels.\n', ...
    mean(abs(row_difference_pixels)),max(abs(row_difference_pixels)));
disp('Small row differences support rectification. Confirm the visual appearance in Figure 5.');
