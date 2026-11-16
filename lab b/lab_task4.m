%% Task 4: Fundamental matrix and eight-point algorithm
clear; clc; close all;
labfolder = 'C:\Users\Oneit\Desktop\TSBB06\lab b';
addpath(labfolder);
im1 = imread(fullfile(labfolder,'dinosaur0.png'));
im2 = imread(fullfile(labfolder,'dinosaur5.png'));
load(fullfile(labfolder,'dino_Ps.mat'),'P');
C1 = P{1}; C2 = P{6};
c1 = null(C1);
e21 = C2*c1;
F21 = liu_crossop(e21)*C2*pinv(C1);

figure(1); set(gcf,'WindowState','maximized');
disp('Select at least 8 matching pairs at different depths. Right-click to finish.');
[y1,y2] = correspondences_select(im1,im2);
assert(size(y1,2)>=8,'The eight-point algorithm needs at least eight pairs.');
y1(end+1,:) = 1;
y2(end+1,:) = 1;

% The supplied function forms A, applies SVD, reshapes F, and enforces rank 2.
F = fmatrix_nn8pa(y1,y2);
% With exactly 8 pairs, A is 8-by-9: its ninth singular value is implicitly zero.

% Compare line distances using the same correspondences for both matrices.
matrices = {F,F21};
labels = {'Estimated F','Camera-derived F21'};
distance1_pixels = zeros(2,size(y1,2));
distance2_pixels = zeros(2,size(y2,2));
for j = 1:2
    l2 = matrices{j}*y1;
    l1 = matrices{j}'*y2;
    for k = 1:size(y1,2)
        l1(:,k) = l1(:,k)/norm(l1(1:2,k));
        l2(:,k) = l2(:,k)/norm(l2(1:2,k));
    end
    figure(j+1); clf;
    subplot(1,2,1); image(im1); axis image; hold on; title(labels{j});
    subplot(1,2,2); image(im2); axis image; hold on; title(labels{j});
    for k = 1:size(y1,2)
        subplot(1,2,1); plot(y1(1,k),y1(2,k),'ro'); drawline(l1(:,k));
        subplot(1,2,2); plot(y2(1,k),y2(2,k),'ro'); drawline(l2(:,k));
    end
    distance1_pixels(j,:) = abs(sum(y1.*l1,1));
    distance2_pixels(j,:) = abs(sum(y2.*l2,1));
    fprintf('\n%s: image 1 distances in pixels\n',labels{j});
    disp(distance1_pixels(j,:));
    fprintf('%s: image 2 distances in pixels\n',labels{j});
    disp(distance2_pixels(j,:));
end
mean_distance_pixels = mean([distance1_pixels distance2_pixels],2)
% Rows correspond to estimated F and camera-derived F21 respectively.
% A smaller error on the fitted points alone does not prove better geometry.

%% Answers: Task 4
% Q: Is there only one significantly small singular value?
% ANSWER: Inspect the values at sufficient precision. Exactly eight pairs
% produce an 8-by-9 A with an implicit ninth zero singular value. Another
% very small value can indicate an unstable or nearly degenerate solution.
A = zeros(size(y1,2),9);
for k = 1:size(y1,2)
    A(k,:) = reshape(y2(:,k)*y1(:,k)',1,9);
end
singular_values = svd(A);
fprintf('\nANSWER - Singular values (without display rounding):\n');
fprintf('%.9g\n',singular_values);
if size(A,1)==8
    disp('There is also an implicit ninth zero singular value because A has 8 rows and 9 columns.');
end
disp('More than one near-null direction makes the estimate ambiguous or unstable.');
disp('Inspect the separation of the smallest values; raw unnormalized values alone cannot establish geometric degeneracy.');

% Q: How can you compare the estimated F with F21?
% ANSWER: Use both matrices to generate epipolar lines for the same points.
% Compare normalized point-to-line distances in pixels, not raw entries,
% since a fundamental matrix is only defined up to a nonzero scale.
disp(' ');
disp('ANSWER - Comparison: Figures 2 and 3 show the lines for F and F21; compare their pixel distances.');

% Q: Is F better or worse than F21? Why?
fprintf('\nANSWER - Mean pixel distance: estimated F %.4f; camera-derived F21 %.4f.\n', ...
    mean_distance_pixels(1),mean_distance_pixels(2));
if mean_distance_pixels(1)<mean_distance_pixels(2)
    disp('Estimated F fits the selected points more closely. It was fitted to these same clicks, including their errors.');
elseif mean_distance_pixels(1)>mean_distance_pixels(2)
    disp('Camera-derived F21 fits the selected points more closely; estimation noise or poor conditioning may affect F.');
else
    disp('Both matrices have the same mean distance on the selected points.');
end
disp('This fitted-point comparison alone does not prove which matrix generalizes better to the actual scene.');
