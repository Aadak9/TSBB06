%% Task 5: Triangulation
clear; clc; close all;
labfolder = 'C:\Users\Oneit\Desktop\TSBB06\lab b';
addpath(labfolder);
load(fullfile(labfolder,'tridata.mat'),'x','y1','y2','im1','im2','C1','C2');

figure(1); image([im1 im2]); axis image;
colormap(gray(256)); hold on;
for k = 1:size(y1,2)
    plot(y1(2,k),y1(1,k),'ro');
    plot(y2(2,k)+size(im1,2),y2(1,k),'ro');
end
hold off;

% Use the existing preparatory-exercise triangulation function.
xrec = zeros(4,size(y1,2));
for k = 1:size(y1,2)
    xrec(:,k) = prepQ12(y1(:,k),y2(:,k),C1,C2);
end
assert(all(isfinite(xrec(:))),'Triangulation returned a nonfinite point.');

figure(2); plot3(x(1,:),x(2,:),x(3,:),'o');
axis equal; grid on; title('True 3-D points');
figure(3); plot3(xrec(1,:),xrec(2,:),xrec(3,:),'o');
axis equal; grid on; title('Reconstructed 3-D points');

% Euclidean reconstruction error for each point, in millimetres.
error3D_mm = sqrt(sum((xrec(1:3,:)-x(1:3,:)).^2,1));
mean_error3D_mm = mean(error3D_mm)
std_error3D_mm = std(error3D_mm)

% Reproject the reconstructed points and convert to Cartesian image coordinates.
projected1 = vgg_get_nonhomg(C1*xrec);
projected2 = vgg_get_nonhomg(C2*xrec);
observed1 = vgg_get_nonhomg(y1);
observed2 = vgg_get_nonhomg(y2);
error1_pixels = sqrt(sum((projected1-observed1).^2,1));
error2_pixels = sqrt(sum((projected2-observed2).^2,1));
mean_reprojection1_pixels = mean(error1_pixels)
std_reprojection1_pixels = std(error1_pixels)
mean_reprojection2_pixels = mean(error2_pixels)
std_reprojection2_pixels = std(error2_pixels)

%% Answers: Task 5
% Q: Do the reconstructed points appear as they should?
% ANSWER: Compare Figures 2 and 3; they should show the same 3-D structure.
% Reconstruction errors below quantify the agreement, but do not replace
% checking the plotted shape yourself.
disp(' ');
disp('ANSWER - Shape: compare the true points in Figure 2 with reconstructed points in Figure 3.');

% Q: What are the mean and standard deviation of the 3-D error?
% ANSWER: Compute the Euclidean distance from each reconstructed point to
% its true location, then take the mean and sample standard deviation.
fprintf('\nANSWER - 3-D error: mean %.4f mm; standard deviation %.4f mm.\n', ...
    mean_error3D_mm,std_error3D_mm);

% Q: What are the mean and standard deviation of the reprojection error?
% ANSWER: Project xrec through each camera, divide by homogeneous scale,
% and measure the Euclidean image distance to the supplied observations.
fprintf('\nANSWER - Image 1 reprojection error: mean %.4f pixels; standard deviation %.4f pixels.\n', ...
    mean_reprojection1_pixels,std_reprojection1_pixels);
fprintf('ANSWER - Image 2 reprojection error: mean %.4f pixels; standard deviation %.4f pixels.\n', ...
    mean_reprojection2_pixels,std_reprojection2_pixels);
