clc;
close all;

%%P olynomial Expansion
k = 0:1:100;
s = sin(k/10);
figure(1);plot(s); title('Original sine wave'); xlabel('Sample index'); ylabel('Amplitude');

x = (-3: 3)';
b0 = ones(7, 1);
b1 = x;
b2 = x.^2;
figure(2);
subplot(4, 1, 1);plot(b0, 'o'); title('Constant basis b_0'); xlabel('Window index'); ylabel('Basis value');
subplot(4, 1, 2);plot(b1, 'o'); title('Linear basis b_1'); xlabel('Window index'); ylabel('Basis value');
subplot(4, 1, 3);plot(b2, 'o'); title('Quadratic basis b_2'); xlabel('Window index'); ylabel('Basis value');

a = exp(-x.^2/4);
figure(2);
subplot(4, 1, 4);plot(a, 'o'); title('Gaussian applicability'); xlabel('Window index'); ylabel('Weight');

f0 = b0.*a; f0 = f0(end:-1:1);
f1 = b1.*a; f1 = f1(end:-1:1);
f2 = b2.*a; f2 = f2(end:-1:1);

h0 = conv(s, f0, 'same');
h1 = conv(s, f1, 'same');
h2 = conv(s, f2, 'same');

G0 = diag(a);
B = [b0 b1 b2];
G = B'*G0*B;
%{
Q: Is the basis orthonormal?
Answer:  No, the basis is not orthonormal. The weighted dot product
% of b0 and b2 is sum(a .* x.^2), which is positive. Since this dot
% product is not zero, b0 and b2 are not orthogonal.
%}


c = inv(G)*[h0;h1;h2];
figure(4);
subplot(3, 1, 1);plot(c(1,:)); title('Local value'); xlabel('Sample index'); ylabel('c_0');
subplot(3, 1, 2);plot(c(2,:)); title('Local slope'); xlabel('Sample index'); ylabel('c_1');
subplot(3, 1, 3);plot(c(3,:)); title('Quadratic coefficient'); xlabel('Sample index'); ylabel('c_2');
%{
Q: The last three plots show the elements of the coordinate vector
c. Are they as expected? Hint: polynomial expansion is related to a Taylor
expansion of a function, as described in preparatory exercise 2. In what
way?
Answer: % Answer: Yes, away from the ends of the signal. The local quadratic
% fit c0 + c1*h + c2*h^2 has the same form as a second-order Taylor
% expansion. Here h is the offset from the center sample k.
%
% c(1,:) estimates the signal value       s(k)     = sin(k/10).
% c(2,:) estimates the first derivative   s'(k)    = 0.1*cos(k/10).
% c(3,:) estimates half the second derivative s''(k)/2 = -0.005*sin(k/10).
%
% These are approximate because the fit uses a window of nearby samples.
% Near the ends, convolution also uses zeros outside the signal, which
% can distort the coefficients.
%}

figure(5);
localsig = s(60-3:60+3);
reconsig = (B*c(:,60))';
diffsig = localsig - reconsig;
subplot(3, 1, 1);plot(localsig); title('Original local samples'); xlabel('Window index'); ylabel('Amplitude');
subplot(3, 1, 2);plot(reconsig); title('Local quadratic reconstruction'); xlabel('Window index'); ylabel('Amplitude');
subplot(3, 1, 3);plot(diffsig); title('Original minus reconstruction'); xlabel('Window index'); ylabel('Residual');
%{
Q: s the reconstructed signal very different or very similar to the
local signal? Explain why.
Answer: The reconstructed signal is very similar to the local signal.
% The sine wave is smooth, and over this short window a quadratic
% polynomial captures its value, slope, and curvature well. The small
% difference remains because the sine wave is not exactly quadratic.
%}


%% Uncertain Data
cert = double(rand(1,101)>0.2);
scert = s.*cert;
figure(6);plot(scert); title('Sine wave with missing samples'); xlabel('Sample index'); ylabel('Amplitude');

x = (-3: 3)';
b0 = ones(7, 1);
b1 = x;
b2 = x.^2;
%{
figure(7);
subplot(4, 1, 1);plot(b0, 'o');
subplot(4, 1, 2);plot(b1, 'o');
subplot(4, 1, 3);plot(b2, 'o');
%}

a = exp(-x.^2/4);
%figure(7);
%subplot(4, 1, 4);plot(a, 'o');

f0 = b0.*a; f0 = f0(end:-1:1);
f1 = b1.*a; f1 = f1(end:-1:1);
f2 = b2.*a; f2 = f2(end:-1:1);

h0 = conv(scert, f0, 'same');
h1 = conv(scert, f1, 'same');
h2 = conv(scert, f2, 'same');

G0 = diag(a);
B = [b0 b1 b2];
G = B'*G0*B;



c = inv(G)*[h0;h1;h2];
figure(7);
subplot(3, 1, 1);plot(c(1,:)); title('Masked signal: local value'); xlabel('Sample index'); ylabel('c_0');
subplot(3, 1, 2);plot(c(2,:)); title('Masked signal: local slope'); xlabel('Sample index'); ylabel('c_1');
subplot(3, 1, 3);plot(c(3,:)); title('Masked signal: quadratic coefficient'); xlabel('Sample index'); ylabel('c_2');


figure(8);
localsig = scert(60-3:60+3);
reconsig = (B*c(:,60))';
diffsig = localsig - reconsig;
subplot(3, 1, 1);plot(localsig); title('Masked local samples'); xlabel('Window index'); ylabel('Amplitude');
subplot(3, 1, 2);plot(reconsig); title('Reconstruction with fixed G'); xlabel('Window index'); ylabel('Amplitude');
subplot(3, 1, 3);plot(diffsig); title('Masked samples minus reconstruction'); xlabel('Window index'); ylabel('Residual');
%{
Q: In what way do the results change? Explain why.
Answer: The coefficient curves become less smooth, and the reconstructed
% signal can deviate more from the original sine wave. Missing samples are
% set to zero, but G0 still gives them their original weights. The fit
% therefore treats those zeros as real measurements. This can pull the
% fitted value toward zero and distort the estimated slope and curvature.
%}

h0 = conv(scert, f0, 'same');
h1 = conv(scert, f1, 'same');

G11 = conv(cert, f0, 'same');
G12 = conv(cert, f1, 'same');
G22 = conv(cert, f2, 'same');

detG = G11.*G22-G12.^2;
c0 = (G22.*h0-G12.*h1)./detG;
c1 = (-G12.*h0+G11.*h1)./detG; figure(9);
subplot(2, 1, 1);plot(c0); title('Certainty-weighted local value'); xlabel('Sample index'); ylabel('c_0'); subplot(2, 1, 2);plot(c1); title('Certainty-weighted local slope'); xlabel('Sample index'); ylabel('c_1');
%{
Q:How do you interpret the computations that are made for c0
and c1?
Answer: At each sample, we multiply the filter outputs [h0; h1]
% by the inverse of the local matrix G. This gives the coefficients
% of the fitted line c0 + c1*h, where h is the offset from the center.
% c0 estimates the signal value at the center, and c1 estimates how
% much the signal changes per sample. The certainty mask gives missing
% samples zero weight in the fit.

Q: Compare this result, (c0, c1), to the previous one where the
certainty function was not taken into account. What it is the difference?
Answer: 

Q: The signals c0 and c1 should not exhibit any border effects. In
the previous task you did get border effects when you computed coordinates
c0 and c1 (see figure 4). Explain this observation.
Answer: Near the ends, convolution pads the signal with zeros.
% In Figure 4, G was fixed and gave those padded positions the same
% weights as real samples, so the zeros distorted the fit.
% Now we compute G from cert at each position. Outside the signal,
% cert is also padded with zeros, so those positions get zero weight.
% The fit uses only the available samples, removing the distortion
% caused by treating the padded zeros as measurements.

Q: Modify the signal by removing a higher rate of the samples
(e.g. 50%) and see how this changes the result. How large rate of samples
can be removed without too much distortion in the measured coordinates
(c0, c1)?
Answer: For this sine wave and seven-sample window, about 30% random
% removal is a reasonable limit if we want estimates at almost every
% position. In my simulation, fewer than 1% of positions had no estimate
% at 30% removal, compared with about 7% at 50% removal.
%
% At 50%, c0 usually stays close to the sine wave where the fit is defined,
% but c1 is more sensitive to which samples are missing. If fewer than
% two samples are available in a window, we cannot fit a unique line.

%}

%% Normalised averaging of images
im = double(imread('Scalespace0.png'));
figure(10);colormap(gray);imagesc(im); title('Original grayscale image'); xlabel('Image column (pixels)'); ylabel('Image row (pixels)');

cert = double(rand(size(im)) > 0.6); imcert = im.*cert;
figure(11);colormap(gray);imagesc(imcert); title('Image with missing pixels'); xlabel('Image column (pixels)'); ylabel('Image row (pixels)');

x = ones(7, 1)*(-3:3);
y = x';
a = exp(-(x.^2+y.^2)/4);
figure(12);mesh(a); title('Applicability kernel'); xlabel('Kernel column'); ylabel('Kernel row');

imlp = conv2(imcert, a, 'same');
figure(13);colormap(gray);imagesc(imlp); title('Unnormalised weighted image'); xlabel('Image column (pixels)'); ylabel('Image row (pixels)');
%{
Q: Did this solve the problem of missing pixels? Explain why.
Answer: Only partly. Low-pass filtering can give a missing pixel a
% value based on its neighbors, but it treats missing pixels as black
% because they were set to zero. This makes areas with many missing
% pixels too dark. Normalised averaging fixes this by dividing the
% weighted sum by the total weight of the available pixels.
%}

G = conv2(cert, a, 'same');
c = imlp./G;
rng(1);
r = rand(size(im));

missing = [60 60 60, 90 95 99, 90 90 90];
width   = [ 1  2  4,  2  2  2,  2  4  6];
radius  = [ 3  3  3,  3  3  3,  3  6  9];

figure(14); clf; colormap(gray);

for k = 1:9
    cert = double(r > missing(k)/100);
    imcert = im.*cert;

    [x, y] = meshgrid(-radius(k):radius(k));
    a = exp(-(x.^2+y.^2)/width(k)^2);

    imlp = conv2(imcert, a, 'same');
    G = conv2(cert, a, 'same');
    c = imlp./G;

    subplot(3,3,k);
    img = imagesc(c, [min(im(:)) max(im(:))]);
    set(img, 'AlphaData', ~isnan(c));
    axis image; axis on; xlabel('Image column (pixels)'); ylabel('Image row (pixels)');
    set(gca, 'Color', 'k');
    title(sprintf('%g%% missing, width %g\n%d x %d, %d NaNs', ...
        missing(k), width(k), size(a,1), size(a,2), nnz(isnan(c))));
end
%{
Q: Is this result better than before? Explain why?
Answer: Yes. Normalised averaging gives a better brightness estimate
% because it divides the weighted sum by the total weight of the known
% pixels. Missing pixels have zero weight, so they no longer pull the
% average toward black. The result can still look blurry because each
% estimate is an average of nearby pixels.

Q: Change the width of the applicability function to smaller or
larger values. What happens?
Answer: A narrower applicability function gives more weight to the
% closest pixels. This preserves edges and fine details better, but the
% estimate can be less reliable when the closest pixels are missing.
%
% A wider applicability function gives farther pixels more influence.
% This makes the estimate less sensitive to missing pixels nearby, but
% it can blur edges by averaging pixels from both sides.

Q: Increase the number of missing pixels. How large percentage
can be removed until you cannot see the picture in the reconstructed image?
Answer: It is some confirmation bias but it is very easy to distinguish the
picture with 90% of the pixels missing.


%}


