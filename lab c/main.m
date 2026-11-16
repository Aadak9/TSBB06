clc;
close all;

%%P olynomial Expansion
k = 0:1:100;
s = sin(k/10);
figure(1);plot(s);

x = (-3: 3)';
b0 = ones(7, 1);
b1 = x;
b2 = x.^2;
figure(2);
subplot(4, 1, 1);plot(b0, 'o');
subplot(4, 1, 2);plot(b1, 'o');
subplot(4, 1, 3);plot(b2, 'o');

a = exp(-x.^2/4);
figure(2);
subplot(4, 1, 4);plot(a, 'o');

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
Answer:
%}


c = inv(G)*[h0;h1;h2];
figure(4);
subplot(3, 1, 1);plot(c(1,:));
subplot(3, 1, 2);plot(c(2,:));
subplot(3, 1, 3);plot(c(3,:));
%{
Q: The last three plots show the elements of the coordinate vector
c. Are they as expected? Hint: polynomial expansion is related to a Taylor
expansion of a function, as described in preparatory exercise 2. In what
way?
Answer:
%}

figure(5);
localsig = s(60-3:60+3);
reconsig = (B*c(:,60))';
diffsig = localsig - reconsig;
subplot(3, 1, 1);plot(localsig);
subplot(3, 1, 2);plot(reconsig);
subplot(3, 1, 3);plot(diffsig);
%{
Q: s the reconstructed signal very different or very similar to the
local signal? Explain why.
Answer: 
%}


%% Uncertain Data
cert = double(rand(1,101)>0.2);
scert = s.*cert;
figure(6);plot(scert);

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
subplot(3, 1, 1);plot(c(1,:));
subplot(3, 1, 2);plot(c(2,:));
subplot(3, 1, 3);plot(c(3,:));


figure(8);
localsig = scert(60-3:60+3);
reconsig = (B*c(:,60))';
diffsig = localsig - reconsig;
subplot(3, 1, 1);plot(localsig);
subplot(3, 1, 2);plot(reconsig);
subplot(3, 1, 3);plot(diffsig);
%{
Q: In what way do the results change? Explain why.
Answer:
%}

h0 = conv(scert, f0, 'same');
h1 = conv(scert, h1, 'same');

G11 = conv(cert, f0, 'same');
G12 = conv(cert, f1, 'same');
G22 = conv(cert, f2, 'same');

detG = G11.*G22-G12.^2;
c0 = (G22.*h0-G12.*h1)./detG;
c1 = (-G12.*h0+G11.*h1)./detG; figure(9);
subplot(2, 1, 1);plot(c0);subplot(2, 1, 2);plot(c1);
%{
Q:How do you interpret the computations that are made for c0
and c1?
Answer:

Q: Compare this result, (c0, c1), to the previous one where the
certainty function was not taken into account. What it is the difference?
Answer: 

Q: The signals c0 and c1 should not exhibit any border effects. In
the previous task you did get border effects when you computed coordinates
c0 and c1 (see figure 4). Explain this observation.
Answer:

Q: Modify the signal by removing a higher rate of the samples
(e.g. 50%) and see how this changes the result. How large rate of samples
can be removed without too much distortion in the measured coordinates
(c0, c1)?
Answer:

%}

%% Normalised averaging of images
im = double(imread('Scalespace0.png'));
figure(10);colormap(gray);imagesc(im);

cert = double(rand(size(im)) > 0.6); imcert = im.*cert;
figure(11);colormap(gray); imagesc(imcert);

x = ones(7, 1)*(-3:3);
y = x';
a = exp(-(x.^2+y)/4);
figure(12);mesh(a);

imlp = conv2(imcert, a, 'same');
figure(13);colormap(gray);imagesc(imlp);
%{
Q: Did this solve the problem of missing pixels? Explain why.
Answer: 
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
    axis image off;
    set(gca, 'Color', 'k');
    title(sprintf('%g%% missing, width %g\n%d x %d, %d NaNs', ...
        missing(k), width(k), size(a,1), size(a,2), nnz(isnan(c))));
end
%{
Q: Is this result better than before? Explain why?
Answer:

Q: Change the width of the applicability function to smaller or
larger values. What happens?
Answer:

Q: Increase the number of missing pixels. How large percentage
can be removed until you cannot see the picture in the reconstructed image?
Answer:


%}


