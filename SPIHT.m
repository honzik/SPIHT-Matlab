function SPIHT(soubor, level, bpp, mode)
% DWT testing function
% parameters:     soubor - input image file                
%                 level - transform depth
%                 bpp - bits per pixel quantifier
%                 mode - 'b' is for base SPIHT, 'd' is for degraded

% residuals coefficient
mult = 10;

% TICK THIS TO PROCESS ONLY PLANE Y
% "grayscales" color image
color = 0;

% init
image = floor(double(rgb2ycbcr(imread(soubor))));

% run for all channels
size_x = size(image,2);
size_y = size(image,1);
planes = size(image,3);

if nargin < 4, mode='b'; end

disp(['image size ' num2str(size_x) ' x ' num2str(size_y)]);

tile=[];
image=image-128;

if color == 0
    image(:,:,2) = zeros(size_y, size_x);
    image(:,:,3) = zeros(size_y, size_x);
    planes=1;
end


% image -> DWT domain
for i=1:planes
    tile(:,:,i) = waveletcdf97(image(:,:,i), level);
end

% dSPIHT encoder
% bpp -> bits count
bpp_full = 8;
bytes = ceil((bpp/bpp_full)*(size_x*size_y*planes));

eff = bytes;
bits = 8*bytes;


% SPIHT coder
disp('');
disp('doing Y plane');
[max_pass, passess, bitstream{1}, timeel] = encodeSPIHT(tile(:,:,1), bits, level, mode);
disp(['ENCODER: Encoding stopped at ' num2str(max_pass - passess) ' after ' num2str(toc) 's, ' num2str(ceil(size(bitstream,2)/8)+1) 'B transmitted']);
timeel = timeel + toc;

disp(['ENCODER: Total time elapsed: ' num2str(timeel) 's.']);  
disp('');
% disp(['doing Cb, Cr planes']);
% for i=2:planes
%     disp(['doing color plane ' num2str(i)]);        
%     [max_pass, passess, bitstream{i}, timeel] = encodeSPIHT(tile(:,:,i), colorBits, level);
% end

% SPIHT decoder
disp('');
disp('doing Y plane');
[max_pass, passess, tile2y, timeel2] = decodeSPIHT(bitstream{1}, mode);
disp(['DECODER: Decoding stopped at ' num2str(max_pass - passess) ' after ' num2str(toc) 's']);
timeel2 = timeel2 + toc;
disp(['DECODER: Total time elapsed: ' num2str(timeel2) 's.']);
disp('');
% disp(['doing Cb, Cr planes']);
% 
% disp(['doing color plane ' num2str(2)]);        
% [max_pass, passess, tile2cb, timeel2] = decodeSPIHT(bitstream{2});
% disp(['doing color plane ' num2str(3)]);        
% [max_pass, passess, tile2cr, timeel2] = decodeSPIHT(bitstream{3});

%tile2(:,:,1) = double(tile2y);
%tile2(:,:,2) = double(tile2cb);
%tile2(:,:,3) = double(tile2cr);

% DWT domain -> image
  
  recon(:,:,1) = waveletcdf97(tile2y, -level);
  recon(:,:,2) = image(:,:,2);
  recon(:,:,3) = image(:,:,3);
%   recon(:,:,2) = waveletcdf97(tile2cb, -level);
%   recon(:,:,3) = waveletcdf97(tile2cr, -level);

recon = recon + 128;
image = ycbcr2rgb(uint8(image+128));
recon = ycbcr2rgb(uint8(recon));


%
% residuals computation
%residuals = abs(recon-image);
%mult = 10;
%residuals = residuals .* mult;
% PSNR + rounding
sum_1 = 0;
for i=1:size_y
    for j=1:size_x
        sum_1 = sum_1 + double((image(i,j,1) - recon(i,j,1))^2);
    end
end

MSE = sum_1/(size_y*size_x);
RMSE = sqrt(double(MSE));
PSNR = 20*log10(255/RMSE)*100;
PSNR = round(PSNR);
PSNR = PSNR/100

% show results
str_orig = ['Original, dimensions=' num2str(size_x) 'x' num2str(size_x) 'px, size=' num2str(size_x*size_y*planes) 'B'];
str_dec = ['DWT-SPIHT result, size=' num2str(bytes) 'B (1:' num2str(round((size_x*size_y*planes)/bytes)) ' / ' num2str(bpp) 'bpp), PSNR=' num2str(PSNR) 'dB'];

figure(1);
subplot(2,1,1);
imshow(uint8(image));
title(str_orig);
    
subplot(2,1,2);
imshow(recon);
title(str_dec);



