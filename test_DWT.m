function test_DWT(soubor, wavelet, level, bpp, zobr)
% DWT testing function
% parameters:     soubor - input image file
%                 wavelet - wavelet name or 'cdf97'
%                 level - transform depth
%                 bpp - bits per pixel quantifier
%                 zobr - result 
%                            0 - only input + output images
%                            1 - more details

% residuals coefficient
mult = 10;

% init
image = floor(double(rgb2ycbcr(imread(soubor))));
% run for all channels
size_x = size(image,1);
size_y = size(image,2);
planes = size(image,3);

tile=[];
% image -> DWT domain
if strcmp(wavelet, 'cdf97')
    for i=1:planes
        tile(:,:,i) = waveletcdf97(image(:,:,i), level);
    end
else
  mode = 'per';
  dwtmode(mode);
  type = wavelet;
  [Lo_D,Hi_D,Lo_R,Hi_R] = wfilters(type);
  for i=1:planes
    [tile(:,:,i)] = makeDwtTile(image(:,:,i), level, Lo_D, Hi_D);
  end
end

% SPIHT encoder
% bpp -> bits count
bpp_full = 8;
bytes = ceil((bpp/bpp_full)*(size_x^2));

colorBytes = bytes/4;
bytes = bytes-2*colorBytes;

eff = bytes;
bits = bytes*8;
colorBits = colorBytes*8;


% SPIHT coder
disp('');
disp(['doing Y plane']);
[max_pass, passess, bitstream{1}, timeel] = encodeSPIHT(tile(:,:,1), bits, level);
disp(['ENCODER: Encoding stopped at ' num2str(max_pass - passess) ' after ' num2str(toc) 's, ' num2str(ceil(size(bitstream,2)/8)+1) 'B transmitted']);
timeel = timeel + toc;
disp(['ENCODER: Total time elapsed: ' num2str(timeel) 's.']);  
disp('');
disp(['doing Cb, Cr planes']);
for i=2:planes
    disp(['doing color plane ' num2str(i)]);        
    [max_pass, passess, bitstream{i}, timeel] = encodeSPIHT(tile(:,:,i), colorBits, level);
end




% SPIHT decoder
disp('');
disp(['doing Y plane']);
[max_pass, passess, tile2y, timeel2] = decodeSPIHT(bitstream{1});
disp(['DECODER: Decoding stopped at ' num2str(max_pass - passess) ' after ' num2str(toc) 's']);
timeel2 = timeel2 + toc;
disp(['DECODER: Total time elapsed: ' num2str(timeel2) 's.']);
disp('');
disp(['doing Cb, Cr planes']);

disp(['doing color plane ' num2str(2)]);        
[max_pass, passess, tile2cb, timeel2] = decodeSPIHT(bitstream{2});
disp(['doing color plane ' num2str(3)]);        
[max_pass, passess, tile2cr, timeel2] = decodeSPIHT(bitstream{3});

%tile2(:,:,1) = double(tile2y);
%tile2(:,:,2) = double(tile2cb);
%tile2(:,:,3) = double(tile2cr);

% DWT domain -> image
if strcmp(wavelet, 'cdf97')
  
  recon(:,:,1) = waveletcdf97(tile2y, -level);
  recon(:,:,2) = waveletcdf97(tile2cb, -level);
  recon(:,:,3) = waveletcdf97(tile2cr, -level);
  
else
    
  recon(:,:,1) = reconstructDwtTile(tile2y, Lo_R, Hi_R, level);
  recon(:,:,2) = reconstructDwtTile(tile2cb, Lo_R, Hi_R, level);
  recon(:,:,3) = reconstructDwtTile(tile2cr, Lo_R, Hi_R, level);
 
end

recon = ycbcr2rgb(uint8(recon));
imshow(recon);

return;

%
% residuals computation
%residuals = abs(recon-image);
%mult = 10;
%residuals = residuals .* mult;
% PSNR + rounding
%sum_1 = 0;

for i=1:size(image,1)
    for j=1:size(image,1)
        sum_1 = sum_1 + (image(i,j) - recon(i,j))^2;
    end
end

MSE = sum_1/(size_x^2);
RMSE = sqrt(MSE);
PSNR = 20*log10(255/RMSE)*100;
PSNR = round(PSNR);
PSNR = PSNR/100

% show results
str_orig = ['Original, dimensions=' num2str(size_x) 'x' num2str(size_x) 'px, size=' num2str(size_x^2) 'B'];
str_tile = ['DWT decomposition, wavelet=' wavelet ', level=' num2str(level)];
str_dec = ['DWT-SPIHT result, size=' num2str(ceil(size(bitstream,2)/8)+1) 'B (1:' num2str(round((size_x^2)/(size(bitstream,2)/8))) ' / ' num2str(bpp) 'bpp), PSNR=' num2str(PSNR) 'dB'];
str_res = ['Absolute difference (values magnified by ' num2str(mult) 'x)'];

if zobr == 1
    figure(5);
    
    subplot(1,2,1);
    imagesc(tile);
    colormap('gray');
    title(str_tile);
    
    subplot(1,2,2);
    imagesc(residuals, [0 255]);
    colormap('gray');
    title(str_res);
    
    figure(6);
    mesh(tile);
    title('3D graph of the decomposition image, level=3');
    
end
figure(1);
subplot(1,2,1);
imagesc(image, [0 255]);
colormap('gray');
title(str_orig);
    
subplot(1,2,2);
imagesc(recon, [0 255]);
colormap('gray');
title(str_dec);



