clc;
clear all;
close all;

%Things that can be changed later/for diff purposes
resizingscale = 0.5; %How much image is resized (1 for full res)
desiredconnectivity = 26; %This is pixel connectivity. In 3d, 26 means that
%pixels are connected when their faces, edges, or corners touch.
removalsize = 100; %this pixel size and under are removed
imageslide = 1; %image slide to look at in each figure

%If you want to see figures, 1, if not, 0
seeFigures = 1;

%% Image pulled from files, outputted to a final one
imageplace = '/Users/PhoebeDainer/Documents/MATLAB/slicephotos/';
imgname = 'SlicesY0334.tif';
outputfolder = '/Users/PhoebeDainer/Documents/MATLAB/barf/';
outputname = 'Fix_This.tif';


%% Loading Image...
retrievedimg =[imageplace,imgname]; %full image + place

firstimg = imread(retrievedimg); %Reads first image in slideshow (i think, since its 
% a multimage. also type tif!
info = imfinfo(retrievedimg); %finds info of image

image1 = zeros(size(firstimg,1),size(firstimg,2),length(info)); % Construct empty image

for k = 1:length(info) %slice by slice into image
    currentImage = imread(retrievedimg,k,'Info',info);
    image1(:,:,k) = currentImage;
end


%% Image changed size here.
resizedimg = imresize3(image1, resizingscale); %resize 3d img


%% Display figure 1
if seeFigures == 1
    figure; imshow(resizedimg(:,:,imageslide)) %shows a single slice
    title('First image')
end


%% Threshold is done (this is to partition image to fore & back)
thresholdImg = multithresh(resizedimg,1);

%Displays threshold
 disp(['The threshold is ',num2str(thresholdImg)]); 

%Creates logical matrix (changes image to binary)
binary_convert = (resizedimg>thresholdImg); % 1 fore, 0 back


%% Remove small thing from image
cleanedimg = bwareaopen(binary_convert, removalsize, desiredconnectivity);

% Fill holes from removal
filledimg = imfill(cleanedimg, desiredconnectivity,'holes');


%% Display figure 2
if seeFigures == 1
    figure; imshow(filledimg(:,:,imageslide)); 
    title('Image edited & Fixed'); 
end


%% Other part remove & dilate by 2, this is from Sohan's code. 
% Close with spherical structuring element. 
tempImage = ones(size(filledimg)) - filledimg;

% Close this image to clean it
sizeSt = 5;
closedimg = -(imclose(tempImage,strel('sphere', sizeSt)));
closedimg = ones(size(closedimg)) + closedimg;
 
% dilate by about 2 pixels (replaced pixels lost..)
dilationsize = 2;
Idilate = imdilate(closedimg, strel('sphere', dilationsize));

%% Display figure 3
if seeFigures == 1
    figure; imshow(closedimg(:,:,imageslide)); 
    title('Closed, dilated image'); 
end


%% Distance map, made for segmentation (avoid oversaturation) 
distmap = bwdist(~closedimg); 

distmap = -distmap; % Invert distance map (binary associations change)
distmap(~closedimg) = Inf; %INF means infinity. 

%Suppresses all maxima in the intensity image whose height is less than
%maxheight.
maxheight = 5;
distmapsuppressed = -imhmax(-distmap, maxheight, desiredconnectivity);
distmapsuppressed(~closedimg) = Inf;


%% Display Figure 4
if seeFigures == 1
    figure; imshow(distmapsuppressed(:,:,imageslide), [-10, 0]); 
    title('Distance map'); 
end



%% Watershed (retrieved from Matlabs help center)
wshed = watershed(distmapsuppressed,desiredconnectivity);
if seeFigures == 1
    figure; imshow(wshed(:,:,imageslide), [0, 1]); 
    title('Watershed lines'); 
end

%% apply watershed to "original" image (imposes one on the other)
wshed16 = uint16(wshed); % convert to 16-bit
filled16 = uint16(filledimg);

appliedimg = immultiply(wshed16, filled16);


%% ------------FINDING & LOGGING DATA OF PARTICLES!! ---------------

rps = regionprops(appliedimg, 'Area', 'PixelList'); %Area gives actual 
%number of pixels in the region, so it technically is volume here.
%PixelList gives x,y,z coords of pixels. 

side1 = size(appliedimg, 1);
side2 = size(appliedimg, 2);
side3 = size(appliedimg, 3);

minvol = 200;
vol = [rps.Area];
indexStay = find(vol > minvol); %indices of particles greater than a 
%certain volume. 
stats = rps(indexStay); %regionproperties of new image (lil particles gone)

%array for particles over min vol
chosenparticles = zeros(side1, side2, side3);

%array for all the data (x, y, z, mean rad, etc.)
ISlength = numel(indexStay);
mult = 2;
particledata = zeros(ISlength, mult);

%FILLING ARRAYS OF ZEROES WITH DATA!
%ex = x, why = y, zee = z HEHEHEHEHEHEHE EHEHEHEH
for p=1:numel(indexStay)
    %The pixels are voxels for a 3d image! These return their locations
    %with 3 coordinates (indices)
    ex = stats(p).PixelList(:,1);
    why = stats(p).PixelList(:,2);
    zee = stats(p).PixelList(:,3);

    
    particledata(p,1) = mean(ex)/resizingscale; %x
    particledata(p,2) = mean(why)/resizingscale; %y
    particledata(p,3) = mean(zee)/resizingscale; %z
    particledata(p,4) = ((numel(ex))/resizingscale^3*3/4/pi)^(1/3); %mean radii
    %a fifth column CAN be added for mean volume, but kinda unnecessary
    
    % this puts this all on the array of zeroes created abovww
    for q=1:numel(ex)
        chosenparticles(why(q), ex(q), zee(q))=p; %X Y Z ARE OUT OF ORDER
        %TO MAKE FINAL FIGURE LOOK NORMAL! DO NOT CHANGE
    end
end





% %% SAVING THE NEW DATA FILE TO OUTPUT FOLDER
% 
% %retrieved from matlab help
% %open file
% cd(outputfolder)
% save('Particledata_1.mat', 'particledata');
% fID = fopen('Particledata_1.txt','w+');
% 
% %write data
% for n=1:size(particledata,1)
%     fprintf(fID,'%f %f %f %f %i\n',particledata(n,3)/5,particledata(n,2)/5,...
%         particledata(n,1)/5,particledata(n,4)/5);
% end
% 
% %close file
% fclose(fID);
% 
% 
% %% put in folder
% chosenparticles = uint16(chosenparticles);
% imwrite(chosenparticles(:,:,1),outputname);
% 
% for r=2:size(chosenparticles,3)
%     imwrite(chosenparticles(:,:,r),outputname,'WriteMode','append')
% end
% 
% disp([num2str(numel(indexStay)),' final particles']);
% 
% % Final image, to show in figure 6
% rgbfinaldata = label2rgb(chosenparticles(:,:,imageslide),'jet','k','shuffle');
% 
% %% Display Figure 6
% if seeFigures == 1 
%     figure; imshow(rgbfinaldata); 
%     title('Final (segmented) image'); 
% end

