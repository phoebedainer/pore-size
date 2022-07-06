clc;
clear all;
close all;

%% Vars
rescaleimg = 1;
conn = 8; %8 for 2d 26 for 3d
removalsize = 300;

%% Image pulled
imageplace = '/Users/PhoebeDainer/Documents/MATLAB/slicephotos/';
imgname = 'SlicesY0123.tif';

%% Image loaded
retrievedimg = [imageplace, imgname];
img = imread(retrievedimg, 'tif');

info = imfinfo(retrievedimg); %finds info of image

image1 = zeros(size(img,1),size(img,2)); % Construct empty image

%Puts image onto this image
currentImage = imread(retrievedimg,'Info',info);
image1(:,:) = currentImage;

%im resized
resizedimg = imresize(image1, rescaleimg);

% fig 1
% figure;
% imshow(resizedimg(:,:)); %shows a single slice

%% Threshold is done
thresh = multithresh(resizedimg,1);
disp(['OTSU threshold is: ', num2str(thresh)]);

%Binarization
binary_convert = (resizedimg>thresh); % 1 fore, 0 back

%% Remove hooles
cleanedimg = bwareaopen(binary_convert, removalsize, conn);
filledimg = imfill(cleanedimg, conn,'holes');

% fig 2
% figure;
% imshow(filledimg(:,:));

%% Distance map, made for segmentation (avoid oversaturation) 
distmap = bwdist(~filledimg); 

distmap = -distmap; % Invert distance map (binary associations change)
distmap(~filledimg) = Inf; %INF means infinity. 

%Suppresses all maxima in the intensity image whose height is less than
%maxheight.
maxheight = 4;
distmapsuppressed = -imhmax(-distmap, maxheight, conn);
distmapsuppressed(~filledimg) = Inf;

%% fig 3
% figure; 
% imshow(distmapsuppressed(:,:), [-10, 0]); 

%% Watershed (retrieved from Matlabs help center)
wshed = watershed(distmapsuppressed, conn); 

%% Fig 4
% figure; 
% imshow(wshed(:,:), [0, 1]); 

%% apply watershed to "original" image (imposes one on the other)
wshed8 = uint8(wshed); % convert to 8-bit
filled8 = uint8(filledimg);

appliedimg = immultiply(wshed8, filled8);

%% CALCULATE AND DISPLAY CENTROIDS + TRIANGULATION!!!!
rps = regionprops(appliedimg, 'Centroid', 'Area');
allCentroids = vertcat(rps.Centroid); %mat of centroids
particleareas = [rps.Area]; %mat of particle areas

figure
imshow(filledimg(:,:))
hold on

fin = numel(allCentroids(:,1)); %length of list of centroids

allxs = zeros(fin, 1); %all x vals vert matrix
allys = zeros(fin, 1); %all y vals vert matrix
allrads = zeros(fin, 1); %all radii vert matrix

for k = 1 : fin
    
    x = allCentroids(k, 1);
    y = allCentroids(k, 2);
    
    allxs(k, 1) = x;
    allys(k, 1) = y;

    %calculate mean radius
    thisarea = particleareas(1, k); 
    r = sqrt(thisarea / pi); % r=mean rad

    allrads(k, 1) = r;

%     %display x, y, r
%     disp(['#', num2str(k),'  x        y        r']);
%     disp(['     ', num2str(x), ' ', num2str(y), ' ', num2str(r)]);

    %plot centroids
    plot(x, y, 'r.', 'MarkerSize', 1)
    hold on

    %create circles with mean radii around centroids
    circle(x,y, r);
    hold on

end

%plots triangulation!
DT = delaunay(allxs,allys);
triplot(DT,allxs,allys);
hold on

%% Getting Area 

DTrows = size(DT,1);
AreaMatrix = zeros(DTrows,1);
% totalpixels = 0;

for i=1:DTrows
    
    a = [allxs(DT(i,1), 1), allys(DT(i,1), 1)];
    b = [allxs(DT(i,2), 1), allys(DT(i,2), 1)];
    c = [allxs(DT(i,3), 1), allys(DT(i,3), 1)];

    pixels = 0;
   
    [rows,cols] = size(filledimg);
    for col = 1:cols
        for row = 1:rows

            value = filledimg(row, col);
            point = [col, row];
            
            if(value == 0) %value == 0 because void is 0
                yesorno = checkInTriangle(a,b,c,point); 
                if(yesorno == 1)
                    pixels = pixels + 1;
                end
            end
        end
    end
    
    areaoftriangleapprox = findTriAreaApprox(a, b, c);

    AreaMatrix(i) = pixels;
%     totalpixels = totalpixels + pixels;
    
end

disp('Pore area matrix(number of void pixels): ')
disp(AreaMatrix); 
% disp('Area matrix added up: ')
% disp(totalpixels);


%function to check if a point is in a triangle. return true if so. 
%Retrieved from Roger Stafford on MATLAB Answers, 9 Apr 2016.
% https://www.mathworks.com/matlabcentral/answers/277984-check-points-inside-triangle-or-on-edge-with-example
function t = checkInTriangle(P1,P2,P3,P)

    P12 = P1-P2; 
    P23 = P2-P3; 
    P31 = P3-P1;
   t = sign(det([P31;P23]))*sign(det([P3-P;P23])) >= 0 & ...
       sign(det([P12;P31]))*sign(det([P1-P;P31])) >= 0 & ...
       sign(det([P23;P12]))*sign(det([P2-P;P12])) >= 0 ;
    
end

%finding triangle area
function trianglearea = findTriAreaApprox(pt1,pt2,pt3)

    pointmatrix =[pt1(1) pt1(2) 1; pt2(1) pt2(2) 1;pt3(1) pt3(2) 1];

    trianglearea = abs(0.5 * det(pointmatrix));
end

%function to draw a circle with x&y center and rad r, retrieved from
%internet
function h = circle(x,y,r)
hold on
th = 0:pi/10:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = plot(xunit, yunit);
hold off
end

