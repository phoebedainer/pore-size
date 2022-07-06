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

%% Getting Area approx

AreaMatrix = zeros(size(DT,1),1);

%GETTING AREA TO REMOVE
for i=1:size(DT,1)

    pt1 = [allxs(DT(i,1)) allys(DT(i,1))];
    pt2 = [allxs(DT(i,2)) allys(DT(i,2))];
    pt3 = [allxs(DT(i,3)) allys(DT(i,3))];
    
   rad1 = allrads(DT(i, 1));
   rad2 = allrads(DT(i, 2));
   rad3 = allrads(DT(i, 3));
 
   
   areaoftriangle = findTriAreaApprox(pt1, pt2, pt3)
   porearea = findTriAreaApprox(pt1, pt2, pt3) - areaToRemoveApprox(rad1,rad2,rad3,pt1,pt2,pt3);
    
   if(porearea < 0)
       porearea = 0
   end

AreaMatrix(i,1) = porearea; 
% AreaMatrix(i, 2) = findTriAreaApprox(pt1,pt2,pt3);
% AreaMatrix(i,3) = areaToRemoveApprox(rad1,rad2,rad3,pt1,pt2,pt3);
end

disp('Pore area matrix: ')
disp(AreaMatrix); %column 1 is good area(one we want), col 2 is total 
% without subtracting circles, and col 3 is area TO subtract.


%% FUNCTIONS (findTriArea, areaToRemove, and circle!)

%function to find area of one triangle (including intersections
%with circles, which need to be removed.) this should be used in a loop
function trianglearea = findTriAreaApprox(pt1,pt2,pt3)

    pointmatrix =[pt1(1) pt1(2) 1; pt2(1) pt2(2) 1;pt3(1) pt3(2) 1];

    trianglearea = abs(0.5 * det(pointmatrix));
end

%function to find the area where a triangle intersects with a sphere.
%
function removedarea = areaToRemoveApprox(rad1,rad2,rad3,pt1,pt2,pt3)

    vecone1=pt2-pt1;
    vecone2=pt3-pt1;

    vectwo1=pt1-pt2;
    vectwo2=pt3-pt2;

    vecthree1 = pt1-pt3; 
    vecthree2 = pt2-pt3;

    %needs rad of sphere around point 1
    cosineth1 = (dot(vecone1,vecone2))/(norm(vecone1)*norm(vecone2));
    th1 = acos(cosineth1);
    th1deg = (th1*(180/pi));

    ar1 = abs(pi*(rad1*rad1)*(th1deg/360));

    %needs rad of sphere around point 2
    cosineth2 = (dot(vectwo1,vectwo2))/(norm(vectwo1)*norm(vectwo2));
    th2 = acos(cosineth2);
    th2deg = (th2*(180/pi));

    ar2 = abs(pi*(rad2*rad2)*(th2deg/360));

    %needs rad of sphere around point 3
    cosineth3 = (dot(vecthree1,vecthree2))/(norm(vecthree1)*norm(vecthree2));
    th3 = acos(cosineth3);
    th3deg = (th3*(180/pi));

    ar3 = abs(pi*(rad3*rad3)*(th3deg/360));

    removedarea = ar1 + ar2 + ar3;
end

%check in tria
function t = checkInTriangle(P1,P2,P3,P)

    P12 = P1-P2; 
    P23 = P2-P3; 
    P31 = P3-P1;
   t = sign(det([P31;P23]))*sign(det([P3-P;P23])) >= 0 & ...
       sign(det([P12;P31]))*sign(det([P1-P;P31])) >= 0 & ...
       sign(det([P23;P12]))*sign(det([P2-P;P12])) >= 0 ;
    
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
