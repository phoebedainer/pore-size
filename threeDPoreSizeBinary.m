clc;
clear all;
close all;
tic
%% loading everything
load('GrainStats_1.mat');

imageplace = '/Users/PhoebeDainer/Documents/MATLAB/segmented/';
imgname = '13_single_binary.tif';

retrievedimg =[imageplace,imgname]; %full image + place
firstslice = imread(retrievedimg, 'tif'); %Reads first image in slides

info = imfinfo(retrievedimg); %finds info of image

binaryimg = zeros(size(firstslice,1),size(firstslice,2),length(info)); % Construct empty image

for k = 1:size(info,1) %slice by slice into image
    currentImage = imread(retrievedimg,k,'Info',info);
    binaryimg(:,:,k) = currentImage;
end

imshow(binaryimg(:, :, 120));


%% finding volume.

%lists of centroid x, y, z, and that corresponding sphere's rad
x = grainstats(:,1);
y = grainstats(:,2);
z = grainstats(:,3);
meanrads = grainstats(:,4);

%triangulation
DT = delaunayTriangulation(x,y,z);


poresizes = zeros(size(DT,1),1);

for i=1:10

    pt1 = [x(DT(i,1)) y(DT(i,1)) z(DT(i,1))];
    pt2 = [x(DT(i,2)) y(DT(i,2)) z(DT(i,2))];
    pt3 = [x(DT(i,3)) y(DT(i,3)) z(DT(i,3))];
    pt4 = [x(DT(i,4)) y(DT(i,4)) z(DT(i,4))];
    
    %find max/min of slices, to only parse through those
    fourzs = [z(DT(i,1)); z(DT(i,2)); z(DT(i,3)); z(DT(i,4))];
    slicemin = floor(min(fourzs));
    slicemax = ceil(max(fourzs));
    
    %find max/min of columns, to only parse through those
    fourxs = [x(DT(i,1)); x(DT(i,2)); x(DT(i,3)); x(DT(i,4))];
    colmin = floor(min(fourxs));
    colmax = ceil(max(fourxs));
    
    %find max/min of rows, to only parse through those
    fourys = [y(DT(i,1)); y(DT(i,2)); y(DT(i,3)); y(DT(i,4))];
    rowmin = floor(min(fourys));
    rowmax = ceil(max(fourys));

    pixels = 0; 

    [length, width, depth] = size(binaryimg);

    %parse through cube, check if each pixel is 0 & in tetra
    for slice = slicemin:slicemax
        for col = colmin:colmax
            for row = rowmin:rowmax
            
                value = binaryimg(row, col, slice);
                point = [col, row, slice]; %column is x, row is y
            
                if(value == 0) %value == 0 because void is 0
                    yesorno = checkInTetra(pt1, pt2, pt3, pt4, point); 
                    if(yesorno == 1)
                        pixels = pixels + 1;
                    end
                end

            end
        end
        
    end
   
    poresizes(i) = pixels; %adds new pixelcount!
    disp(['New pixelcount added.', num2str(i)]);
end

%removes all zero values left over in poresizes & reshapes matrix
poresizes(poresizes == 0) = [];

disp(poresizes(:,1));

%% Create a histogram of pore sizes, converts pixels to cubic microns/mm/etc.
pixelconvert = 12.5^3;
    figure;
    histogram(poresizes*(pixelconvert));
    title('Pore size distribution');
    
    xlab = 'unknown';
    if(pixelconvert == (12.5^3))
        xlab = 'microns^3';
    elseif(pixelconvert == (0.0125^3))
        xlab = 'millimeters^3';
    end

    xlabel(['Pore volume, in ', xlab])
    ylabel('Frequency')

toc

%% Functions
%check if something is in tetrahedron. 0 if not, 1 if yes.
function t = checkInTetra(pt1, pt2, pt3, pt4, pt)
    m0 = [pt1 1; pt2 1; pt3 1; pt4 1];
    m1 = [pt 1; pt2 1; pt3 1; pt4 1];
    m2 = [pt1 1; pt 1; pt3 1; pt4 1];
    m3 = [pt1 1; pt2 1; pt 1; pt4 1];
    m4 = [pt1 1; pt2 1; pt3 1; pt 1];

    d0 = det(m0);
    d1 = det(m1);
    d2 = det(m2);
    d3 = det(m3);
    d4 = det(m4);

    if(isneg(d0) && isneg(d1) && isneg(d2) && isneg(d3) && isneg(d4))
        t = 1;
    elseif(isneg(d0) == 0 && isneg(d1) == 0 && isneg(d2) == 0 && isneg(d3) == 0 && isneg(d4) == 0)
        t = 1;
    else
        t = 0;
    end

end

%check if something is negative. 1 if yes (negative), 0 if no (positive).
function answer = isneg(value)
    if value < 0
        answer = 1;
    else
        answer = 0;
    end
end
