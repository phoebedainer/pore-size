clc;
clear all;
close all;

load('GrainStats_1.mat');

%% varchange
% do you want to see? 0 if no, 1 if yes. 
%fig 1
displaytetra = 0;
displayspheres = 0;
displaycentroids = 1;

%fig 2
displayhistogram = 1;

%changes scale of histogram (when converting from pixels -> cubic mm or whatev)
pixelconvert = 12.5^3;   %12.5^3 for microns^3, 0.0125^3 for mm^3

%% displaying, grabbing data
%lists of centroid x, y, z, and rad
x = grainstats(:,1);
y = grainstats(:,2);
z = grainstats(:,3);
meanrads = grainstats(:,4);

%display centroids
if(displaycentroids == 1)
    scatter3(x,y,z,20,'filled', 'MarkerFaceColor', 'r');
    hold on
end

%triangulation, shows tetramesh
DT = delaunayTriangulation(x,y,z);

if(displaytetra == 1)
    tetramesh(DT,'FaceAlpha',0);
    axis equal
    hold on
end

%displaying spheres.
fin = numel(x);

if(displayspheres == 1)
    for i=1:fin
        [a, b, c] = sphere;

        a = a * meanrads(i);
        b = b * meanrads(i);
        c = c * meanrads(i);

        surf(x(i) + a, y(i) + b, z(i) + c)
        axis equal
        hold on
    end
end

%% Finding volume of pores

poresizes = zeros(size(DT, 1), 1);

%Find volume of each tetrahedron.
for i=1:11

    pt1 = [x(DT(i,1)) y(DT(i,1)) z(DT(i,1))];
    pt2 = [x(DT(i,2)) y(DT(i,2)) z(DT(i,2))];
    pt3 = [x(DT(i,3)) y(DT(i,3)) z(DT(i,3))];
    pt4 = [x(DT(i,4)) y(DT(i,4)) z(DT(i,4))];

    rad1 = meanrads(DT(i, 1));
    rad2 = meanrads(DT(i, 2));
    rad3 = meanrads(DT(i, 3));
    rad4 = meanrads(DT(i, 4));
 
    removethis = spheres_vol(rad1, rad2, rad3, rad4, pt1, pt2, pt3, pt4);
    
    totaltetravol = tetra_total(pt1, pt2, pt3, pt4);

    volumeIwant = totaltetravol - removethis;

    if(volumeIwant > 0)
        poresizes(i)=volumeIwant;
    end
end

%removes all zero values left over in poresizes & reshapes matrix
poresizes(poresizes == 0) = [];

disp(poresizes);

%% Create a histogram of pore sizes, converts pixels to cubic microns/mm/etc.
if(displayhistogram == 1)
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
end

%% Functions

%vol of full tetrahedron
function tetravol=tetra_total(pt1, pt2, pt3, pt4)

    xyzmat=[pt1(1) pt1(2) pt1(3) 1; pt2(1) pt2(2) pt2(3) 1; pt3(1) pt3(2) pt3(3) 1; pt4(1) pt4(2) pt4(3) 1]; %matrix of x, y, z vals
    tetravol = abs(((1/6)*det(xyzmat)));

end

function removedvol=spheres_vol(rad1,rad2,rad3,rad4,pt1, pt2, pt3, pt4)

    %vectors from point 1 out
    pt1vec1 = pt2 - pt1;
    pt1vec2 = pt3 - pt1;
    pt1vec3 = pt4 - pt1;

    %vectors from point 2 out
    pt2vec1 = pt1 - pt2;
    pt2vec2 = pt3 - pt2;
    pt2vec3 = pt4 - pt2;
   
    %vectors from point 3 out
    pt3vec1 = pt1 - pt3;
    pt3vec2 = pt2 - pt3;
    pt3vec3 = pt4 - pt3;
    
    %vectors from point 4 out
    pt4vec1 = pt1 - pt4;
    pt4vec2 = pt2 - pt4;
    pt4vec3 = pt3 - pt4;


    %Removed vol SPHERE 1
    top1 = (abs(pt1vec1 .* pt1vec2 .* pt1vec3));
    bottom1 = (abs(pt1vec1) .* abs(pt1vec2) .* abs(pt1vec3)) + ((dot(pt1vec1,pt1vec2)).*(abs(pt1vec3))) + ((dot(pt1vec1,pt1vec3)).*(abs(pt1vec2))) + ((dot(pt1vec2,pt1vec3)).*(abs(pt1vec1)));
    
    solidangle1 = 2*(atan(top1/bottom1));
    curved_SA_1 = abs((solidangle1 * (rad1^2)));

    removedvol1 = abs(((rad1 * curved_SA_1)/3));

    %Removed vol SPHERE 2
    top2 = (abs(pt2vec1 .* pt2vec2 .* pt2vec3));
    bottom2 = (abs(pt2vec1) .* abs(pt2vec2) .* abs(pt2vec3)) + ((dot(pt2vec1,pt2vec2)).*(abs(pt2vec3))) + ((dot(pt2vec1,pt2vec3)).*(abs(pt2vec2))) + ((dot(pt2vec2,pt2vec3)).*(abs(pt2vec1)));
    
    solidangle2 = 2*(atan(top2/bottom2));
    curved_SA_2 = abs((solidangle2 * (rad2^2)));

    removedvol2 = abs(((rad2 * curved_SA_2)/3));

    %Removed vol SPHERE 3
    top3 = (abs(pt3vec1 .* pt3vec2 .* pt3vec3));
    bottom3 = (abs(pt3vec1) .* abs(pt3vec2) .* abs(pt3vec3)) + ((dot(pt3vec1,pt3vec2)).*(abs(pt3vec3))) + ((dot(pt3vec1,pt3vec3)).*(abs(pt3vec2))) + ((dot(pt3vec2,pt3vec3)).*(abs(pt3vec1)));

    solidangle3 = 2*(atan(top3/bottom3));
    curved_SA_3 = abs((solidangle3 * (rad3^2)));

    removedvol3 = abs(((rad3 * curved_SA_3)/3));

    %Removed vol SPHERE 4
    top4 = (abs(pt4vec1 .* pt4vec2 .* pt4vec3));
    bottom4 = (abs(pt4vec1) .* abs(pt4vec2) .* abs(pt4vec3)) + ((dot(pt4vec1,pt4vec2)).*(abs(pt4vec3))) + ((dot(pt4vec1,pt4vec3)).*(abs(pt4vec2))) + ((dot(pt4vec2,pt4vec3)).*(abs(pt4vec1)));
   
    solidangle4 = 2*(atan(top4/bottom4));
    curved_SA_4 = abs((solidangle4 * (rad4^2)));
    
    removedvol4 = abs(((rad4 * curved_SA_4)/3));


    %All removed vols!
    removedvol = removedvol1 + removedvol2 + removedvol3 + removedvol4;
    
end