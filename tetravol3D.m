clc
clearvars
close all

rad = 0.1;

rng("default");
x = rand([10, 1]);
y = rand([10, 1]);
z = rand([10, 1]);
DT = delaunayTriangulation(x,y,z);

tetramesh(DT,'FaceAlpha',0.2); %3d delaunay triangulation
hold on

% start sphere
[a, b, c] = sphere;
 a = a * rad;
 b = b * rad;
 c = c * rad;

%output sphere at every vertex
for i = 1:length(x)
    surf(x(i) + a, y(i) + b, z(i) + c) %displacement
    hold on
end

%-----------FINDING VOLUME-----------%

tetravolumes = zeros(size(DT, 1), 1);

%Find volume of each tetrahedron.
for i=1:size(DT, 1)

    pt1=[x(DT(i,1)) y(DT(i,1)) z(DT(i,1))];
    pt2=[x(DT(i,2)) y(DT(i,2)) z(DT(i,2))];
    pt3=[x(DT(i,3)) y(DT(i,3)) z(DT(i,3))];
    pt4=[x(DT(i,4)) y(DT(i,4)) z(DT(i,4))];

    these_spheres_vol = spheres_vol(rad, pt1, pt2, pt3, pt4);
    
    this_tetra_vol = tetra_total(pt1, pt2, pt3, pt4);

    volumeIwant = this_tetra_vol - these_spheres_vol; %volume of a tetrahedron!!!! yum yum
    tetravolumes(i)=volumeIwant;
end

%MATRIX OF ALL THE TETRAHEDRONS' VOLUMES.
disp(tetravolumes);

function tetravol=tetra_total(pt1, pt2, pt3, pt4)
%VOLUME OF TETRAHEDRON TOTAL (breaks if more than 1..) 

    xyzmat=[pt1(1) pt1(2) pt1(3) 1; pt2(1) pt2(2) pt2(3) 1; pt3(1) pt3(2) pt3(3) 1; pt4(1) pt4(2) pt4(3) 1]; %matrix of x, y, z vals
    tetravol = abs(((1/6)*det(xyzmat)));

end



function removedvol=spheres_vol(rad, pt1, pt2, pt3, pt4)

    %vectors.
    a = pt2 - pt1;
    b = pt3 - pt1;
    c = pt4 - pt1;
    
    %----finding solid angle----
    top = (abs(a .* b .* c));
    bottom = (abs(a) .* abs(b) .* abs(c)) + ((a.*b).*(abs(c))) + ((a.*c).*(abs(b))) + ((b.*c).*(abs(a))) + ((b.*c).*(abs(a)));
    solidangle = 2*(atan(top/bottom));
    %-------------------------

    %finding curved surface area-----
    curved_area = abs((solidangle * (rad^2)));
    %-------

    %finding vol w/ area!
    removedvol = abs(((rad * curved_area)/3));
    
end



