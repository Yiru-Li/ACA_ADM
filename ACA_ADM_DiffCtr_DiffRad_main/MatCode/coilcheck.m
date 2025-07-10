function coilcheck(p,tri_Scalp,Anor,ErrCName)
% clc; clear all;
% close all;
% 
% 
% %%example setup parameters
% mu0=1.25663706*10^-6;
% eps0=8.85418782*10^-12;
% dadt=6.664324407237550e+07;%coil current 
% N=[17,17,2]; % the number of dipoles along x,y and z 17 by 17 by 2 is recommended (see publication
% %for futher guidance)
% FEMord=1;% order of the numerical approximation suggested 1 or 2 
% th_hair=.005;%hair thickness or distance of coil from scalp
% scth=.05;%width of the square desired scalp coil position search space
% load exampletest.mat p te2p conductivity;
% %load input data structures and make turn conducitivity into a tensor conductivity
% %generate coil
% % subplot(1,2,1),
% [rs ks]=genfig8(.056/2,.087/2,.006,9);%magstim specs
[t2pcoil,pcoil]=figure8coilmodel(9,.056,.087,.006,.0018,.001);
%Form a simplified cube for the coil
minx = min(min(pcoil(:,1)),min(pcoil(:,2))); maxx = max(max(pcoil(:,1)),max(pcoil(:,2)));
miny = minx; maxy = maxx;
minz = min(pcoil(:,3)); maxz = max(pcoil(:,3));

cube_vertices = [minx,miny,minz;
                 maxx,miny,minz;
                 minx,maxy,minz;
                 maxx,maxy,minz;
                 minx,miny,maxz;
                 maxx,miny,maxz;
                 minx,maxy,maxz;
                 maxx,maxy,maxz];
T1=delaunay(cube_vertices);
[tcoil]=surftri(cube_vertices,T1);
coil.faces=tcoil;

fid = fopen(ErrCName, 'w');
scalp.faces = tri_Scalp;
scalp.vertices = p';
verbose=0;
Nc=size(Anor,3); idc=floor(Nc/10);
for coil_pos=1:Nc
    if coil_pos==410
        coil_pos;
    end
    for coil_ang=1:1
        if mod(coil_pos,idc)==1 | coil_pos==Nc
            disp(['coil_pos=',num2str(coil_pos),'/',num2str(size(Anor,3)),',coil_ang=',num2str(coil_ang)])
        end
        coil.vertices = rotation_matrix(cube_vertices,Anor,coil_pos,coil_ang);

        %check for intersection
        [intMatrix, intSurface]=SurfaceIntersection(scalp, coil);
        if nnz(intMatrix(:))>0
            fwrite(fid,[num2str(coil_pos),',',num2str(coil_ang),',',newline]);
        end
%         if verbose
%             figure('visible','off');
%             trisurf(scalp.faces,scalp.vertices(:,1),scalp.vertices(:,2),scalp.vertices(:,3),'edgealpha',0,'facealpha',.2)
%             hold on
%             trisurf(coil.faces,coil.vertices(:,1),coil.vertices(:,2),coil.vertices(:,3),'edgealpha',0,'facealpha',.2)
%             if nnz(intMatrix(:))>0
%                 trisurf(intSurface.faces,intSurface.vertices(:,1),intSurface.vertices(:,2),intSurface.vertices(:,3))
%             end
%             pcoilhead = rotation_matrix(pcoil',Anor,coil_pos,coil_ang);
%             trisurf(t2pcoil,pcoilhead(:,1),pcoilhead(:,2),pcoilhead(:,3),'facecolor',[184 115 51]/256,'edgealpha',0,'facealpha',1);
%             hold off
%             axis equal; light; lighting gouraud;
%             lightangle(gca,40,40);
%             headname='Scalp_Coil_Intersection';
%             saveas(gcf,[headname,'.png']);
%             saveas(gcf,[headname,'.fig']);
%         end
    end
end
fclose(fid);


function rotated_points = rotation_matrix(vertices,translation,new_position,rotation_angle)
    point_rot = [];
    point_rot(1,:)=vertices(:,1)*cos(rotation_angle/180*pi)+vertices(:,2)*sin(rotation_angle/180*pi);
    point_rot(2,:)=-vertices(:,1)*sin(rotation_angle/180*pi)+vertices(:,2)*cos(rotation_angle/180*pi);
    point_rot(3,:)=vertices(:,3);
    point_rot(4,:)=1;
    rotated_points = (translation(1:3,:,new_position)*point_rot)';

