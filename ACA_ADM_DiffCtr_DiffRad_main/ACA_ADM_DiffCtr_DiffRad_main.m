clc; clear all;
close all;
addpath('MatCode');

mu0=1.25663706*10^-6;
eps0=8.85418782*10^-12;
dadt=6.664324407237550e+07;%coil current 
N=[17,17,2]; % the number of dipoles along x,y and z 17 by 17 by 2 is recommended (see publication
%for futher guidance)
FEMord=1;% order of the numerical approximation suggested 1 or 2 
th_hair=.005;%hair thickness or distance of coil from scalp
scth=.05;%width of the square desired scalp coil position search space
load exampletest.mat p te2p conductivity;

load CoilPositions % ROI center list
iroic=32; % ROI center index in the list
roictr=p(:,idxCoilP(iroic)); %ROI center
load(['CaseCoils/CoilPosistion_Pos',num2str(iroic)]);  % Coil positions corresponding to the ROI center

%generate coil
[rs ks]=genfig8(.056/2,.087/2,.006,9);%magstim specs
% [t2pcoil,pcoil]=figure8coilmodel(9,.056,.087,.006,.0018,.001);
rs=rs';ks=ks';
%turn conductivity into conductivity tensor
[~,conductivity]=ndgrid(1:9,conductivity(:));
conductivity([2,3,4,6,7,8],:)=0;

idxx=[1:10];
roilist=5*[1:40]*1e-3; % m
nk0=[20 25 30 30 30 40 50 50 60 60];  % maximum rank set to terminate the iteration, for radii<5 cm, 60 is enough to reach 2% error
for iroi=idxx % loop for diffferent ROI radii
    roirad=roilist(iroi);
    disp(['Runnning ROI Center #',num2str(iroic),' ',mat2str(round(roictr'*1e3,2)),'mm Diameter: ', num2str(roirad*2000),' mm ...']);
    xx(conductivity(1,:)==.1260)=1; % White
    xx(conductivity(1,:)==.2760)=1; % Grey
    % generating ROI
    [x y z]=ndgrid(0:pi/50:2*pi,0:pi/50:pi,0:roirad/50:roirad);
    xp=z(:).*cos(x(:)).*sin(y(:))+roictr(1);
    yp=z(:).*sin(x(:)).*sin(y(:))+roictr(2);
    zp=z(:).*cos(y(:))+roictr(3);
    clear x y z
    TR=triangulation(te2p',p');
    teid=pointLocation(TR,xp,yp,zp);
    teid=teid(isnan(teid)==0);
    teid=unique(teid);
    teid=teid(conductivity(1,teid)==.276);
    clear xp yp zp;
    that=zeros([3 numel(teid)]);
    that(2,:)=1;
    
    % scalp, ROI, Matters surfaces
    tri_Scalp=surftri(p',te2p');
    tri_ROI=surftri(p',te2p(:,teid)');
    tri_matters=surftri(p',te2p(:,xx(:)==1)');
    tri_GM=surftri(p',te2p(:,conductivity(1,:)==0.276)');
    
    % plot ROI in head
    figure;
    trisurf(tri_Scalp,p(1,:)',p(2,:)',p(3,:)','edgealpha',0,'facealpha',.1,'facecolor',[150 114 100]/255); hold on;
    trisurf(tri_matters,p(1,:)',p(2,:)',p(3,:)','edgealpha',0,'facealpha',1,'facecolor',[242,174,177]/255);
    trisurf(tri_ROI,p(1,:)',p(2,:)',p(3,:)','edgealpha',0,'facealpha',1,'facecolor','red'); hold on;
    axis equal; light; lighting gouraud; lightangle(gca,0,50); view([0,90]);
    set(gcf,'Position',[680 120 700 850]);
    
    %[pp,Anor,tri]=findSurf_cluster(te2p,p,roictr,scth,th_hair);
    teidN=teid; AnorN=Anor; ppN=pp;
    disp(['Nte: ',num2str(size(teid,1)),' Ncoil: ',num2str(size(Anor,3))]);
    pause(0.1);
    
    %% ACA_ADM
    omega=dadt; nk=nk0(iroi); % max rank used
    disp(['ACA Only U V ...']);
    fname=['ACA_UVMisfit_Pos',num2str(iroic),'_d',num2str(2000*roirad),'mm.mat'];
    tic
    [Ux,Vx,MisfitUV,Ik,Jk]=ACA_ADM_ErrUVR(fname,nk,te2p,p,ppN,AnorN,conductivity,teid,that,rs,ks,omega,scth,th_hair,N,FEMord);
    toc
    save(fname, 'Ux','Vx','MisfitUV','Ik','Jk','teid','pp','Anor','-v7.3');
end
disp(['Completed.']);

function [pp,Anor,tri]=findSurf_cluster(te2p,p,roicen,scth,th_hair)

[tri]=surftri(p',te2p');
if scth~=0;
    pcen=(p(:,tri(:,1))+p(:,tri(:,2))+p(:,tri(:,3)))/3;
%         roicen=sum(p(:,te2p(1,teid)),2)/numel(te2p(1,teid));

    err=abs(pcen(1,:)-roicen(1))+abs(pcen(2,:)-roicen(2))+abs(pcen(3,:)-roicen(3));
    roicen=pcen(:,err==min(err));
    tri=tri(pcen(1,:)-roicen(1)<scth,:);
    pcen=pcen(:,pcen(1,:)-roicen(1)<scth);
    tri=tri(pcen(1,:)-roicen(1)>-scth,:);
    pcen=pcen(:,pcen(1,:)-roicen(1)>-scth);
    tri=tri(pcen(2,:)-roicen(2)<scth,:);
    pcen=pcen(:,pcen(2,:)-roicen(2)<scth);
    tri=tri(pcen(2,:)-roicen(2)>-scth,:);
    pcen=pcen(:,pcen(2,:)-roicen(2)>-scth);
    tri=tri(pcen(3,:)-roicen(3)<scth,:);
    pcen=pcen(:,pcen(3,:)-roicen(3)<scth);
    tri=tri(pcen(3,:)-roicen(3)>-scth,:);
    pcen=pcen(:,pcen(3,:)-roicen(3)>-scth);
end
[pp,~,tri]=unique(tri(:));
tri=reshape(tri,[numel(tri)/3 3]);
pp=p(:,pp)';
np=numel(pp)/3;
%
 v1=pp(tri(:,2),:)-pp(tri(:,1),:);
 v2=pp(tri(:,3),:)-pp(tri(:,1),:);
 nhat=cross(v1,v2,2);
 for i=1:numel(nhat)/3
     nhat(i,:)=nhat(i,:)/norm(nhat(i,:));
 end
 clear v1 v2 
%extract point normals and local coordinate systems
Anor=zeros([4 4 np]);
v=zeros([3 3]);
nhatp=zeros(size(pp));

for i=1:numel(tri)/3
    v(:,1)=(pp(tri(i,2),:)-pp(tri(i,3),:));
    v(:,1)=v(:,1)/norm(v(:,1));
    v(:,2)=(pp(tri(i,3),:)-pp(tri(i,1),:));
    v(:,2)=v(:,2)/norm(v(:,2));
    v(:,3)=(pp(tri(i,1),:)-pp(tri(i,2),:));
    v(:,3)=v(:,3)/norm(v(:,3));
    for j=1:3
        nhatp(tri(i,j),:)=nhatp(tri(i,j),:)+...
            acos(-sum(v(:,mod(j,3)+1).*v(:,mod(j+1,3)+1)))...
            *nhat(i,:);
    end
end

for i=1:np
    nhatp(i,:)=nhatp(i,:)/norm(nhatp(i,:));
end

thatp1=zeros(size(nhatp));
thatp1(:,1)=1-nhatp(:,1).*nhatp(:,1);
thatp1(:,2)=-nhatp(:,1).*nhatp(:,2);
thatp1(:,3)=-nhatp(:,1).*nhatp(:,3);
for i=1:np
    thatp1(i,:)=thatp1(i,:)/norm(thatp1(i,:));
end
thatp2=cross(nhatp,thatp1);
pp=pp+nhatp*th_hair;
Anor(1:3,1,:)=reshape(thatp1',[3 1 np]);
Anor(1:3,2,:)=reshape(thatp2',[3 1 np]);
Anor(1:3,3,:)=reshape(nhatp',[3 1 np]);
Anor(1:3,4,:)=reshape(pp',[3 1 np]);
Anor(4,4,:)=1;
end

