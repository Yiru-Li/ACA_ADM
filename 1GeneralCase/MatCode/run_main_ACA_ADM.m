clc; clear all;
close all;
if ~isempty(dir('ACA_ADM.log'))
    diary off;
    copyfile ACA_ADM.log ACA_ADM_copy.log    
    delete('ACA_ADM.log');
end
diary ACA_ADM.log

%%example setup parameters
mu0=1.25663706*10^-6;
eps0=8.85418782*10^-12;
dadt=6.664324407237550e+07;%coil current 
N=[17,17,2]; % the number of dipoles along x,y and z 17 by 17 by 2 is recommended (see publication
%for futher guidance)
FEMord=1;% order of the numerical approximation suggested 1 or 2 
th_hair=.005;%hair thickness or distance of coil from scalp
scth=.05;%width of the square desired scalp coil position search space
load exampletest.mat p te2p conductivity;
%load input data structures and make turn conducitivity into a tensor conductivity
%generate coil
% subplot(1,2,1),
[rs ks]=genfig8(.056/2,.087/2,.006,9);%magstim specs
[t2pcoil,pcoil]=figure8coilmodel(9,.056,.087,.006,.0018,.001);
% axis equal
rs=rs';ks=ks';
%turn conductivity into conductivity tensor
[~,conductivity]=ndgrid(1:9,conductivity(:));
conductivity([2,3,4,6,7,8],:)=0;
%generate ROI
roilist=[1 2.5 5 10 15 20 25 30 35 40 45 50]*1e-3; % m
roictr=[0 0 -.015];
xyz={'x','y','z'};
for dir00=1:1
for iroi=1:1
    roirad=roilist(iroi);
    xx(conductivity(1,:)==.1260)=1; % White
    xx(conductivity(1,:)==.2760)=1; % Grey
    [x y z]=ndgrid(0:pi/150:2*pi,0:pi/150:pi,0:roirad/150:roirad);
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
    clear teid
    teid=1:numel(te2p)/4;
    teid=teid(conductivity(1,teid)==.276);
    
    that=zeros([3 numel(teid)]);
    that(2,:)=1;

    tri=surftri(p',te2p');
    [ROItri,node5]=surftri(p',te2p(:,teid)');
    tri_matters=surftri(p',te2p(:,xx(:)==1)');
    tri_matters2=surftri(p',te2p(:,conductivity(1,:)==0.126)');

    figure;
    trisurf(tri,p(1,:)',p(2,:)',p(3,:)','edgealpha',0,'facealpha',.1,'facecolor',[150 114 100]/255)
    hold on
    trisurf(tri_matters,p(1,:)',p(2,:)',p(3,:)','edgealpha',0,'facealpha',1,'facecolor',[242,174,177]/255)
    trisurf(ROItri,p(1,:)',p(2,:)',p(3,:)','edgealpha',0,'facealpha',1,'facecolor','red'); hold on;
    axis equal
    light
    lighting gouraud

    [pp,Anor,tri]=findSurf(te2p,p,teid,scth,th_hair);
    pcoil(:,4)=1;
    idx = 12; view([0,100]); lightangle(gca,40,40);
    pcoilhead=(Anor(1:3,:,idx)*(pcoil'))';
    % trisurf(t2pcoil,pcoilhead(:,1),pcoilhead(:,2),pcoilhead(:,3),'facecolor',[184 115 51]/256,'edgealpha',0,'facealpha',1);

    xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
    set(gca,'fontsize',12);
    axis off;
    headname=['HeadModel_WholeBrain'];
    saveas(gcf,[headname,'.png']);
    saveas(gcf,[headname,'.fig']);
    hold on; plot3(pp(:,1),pp(:,2),pp(:,3),'k.','markersize',12);
%     legend('Scalp','Grey Matter','ROI','Coil Position');

    headname=['HeadModel_Coil_WholeBrain'];
    saveas(gcf,[headname,'.png']);
    saveas(gcf,[headname,'.fig']);
    teidN=teid;
    AnorN=Anor;
    ppN=pp;
%     disp(['D_{ROI} = ',num2str(roirad*200), ' cm  #teidN = ',num2str(length(teidN)), ...
%          ', #coils = ',num2str(size(AnorN,3))]);    
    disp(['Whole Head ...']);
    pause(2);
    
    omega=dadt;
    disp(['ACA Only U V ...']);
    tic
    [Ux,Vx,MisfitUV,Ik,Jk]=ACA_ADM_ErrUVR(te2p,p,tri,pp,Anor,conductivity,teid,that,rs,ks,omega,scth,th_hair,N,FEMord);
    toc
    fname=['ACA_UVMisfit_WholeBrain','.mat'];
    save(fname, 'Ux','Vx','MisfitUV','Ik','Jk');
    U{iroi}=Ux; V{iroi}=Vx;  MisfitUV0{iroi}=MisfitUV;
    figure; 
    semilogy(MisfitUV,'k.-','linewidth',2,'markersize',20); grid on;
    xlabel('Rank'); ylabel('L^2 Norm Error');
    title('Convergence Curve');
    set(gca,'fontsize',20);
    set(gcf,'position',[15,300,700,500]);    
    compname=['ConvergenceCurve_WholeBrain'];
    saveas(gcf,[compname,'.png']);
    saveas(gcf,[compname,'.fig']);
    
    
%     disp(['Constructing the matrix ...']);
%     tic
%     ACA_MatrixUV=Ux*Vx;
%     toc

    % Calculating Amplitude
%     for ix=1:size(ACA_MatrixUV,2)/3
% %         absROI2coilM(:,ix) = sqrt(ROI2coil_Mat(:,1+3*(ix-1)).^2+ROI2coil_Mat(:,2+3*(ix-1)).^2+ROI2coil_Mat(:,3+3*(ix-1)).^2);
%         absACAMatrixUV(:,ix) = sqrt(ACA_MatrixUV(:,1+3*(ix-1)).^2+ACA_MatrixUV(:,2+3*(ix-1)).^2+ACA_MatrixUV(:,3+3*(ix-1)).^2);
%     end
%     
%     ntet=size(teid,1);
%     absADMm=[];
%     absACA3=reshape(absACAMatrixUV,360,[],ntet);
%     absACAm=squeeze(max(absACA3));
%     
%     fname=['absMax','_d',num2str(roilist(iroi)*2e3),'mm.mat'];
%     save(fname, 'absACAm', 'absADMm', 'ROItri', 'p', 'node5');
% % 
% %     load(fname);
%     maxv0=find(absACAm==max(max(absACAm)));
%     ncoil=size(Anor,3);
%     idx=mod(maxv0(1),ncoil)
% % 
%     figure; 
%     trisurf(ROItri,p(1,:)',p(2,:)',p(3,:)',absACAm(idx,node5),'edgealpha',0,'facealpha',.8);
%     title('ACA Result'); 
%     grid on; axis equal;
%     xlabel('x (m)');ylabel('y (m)'); zlabel('z (m)');
%     set(gca,'fontsize', 15);
%     colorbar;
%     set(gcf,'position',[10,300,610,600])
%     
%     compname=['Comparison_WholeBrain'];
%     saveas(gcf,[compname,'.png']);
%     saveas(gcf,[compname,'.fig']);
%     % save CmpROI_ACA.mat ROI2coil_Mat ACA_Matrix
%     clear absACAMatrixUV Ux Vx ACA_MatrixUV MisfitUV absACA3 absACAm maxv0
end
save(['UV_Mat','.mat'],'U','V','MisfitUV0');
end

copyfile ACA_ADM.log ACA_ADM0.log
diary off
