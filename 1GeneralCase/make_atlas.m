clc; clear all;
close all;
addpath('MatCode');
diary([char(datetime, 'yyyyMMddHHmmss') '.log'])
sublist=[10171,10206,10228,10249,10280,10304,10316,10321,10339];

roilist=[1 2.5 5 10 15 20 25 30 35 40 45 50]*1e-3; % m
roictr0={[-0.02 -0.02 0.08],[-0.02 -0.02 0.04],[-0.02 -0.02 0.04], ...
    [-0.02 -0.02 0.04],[-0.02 -0.02 0.06],[-0.02 -0.02 0.06], ...
    [-0.02 -0.02 0.04],[-0.02 -0.02 0.06],[-0.02 -0.02 0.05]};
scth0=[0.02 .03 .04 .05 .06 .07 .08 .09 .1 .12 .14]; %width of the square desired scalp coil position search space
for k=1%:9
    subn=sublist(k);
    disp(['Running Head sub-',num2str(subn),' ...']);
    roictr=roictr0{k};
    %%example setup parameters
    mu0=1.25663706*10^-6;
    eps0=8.85418782*10^-12;
    dadt=6.664324407237550e+07;%coil current
    N=[17,17,2]; % the number of dipoles along x,y and z 17 by 17 by 2 is recommended (see publication
    %for futher guidance)
    FEMord=1; % order of the numerical approximation suggested 1 or 2
    th_hair=.005; %hair thickness or distance of coil from scalp
    %load input data structures and make turn conducitivity into a tensor conductivity
    load(['HeadModels/example_sub-',num2str(subn),'.mat']);
    % conductivity = tetrahedra conductivity (Nx1)
    % p = nodes (3xA)
    % reg = tetrahedra region (Nx1)
    % t2p = triangle connectivity (3xB)
    %te2p = tetrahedra connectivity (4xN)
    p=p*1e-3;
    TR=triangulation(te2p',p');
    %generate coil
    [rs, ks]=genfig8(.056/2,.087/2,.006,9);%magstim specs
    rs=rs';ks=ks';
    %turn conductivity into conductivity tensor
    [~,conductivity]=ndgrid(1:9,conductivity(:));
    conductivity([2,3,4,6,7,8],:)=0;

    % select whole grey matter as ROI
    teid=[1:numel(te2p)/4]'; % all tetrahedra
    teid=teid(conductivity(1,teid)==.276 | conductivity(1,teid)==.126); % tetrahedra corresponding to GM/WM
    Nt=size(te2p,2);
    xx = [];
    xx(conductivity(1,:)==.1260)=1; % White
    xx(conductivity(1,:)==.2760)=1; % Grey
    teid0=(1:Nt)'; teid0=teid0(xx==1);
    load(['CaseCoils/CoilPosistion_sub-',num2str(subn),'.mat']);
    % pp = coil positions
    % Anor = coil position with normal vector
    disp(['Whole Head ...']);
    that=zeros([3 numel(teid)]);
    that(2,:)=1;

    tri_Scalp=surftri(p',te2p'); % Outer surface
    [tri_ROI,node5]=surftri(p',te2p(:,teid)'); % ROI surface

    xx(conductivity(1,:)==.1260)=1; % White
    xx(conductivity(1,:)==.2760)=1; % Grey
    tri_GRM=surftri(p',te2p(:,xx(:)==1)'); % GM/WM surface
    tri_WHM=surftri(p',te2p(:,conductivity(1,:)==0.126)'); % WM surface

%     figure;
%     % Ourter surface
%     % trisurf(%tri,p(1,:)',p(2,:)',p(3,:)','edgealpha',0,'facealpha',.1,'facecolor',[150 114 100]/255); hold on;
%     % GM/WM surface
%     trisurf(tri_GRM,p(1,:)',p(2,:)',p(3,:)','edgealpha',0,'facealpha',1,'facecolor',[242,174,177]/255); hold on;
%     % ROI surface
%     trisurf(tri_ROI,p(1,:)',p(2,:)',p(3,:)','edgealpha',0,'facealpha',1,'facecolor','red'); hold on;
%     % Center of coil locations
%     plot3(roictr(1),roictr(2),0,'w.','markersize',30);
%     axis equal; light; lighting gouraud;
%     view([0,90]); lightangle(gca,40,40);
%     xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
%     set(gca,'fontsize',12);
%     set(gcf,'position',[6 64 980 926]);
%     plot3(pp(:,1),pp(:,2),pp(:,3),'k.');
%     headname=['HeadModel_WholeBrain_sub-',num2str(subn)];
%     saveas(gcf,[headname,'.png']);

    teidN=teid; AnorN=Anor; ppN=pp;


    disp(['Nte: ',num2str(size(teid,1)),' NCoil: ',num2str(size(Anor,3))]);
    pause(0.1);

    omega=dadt;
    disp(['ACA Only U V ...']);
    fname=['ACA_UVMisfit_sub-',num2str(subn),'_WholeBrain.mat'];
%     for nk=10:200
        tic
        [Ux,Vx,MisfitUV,Ik,Jk]=ACA_ADM_ErrUVR(fname,nk,te2p,p,ppN,AnorN,conductivity,teid,that,rs,ks,omega,scth0(end),th_hair,N,FEMord);
        % tri no use in the function
        toc
%         save(fname, 'Ux','Vx','MisfitUV','Ik','Jk','-v7.3');
%     end
%     figure; %('visible','off');
%     semilogy(MisfitUV,'k.-','linewidth',2,'markersize',20); grid on;
%     xlabel('Rank'); ylabel('L^2 Norm Error');
%     title('Convergence Curve_sub-',num2str(subn));
%     set(gca,'fontsize',20);
%     set(gcf,'position',[15,300,700,500]);
%     compname=['ConvergenceCurve_WholeBrain_sub-',num2str(subn)];
%     saveas(gcf,[compname,'.png']);
end
diary off
