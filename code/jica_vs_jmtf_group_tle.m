% EEG-fMRI data from 10 temporal lobe epilepsy patients
% average EEG waveforms and GLM-based activation maps are taken per
% patient. The activation maps are masked based on the average activation
% map

% jointICA, t-jointICA, CMTF and CMTF with fixed patient signatures is
% compared. The patient signature is fixed based on the spike and slow wave
% amplitude of each patient.

% saves signatures, plots EEG temporal signatures and plots the distribution
% of the fMRI voxel signatures

function jica_vs_jmtf_group_tle

% set path
addpath(genpath('../'));
data_path=['../data'];  %%% same as in /esat/biomeddata/_Biomed_INVENTORY/Data/EEG-fMRI/Epilepsy
results_path=['../results'];
mkdir([results_path '/jica']);
mkdir([results_path '/jmtf']);



%% PREPARE DATA
left=[1 1 0   0   0  1  1  0   1   0]; %% invert left temporal lobe data


avg_act=0;
fMRI_matrix=[];
% lot GLM-based activation maps
for i=1:10%numel(pats)
    
    try
        vv=spm_vol([data_path '/patient' num2str(i) '/spmT_0001.nii']);
    catch
        vv=spm_vol([data_path '/patient' num2str(i) '/spmT_0001.hdr']);
    end
    vv=spm_read_vols(vv);
    
    mid1=(size(vv,1)-1)/2;
    mid2=mid1+2;
    if left(i) %% mirror the actications maps of the left ones
        for n=1:size(vv,2)
            for j=1:size(vv,3)
                for k=-mid1:1:mid1
                    vvm(mid1+1-k,n,j)=vv(mid2-1+k,n,j);
                end
            end
        end
        vv=vvm;
    end
    vv(find(isnan(vv(:))))=0;
    fMRI_matrix(i,:)=vv(:);
    avg_act=avg_act+vv;
end
avg_act=reshape(zscore(avg_act(:)/10),[79 95 68]);

% save the average activation
 try
        vv=spm_vol([data_path '/patient1/spmT_0001.nii']);
    catch
        vv=spm_vol([data_path '/patient1/spmT_0001.hdr']);
    end
vv.fname=[results_path '/avg_mirrored_spm.img'];
spm_write_vol(vv,avg_act);
mask_ind=find(abs(avg_act)>3);
fMRI_matrix=fMRI_matrix(:,mask_ind);


% create EEG tensor from spikes

for i=1:10%numel(pats)
    
    load([data_path '/patient' num2str(i) '/average_spike.mat']);
    if left(i) %% mirror the channels of left TLE patients
        
        allspikes(i,:,:)=spikes_avg2d;
    else
        allspikes(i,:,:)=spikes_avg2d([12:19 21 10 11 1:8 20 9],: ); % take channels which are common among all (internation 10-20 with 19 electrodes + 2 sphenoidal)
    end
    
end
clear EEG_tensor;
sixth=size(allspikes,3)/6; % to reduce the baseline period, so the epoch is [-330 660]ms
for i=1:10%numel(pats)
    
    EEG_tensor(i,:,:)=allspikes(i,:,sixth+1:4*sixth);
end

%% FUSION



%% JMTF no constraints
R=2;


clear model;
for iter=1:5
    size_tens=size(EEG_tensor);
    clear options;
    options.Imag=0; options.Real=@randn;
    init1=cpd_rnd(size_tens,R,options);
    init2=cpd_rnd(size(fMRI_matrix,2),R,options);
    clear options;
    
    model.variables.u = init1{2}; % channel signature
    model.variables.v = init1{3}; % time course signature
    model.variables.z = init1{1};%  patient variability, factors are shared
    
    
    
    model.variables.w = init2{1}; % voxel signature
    
    
    model.factors.U = 'u';
    model.factors.V = 'v';
    model.factors.Z = 'z';
    model.factors.W = 'w';
    
    
    model.factorizations.xfac.data = EEG_tensor;
    model.factorizations.xfac.cpd = {'Z','U','V'};
    
    model.factorizations.yfac.data = fMRI_matrix;
    model.factorizations.yfac.cpd = {'Z','W'};
    
    
    
    options.MaxIter=1000;
    options.Display = 5;
    options.TolFun = 1e-20;
    options.TolX = 1e-10 ;
    
    [solution{R,iter},out{R,iter}] = sdf_nls(model,options);
    fval(iter)=out{R,iter}.fval(end);
    mspikes=squeeze(mean(EEG_tensor,1));
    
    % if time course of sources are flipped, flip back
    for i=1:R
        ctc=corr(-solution{R,iter}.factors.V(:,i),mspikes(9,:)');
        csp=corr(-solution{R,iter}.factors.U(:,i),mspikes(:,84));
        
        if ctc>0 && csp>0
            solution{R,iter}.factors.U(:,i)=-solution{R,iter}.factors.U(:,i);
            solution{R,iter}.factors.W(:,i)=-solution{R,iter}.factors.W(:,i);
        elseif ctc>0 && csp<0
            % do nothing, this is correct
        elseif ctc<0 && csp>0
            solution{R,iter}.factors.V(:,i)=-solution{R,iter}.factors.V(:,i);
            solution{R,iter}.factors.U(:,i)=-solution{R,iter}.factors.U(:,i);
        elseif ctc<0 && csp<0
            solution{R,iter}.factors.V(:,i)=-solution{R,iter}.factors.V(:,i);
            solution{R,iter}.factors.W(:,i)=-solution{R,iter}.factors.W(:,i);
        end
    end
    
    
    
    
end


[~,ix]=min(fval);
sol_jmtf_noc=solution{R,ix};


sources=solution{R,ix}.factors.V;
plotTemporalSignature(sources,'CMTF',EEG_tensor)


%% JMTF with fixed signature


clear model;
for iter=1:5
    size_tens=size(EEG_tensor);
    clear options;
    options.Orth=1; options.Imag=0; options.Real=@randn;
    init1=cpd_rnd(size_tens,R,options);
    init2=cpd_rnd(size(fMRI_matrix,2),R,options);
    clear options;
    
    model.variables.u = init1{2}; % channel signature
    model.variables.v = init1{3}; % time course signature
    model.variables.z =  [EEG_tensor(:,9,84)  EEG_tensor(:,9,163)];% rand(size(EEG_tensor,1),R);%init1{1}; % patient variability, factors are shared
    
    
    
    model.variables.w = init2{1}; % voxel signature
    
    
    model.factors.U = 'u';
    model.factors.V = 'v';
    model.factors.Z = [EEG_tensor(:,9,84)  EEG_tensor(:,9,163)];%{'z'}
    model.factors.W = 'w';
    
    
    model.factorizations.xfac.data = EEG_tensor;
    model.factorizations.xfac.cpd = {'Z','U','V'};
    
    model.factorizations.yfac.data = fMRI_matrix;
    model.factorizations.yfac.cpd = {'Z','W'};
    
    
    
    options.MaxIter=1000;
    options.Display = 5;
    options.TolFun = 1e-20;
    options.TolX = 1e-10 ;
    
    [solution{R,iter},out{R,iter}] = sdf_nls(model,options);
    fval(iter)=min(out{R,iter}.fval(end));
    mspikes=squeeze(mean(EEG_tensor,1));
    
    % if time course of sources are flipped, flip back
    for i=1:R
        ctc=corr(-solution{R,iter}.factors.V(:,i),mspikes(9,:)');
        csp=corr(-solution{R,iter}.factors.U(:,i),mspikes(:,84));
        
        if ctc>0 && csp>0
            solution{R,iter}.factors.U(:,i)=-solution{R,iter}.factors.U(:,i);
            solution{R,iter}.factors.W(:,i)=-solution{R,iter}.factors.W(:,i);
        elseif ctc>0 && csp<0
            % do nothing, this is correct
        elseif ctc<0 && csp>0
            solution{R,iter}.factors.V(:,i)=-solution{R,iter}.factors.V(:,i);
            solution{R,iter}.factors.U(:,i)=-solution{R,iter}.factors.U(:,i);
        elseif ctc<0 && csp<0
            solution{R,iter}.factors.V(:,i)=-solution{R,iter}.factors.V(:,i);
            solution{R,iter}.factors.W(:,i)=-solution{R,iter}.factors.W(:,i);
        end
    end
    
    
    
    save([results_path '/jmtf/resultsR' num2str(R) 'spmt_fix'], 'solution','out' );
    
        vv=spm_vol([data_path '/patient1/spmT_0001.hdr']);
    
    for i=1:R
        
        vv.fname=[results_path '/jmtf/spmt_fix_R' num2str(R) '_s' num2str(i) '_iter' num2str(iter) '.img'];
        
        fmri_comp=zeros(vv.dim);
        fmri_comp(mask_ind)=zscore(solution{R,iter}.factors.W(:,i));
        spm_write_vol(vv,fmri_comp);
        
    end
    
end


[m,ix]=min(fval);
sol_jmtf_fix=solution{R,ix};

sources=solution{R,ix}.factors.V;
plotTemporalSignature(sources,'restricted CMTF',EEG_tensor)





%% jointICA

input=[fMRI_matrix resample(squeeze(EEG_tensor(:,9,:))',size(fMRI_matrix,2),size(EEG_tensor,3))'];
fMRI=fMRI_matrix;
EEG=EEG_tensor;

nic=2;




[iq, A, W, S, sR]=icasso(input,50,'lastEig',nic,'vis','off');

mixing_matrix=A;
sources=S;
sources=sources';


sources_e=sources(size(fMRI,2)+1:end,:);
sources_f=sources(1:size(fMRI,2),:);

for ix=1:nic
    invert=sign(corr(resample(sources_e(:,ix),size(EEG_tensor,3),size(sources_e,1)),squeeze(mean(EEG_tensor(:,9,:),1))));
    sources_e(:,ix)=sources_e(:,ix)*invert;
    sources_f(:,ix)=sources_f(:,ix)*invert;
end


   vv=spm_vol([data_path '/patient1/spmT_0001.hdr']);
for i=1:size(sources_e,2);
    vv.fname=[results_path '/jica/spmt_R' num2str(nic) '_s' num2str(i) '.img'];
    fmri_comp=zeros(vv.dim);
    fmri_comp(mask_ind)=zscore(sources_f(:,i));
    spm_write_vol(vv,reshape(fmri_comp,[79 95 68]));
end



sources=(resample(sources_e,size(EEG_tensor,3),size(sources_e,1)));
plotTemporalSignature(sources,'jointICA',EEG_tensor)




%% t-jointICA


nic=2;




rEEG=reshape(permute(EEG_tensor,[1,3,2]),size(EEG_tensor,1),size(EEG_tensor,2)*size(EEG_tensor,3));
rEEG=resample(rEEG',2,1)';
multi_input=[fMRI_matrix rEEG]; %% they are very similar size

% fastICA
[iq, A, W, S, sR]=icasso(multi_input,50,'lastEig',nic,'approach','symm','vis','off'); % OK
multi_mixing_matrix=A;
multi_sources=S';

multi_sources_e=multi_sources(size(fMRI_matrix,2)+1:end,:); multi_sources_e=resample(multi_sources_e,1,2);
multi_sources_e=permute(reshape(multi_sources_e',size(multi_sources_e,2),size(EEG_tensor,3),size(EEG_tensor,2)),[1,3,2]);
multi_sources_f=multi_sources(1:size(fMRI_matrix,2),:);




for i=1:nic
    invert=sign(corr(squeeze(multi_sources_e(i,9,:)),squeeze(mean(EEG_tensor(:,9,:),1))));
    multi_sources_e(:,i)=multi_sources_e(:,i)*invert;
    multi_sources_f(:,i)=multi_sources_f(:,i)*invert;
end



sources=squeeze(multi_sources_e(:,9,:));
plotTemporalSignature(sources,'t-jointICA',EEG_tensor)




%% visualization
% order sources 1.spike 2. slow wave
sources_e=resample(sources_e,size(EEG_tensor,3),size(sources_e,1));
clear tmp;
tmp=corr(sources_e,sol_jmtf_fix.factors.V); tmp=tmp(1,:); [~,order_jica]=sort(abs(tmp),'descend');
tmp=corr(squeeze(multi_sources_e(:,9,:))',sol_jmtf_fix.factors.V); tmp=tmp(1,:); [~,order_mjica]=sort(abs(tmp),'descend');
tmp=corr(sol_jmtf_noc.factors.V,sol_jmtf_fix.factors.V); tmp=tmp(1,:); [~,order_jmtf]=sort(abs(tmp),'descend');

%% distribution plot: voxel distibution of sources in each ROI


% get ROIs

vv=spm_vol([data_path '/ROI/rdDMN_cort_subcort.nii']);
dDMN=spm_read_vols(vv); dDMN(abs(avg_act)<3)=0;


vv=spm_vol([data_path '/ROI/rRECN_cort_subcort.nii']);
recn=spm_read_vols(vv); recn(abs(avg_act)<3)=0;



vv=spm_vol([data_path '/ROI/rLECN_cort_subcort.nii']);
lecn=spm_read_vols(vv); lecn(abs(avg_act)<3)=0;


ecn=lecn|recn;



vv=spm_vol([data_path '/ROI/rvDMN_cort_subcort.nii']);
vDMN=spm_read_vols(vv); vDMN(abs(avg_act)<3)=0;

Sf1=spm_vol([data_path '/ROI/s1.nii']);
Sf1=spm_read_vols(Sf1); Sf1=(Sf1>3);

Sf2=spm_vol([data_path '/ROI/s2.nii']);
Sf2=spm_read_vols(Sf2); Sf2=(Sf2>3);

Sf3=spm_vol([data_path '/ROI/s3.nii']);
Sf3=spm_read_vols(Sf3); Sf3=(Sf3>3);

clear tmp;
tmp{1}=find(Sf1==1);
tmp{2}=find(Sf2==1);
tmp{3}=find(Sf3==1);
tmp{4}=find(dDMN==1);
tmp{5}=find(ecn==1);

figure
clear fmri_comp
subplot(141)
for i=1:2
    fmri_comp{i}=zeros(vv.dim);
    fmri_comp{i}(mask_ind)=zscore(sources_f(:,order_jica(i)));
end
compareDistribPlot(fmri_comp, tmp)
title('jointICA')
ylim([-3.2 3.2])


subplot(142)
for i=1:2
    fmri_comp{i}=zeros(vv.dim);
    fmri_comp{i}(mask_ind)=zscore(multi_sources_f(:,order_mjica(i)));
end
compareDistribPlot(fmri_comp, tmp)
title('t-jointICA')
ylim([-3.2 3.2])


subplot(143)
for i=1:2
    fmri_comp{i}=zeros(vv.dim);
    fmri_comp{i}(mask_ind)=zscore(sol_jmtf_noc.factors.W(:,order_jmtf(i)));
end
compareDistribPlot(fmri_comp, tmp)
title('jMTF')
ylim([-3.2 3.2])


subplot(144)
for i=1:2
    fmri_comp{i}=zeros(vv.dim);
    fmri_comp{i}(mask_ind)=zscore(sol_jmtf_fix.factors.W(:,i));
end
compareDistribPlot(fmri_comp, tmp)
title('fixed jMTF')
ylim([-3.2 3.2])


function plotTemporalSignature(sources,figure_title,EEG_tensor)


figure

for i=1:size(sources,2);
    subplot(size(sources,2),1,i)
    plot(squeeze(mean(EEG_tensor(:,9,:),1)),'Color',[0.5 0.5 0.5])
    hold on
    tmp_source=sources(:,i);
    [ymax,iy]=max(abs(tmp_source));
    ymax=ymax*sign(tmp_source(iy));
    
    if iy<94
        yrat=mean(EEG_tensor(:,9,84))/ymax;
    elseif iy<132
        yrat=mean(EEG_tensor(:,9,103))/ymax;
    else
        yrat=mean(EEG_tensor(:,9,162))/ymax;
    end
    plot(tmp_source*yrat,'k')
    ylabel(['source ' num2str(i)])
    xlim([0 250])
    ylim([-max(abs(tmp_source*yrat)) max(abs(tmp_source*yrat))])
    set(gca,'XTick',[50 175])
    set(gca,'XTickLabel',{'0';'500'})
end
xlabel('time (ms)')
title(figure_title)

function compareDistribPlot(fmri_comp, tmp)



data{1}=fmri_comp{1}(tmp{1});
data{2}=fmri_comp{2}(tmp{1});
data{3}=[];

data{4}=fmri_comp{1}(tmp{2});
data{5}=fmri_comp{2}(tmp{2});
data{6}=[];

data{7}=fmri_comp{1}(tmp{3});
data{8}=fmri_comp{2}(tmp{3});
data{9}=[];

data{10}=fmri_comp{1}(tmp{4});
data{11}=fmri_comp{2}(tmp{4});
data{12}=[];

data{13}=fmri_comp{1}(tmp{5});
data{14}=fmri_comp{2}(tmp{5});

distributionPlot(data,'color',{'b','r','k','b','r','k','b','r','k','b','r','k','b','r'},'showMM',0,'globalNorm',0);
hold on
for i=1:numel(data)
    if max(data{i})>2
        plot(i,3,'*','Color',[1 1 1]/(numel(find(data{i}>2)))/10)
    end
    if min(data{i})<-2
        plot(i,-3,'*','Color',[1 1 1]/(numel(find(data{i}<-2)))/10)
    end
end


set(gca,'XTick',1.5:3:15);
set(gca,'XTickLabel',{'rTL';'lTL';'OL';'DMN';'ECN'})


