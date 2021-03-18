% test reproducability: rerun the decomposition with decreasing number of
% patient data included in the analysis (6-10 patients) and save the
% signatures. the results will be plotted using fusion_reproducability



% set path
addpath(genpath('../'));
data_path=['../data'];  %%% same as in /esat/biomeddata/_Biomed_INVENTORY/Data/EEG-fMRI/Epilepsy
results_path=['../results'];


mkdir([results_path '/reprod']);



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



saveEEG_tensor=EEG_tensor;
savefMRI_matrix=fMRI_matrix;
% test reproducability agains different patients included -> save results
% -> visualize results using fusion_reproducability.m
for kkk=1:4 %% take different sizes of subsets
    
    
    sel = nchoosek(1:10,10-kkk);  %% all possible subsets of this size
    for sss=1:size(sel,1)
        EEG_tensor=saveEEG_tensor;
        fMRI_matrix=savefMRI_matrix;
        EEG_tensor=EEG_tensor(sel(sss,:),:,:);
        fMRI_matrix=fMRI_matrix(sel(sss,:),:);
        
        
%% JMTF no contraints
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
                model.variables.z = init1{1};% [EEG_tensor(:,9,84)  EEG_tensor(:,9,163)];% rand(size(EEG_tensor,1),R);%init1{1}; % patient variability, factors are shared
                
                
                
                model.variables.w = init2{1}; % voxel signature
                
                
                model.factors.U = 'u';
                model.factors.V = 'v';
                model.factors.Z = 'z';%[EEG_tensor(:,9,84)  EEG_tensor(:,9,163)];%{'z'}
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
                
                
               
             
               
                    save([ results_path '/reprod/jmtf_select' num2str(sss) 'case' num2str(kkk)], 'solution','out' );
               
                
            end
       
        
    
      
        
%% JMTF with fixed signature
        R=2;
           
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
                
                
               
                    save([ results_path '/reprod/jmtffix_select' num2str(sss) 'case' num2str(kkk)], 'solution','out' );
               
                
            end
        
        
       
        
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
        
  
            save([ results_path '/reprod/jica_select' num2str(sss) 'case' num2str(kkk)], 'sources_e','sources_f' );
     
        
        
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
     
        
      
            save([results_path '/reprod/mjica_select' num2str(sss) 'case' num2str(kkk)],  'multi_sources_e','multi_sources_f' );
     
        
     
        
        
    end
end


fusion_reproducability


