%% load the signatues obtained with different number of patients' data included in the fusion and plot the reproducability of the results, quantified as the average pairwise correlation of the signatures

for kkk=1:4
    sel = nchoosek(1:10,10-kkk);
    for sss=1:min(size(sel,1),100)  %% do not compute all the 20000+ possibilities
        
        
        tmp_jica=load(['jica_select' num2str(sss) 'case' num2str(kkk) '.mat']);
        tmp_mjica=load(['mjica_select' num2str(sss) 'case' num2str(kkk) '.mat']);
        tmp_jmtffix=load(['jmtffix_select' num2str(sss) 'case' num2str(kkk) '.mat']);
        tmp_jmtfnoc=load(['jmtf_select' num2str(sss) 'case' num2str(kkk) '.mat']);
        
       
        jica_fmri(sss,:,:)=tmp_jica.sources_f;
        mjica_fmri(sss,:,:)=tmp_mjica.multi_sources_f;
        
        jica_eeg(sss,:,:)=tmp_jica.sources_e;
        mjica_eeg(sss,:,:)=resample(squeeze(tmp_mjica.multi_sources_e(:,9,:))',11923,251);
        
        
        for i=1:5
            fval1(i)=tmp_jmtffix.out{2,i}.fval(end);
            fval2(i)=tmp_jmtfnoc.out{2,i}.fval(end);
        end
        
        [~,ix1]=min(fval1);
        [~,ix2]=min(fval2);
        
        jmtffix_fmri(sss,:,:)=tmp_jmtffix.solution{2,ix1}.factors.W;
        jmtfnoc_fmri(sss,:,:)=tmp_jmtfnoc.solution{2,ix2}.factors.W;
        
        jmtffix_eeg(sss,:,:)=tmp_jmtffix.solution{2,ix1}.factors.V;
        jmtfnoc_eeg(sss,:,:)=tmp_jmtfnoc.solution{2,ix2}.factors.V;
        
    end
    
    sel = nchoosek(1:size(jica_eeg,1),2);
    clear c1 c2 c3 c4
    for xx=1:size(sel,1)
        tmp=abs(corr(squeeze(jica_fmri(sel(xx,1),:,:)),squeeze(jica_fmri(sel(xx,2),:,:))));
        c1(xx)=max(mean(diag(tmp)),mean(diag(fliplr(tmp)))); %% take the diagonal or the antidiagonal, as maybe the order of the sources are not the same
        tmp=abs(corr(squeeze(mjica_fmri(sel(xx,1),:,:)),squeeze(mjica_fmri(sel(xx,2),:,:))));
        c2(xx)=max(mean(diag(tmp)),mean(diag(fliplr(tmp))));
        tmp=abs(corr(squeeze(jmtffix_fmri(sel(xx,1),:,:)),squeeze(jmtffix_fmri(sel(xx,2),:,:))));
        c3(xx)=max(mean(diag(tmp)),mean(diag(fliplr(tmp))));
        tmp=abs(corr(squeeze(jmtfnoc_fmri(sel(xx,1),:,:)),squeeze(jmtfnoc_fmri(sel(xx,2),:,:))));
        c4(xx)=max(mean(diag(tmp)),mean(diag(fliplr(tmp))));
    end
    %figure
    %plot(c1); hold on; plot(c2); plot(c3); plot(c4);
    cjica_fmri(kkk)=mean(c1);
    cmjica_fmri(kkk)=mean(c2);
    cjmtffix_fmri(kkk)=mean(c3);
    cjmtfnoc_fmri(kkk)=mean(c4);
    
    clear c1 c2 c3 c4
    for xx=1:size(sel,1)
        tmp=abs(corr(squeeze(jica_eeg(sel(xx,1),:,:)),squeeze(jica_eeg(sel(xx,2),:,:))));
        c1(xx)=max(mean(diag(tmp)),mean(diag(fliplr(tmp)))); %% take the diagonal or the antidiagonal, as maybe the order of the sources are not the same
        tmp=abs(corr(squeeze(mjica_eeg(sel(xx,1),:,:)),squeeze(mjica_eeg(sel(xx,2),:,:))));
        c2(xx)=max(mean(diag(tmp)),mean(diag(fliplr(tmp))));
        tmp=abs(corr(squeeze(jmtffix_eeg(sel(xx,1),:,:)),squeeze(jmtffix_eeg(sel(xx,2),:,:))));
        c3(xx)=max(mean(diag(tmp)),mean(diag(fliplr(tmp))));
        tmp=abs(corr(squeeze(jmtfnoc_eeg(sel(xx,1),:,:)),squeeze(jmtfnoc_eeg(sel(xx,2),:,:))));
        c4(xx)=max(mean(diag(tmp)),mean(diag(fliplr(tmp))));
    end
    %figure
    %plot(c1); hold on; plot(c2); plot(c3); plot(c4);
    cjica_eeg(kkk)=mean(c1);
    cmjica_eeg(kkk)=mean(c2);
    cjmtffix_eeg(kkk)=mean(c3);
    cjmtfnoc_eeg(kkk)=mean(c4);
    
    
end


figure
subplot(211)
bar([cjica_fmri' cmjica_fmri' cjmtfnoc_fmri' cjmtffix_fmri'])
set(gca,'FontSize',15)
set(gca,'XTickLabel',9:-1:6);
%xlabel('# of patients')
ylabel('mean correlation')
legend({'jointICA';'t-jointICA';'CMTF';'restricted CMTF'},'location','NorthEastOutside')
title('fMRI voxel signature')
ylim([0.5 1])
xlim([0 5])

subplot(212)
bar([cjica_eeg' cmjica_eeg' cjmtfnoc_eeg' cjmtffix_eeg'])
set(gca,'FontSize',15)
set(gca,'XTickLabel',9:-1:6);
xlabel('# of patients')
ylabel('mean correlation')
legend({'jointICA';'t-jointICA';'CMTF';'restricted CMTF'},'location','NorthEastOutside')
title('EEG temporal signature')
ylim([0.5 1])
xlim([0 5])


