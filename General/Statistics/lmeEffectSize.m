function [eta2,CI_eta2] = lmeEffectSize(tbl,Fobs,df1,df2)
% calculated effect size and confidecne interval for linear mixed effects
% model

% define significance
alpha = [2.5 97.5];

% define lamda
lam_hat = Fobs*df1;

% calculate effect size
eta2 = lam_hat/(lam_hat+df2);

% initialize bootstrap
nBoot = 100;
eta2_boot = nan(nBoot,1);
rng(42)

for b = 1:nBoot
    % cluster bootstrap
    ids = unique(tbl.Subject);
    resample_ids = ids(randi(numel(ids),numel(ids),1));
    tbl_b = innerjoin(table(resample_ids,'VariableNames',{'Subject'}),tbl);

    % perform LME
    try
        lme_b = fitlme(tbl_b, 'Y ~ X*Time + (1|Subject)', ...
            'FitMethod','REML','DummyVarCoding','effects');
        an_b = anova(lme_b,'DFMethod','Satterthwaite');
    catch ME
        % handle rank error
        if (strcmp(ME.identifier,...
                'stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:MustBeFullRank_X'))
            continue
        else
            rethrow(ME);
        end
    end

    % calculate current eta2
    idxX = strcmp(an_b.Term,'X');
    F  = an_b.FStat(idxX);
    df1 = an_b.DF1(idxX);
    df2 = an_b.DF2(idxX);
    eta2_boot(b) = (F*df1)/(F*df1+df2);
end

CI_eta2 = prctile(eta2_boot,alpha);

end

