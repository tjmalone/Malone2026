function p_value = ttestMan(m1,s1,n1,m2,s2,n2)
%% ttestMan
% performs a two-tailed unpaired student's t-test from the mean, standard
% deviation,and n. For a paired t-ttest, ttest() should be used. This
% function can be used when pre-averaging across groups of unequal number.

% pooled variance
sp = ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2);
SE = (sp*(1/n1 + 1/n2))^0.5;

% t-statistic
t_stat = (m1-m2)/SE;

% calculate degrees of freedom
df = n1+n2-2;

p_value = 2*(1-tcdf(abs(t_stat),df));

end

