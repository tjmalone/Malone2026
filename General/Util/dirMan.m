%%

a = 210604;
b = 210604;
st = 1;

dates = a:st:b;

for i = 1:length(dates)
    cur = num2str(dates(i));
    mkdir(cur)
    
    cd(cur)
    
    mkdir('loc1')
    mkdir('loc2')
    mkdir('loc2')
    
    cd ..
end

