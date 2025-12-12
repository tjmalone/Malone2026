function B = matOffset(A,offsets)
%% matOffset.m
% Uniformly shifts matrix along n dimensions. Missing data are filled in
% with zeros. Variable offsets must be part of a loop.

sz=size(A);

N=length(sz);

if length(offsets)<N
    offsets(N)=0;
end

B=zeros(sz);

indices=cell(3,N);

for ii=1:N
    for ss=[1,3]
        idx=(1:sz(ii))+(ss-2)*offsets(ii);
        idx(idx<1)=[];
        idx(idx>sz(ii))=[];
        indices{ss,ii}=idx;
    end
end

src_indices=indices(1,:);
dest_indices=indices(3,:);
B(dest_indices{:}) = A(src_indices{:});

end