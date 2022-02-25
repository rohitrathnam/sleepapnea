%% train/test per subject

train_set = [];
test_set = [];

for i=1:5
    [train,test] = datasplit(i);
    drop = [9;11;31;32;33;34;35;40;49;57];
    train(:,drop)=[];
    test(:,drop)=[];
    train_set = [train_set; train];
    test_set = [test_set; test];
end

%%
a=0;
for i=1:20:100
    yfit = boosttree.predictFcn(test_set(i:i+19,1:end-1));
    err = yfit-test_set(i:i+19,end);
    err=abs(err);
    a = a+(length(err)-sum(err))/length(err)*100
end
disp(a/5)