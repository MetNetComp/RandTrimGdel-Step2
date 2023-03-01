function [outputArg1,outputArg2] = example()
initCobraToolbox;
load('gDel-minRN-strategy.mat');
load('e_coli_core.mat');
model=e_coli_core;
[gvalueList] = randTrimGdelStep2(model,'succ_e',gvalue,5)

save('example.mat');
end

