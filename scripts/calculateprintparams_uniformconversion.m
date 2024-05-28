function cp = calculateprintparams_uniformconversion...
                (formulation,printer,targetthickness,targetconversion,...
                inputdeltas,inputtaus,inputgsvs,refresh,zinterval)

generalInfo = struct('formulation',formulation,...
                     'printer', printer,...
                     'targetconversion',targetconversion,...
                     'refresh',logical(refresh),...
                     'zinterval',zinterval);

printcompletion = false;
currentNumberofLayers = 0;

[inputmat, inputcells, m] = generatecombinations...
    (inputdeltas,inputtaus,inputgsvs);

mnol = floor(targetthickness/zinterval*1000); % maxNumberofLayers
writetomat = zeros(mnol,3); % delta, tau, gsv

[generalcells{1:m,1}] = deal(generalInfo);

tic
while ~printcompletion

    currentNumberofLayers = currentNumberofLayers+1;
    
    cnol = currentNumberofLayers;

    [cnolcells{1:m,1}] = deal(cnol);
    
    [writetocells{1:m,1}] = deal(writetomat);

    ase = cellfun(@calculateASE,generalcells,...
                  cnolcells,inputcells,writetocells);

    [~,ind] = min(ase(:));
    
    layer1 = mnol-cnol+1;
    layers = layer1:mnol;

    writetomat(layer1,:) = inputmat(ind,:);

    currentthickness = sum(writetomat(:,1),1)/1000;

    if currentthickness == targetthickness

        printcompletion = true;

    elseif currentthickness > targetthickness

        delta = targetthickness*1000-...
            sum(writetomat(layer1+1:mnol,1),1);

        [inputmat, inputcells, m] = generatecombinations...
                (delta,inputtaus,inputgsvs);
        
        [generalcells{1:m,1}] = deal(generalInfo);

        [cnolcells{1:m,1}] = deal(cnol);

        [writetocells{1:m,1}] = deal(writetomat);
    
        ase = cellfun(@calculateASE,generalcells,...
                      cnolcells,inputcells,writetocells);

        [~,ind] = min(ase(:));
    
        writetomat(layer1,:) = inputmat(ind,:);

        printcompletion = true;
    end
end
toc

deltas = writetomat(layers,1);
taus = writetomat(layers,2);
gsvs = writetomat(layers,3);

cp = zconversionprofile_matrixcalculation....
        (formulation,printer,...
         cnol,deltas,taus,gsvs,...
         refresh,zinterval,[],0);

cp.targetConversionProfile.z = cp.prediction.z;
cp.targetConversionProfile.conversion = targetconversion*...
    ones(size(cp.prediction.z));



%%
function [combinationmat, combinationcells, m] = ...
            generatecombinations(deltas,taus,gsvs)

ld = length(deltas);
lt = length(taus);
lg = length(gsvs);

m = ld*lt*lg;

darray = repmat(deltas',lt*lg,1);
tarray = repmat(taus',lg*ld,1);
garray = repmat(gsvs',ld*lt,1);

combinationmat = [darray tarray garray];
combinationcells = mat2cell(combinationmat,ones(m,1));



function ase = calculateASE(generalInfo,currentNumberofLayers,...
                            inputparameters,writetomat)

formulation = generalInfo.formulation;
printer = generalInfo.printer;
targetconversion = generalInfo.targetconversion;
refresh = generalInfo.refresh;
zinterval = generalInfo.zinterval;

mnol = size(writetomat,1);
cnol = currentNumberofLayers;

layer1 = mnol-cnol+1;
layers = layer1:mnol;

writetomat(layer1,:) = inputparameters;

deltas = writetomat(layers,1);
taus = writetomat(layers,2);
gsvs = writetomat(layers,3);

cpstruct = zconversionprofile_matrixcalculation....
                (formulation,printer,...
                cnol,deltas,taus,gsvs,...
                refresh,zinterval,1,0);

cp = cpstruct.prediction.conversionProfile;
cp0 = targetconversion*ones(size(cp));
ase = (cp-cp0)'*(cp-cp0)/size(cp,1);
