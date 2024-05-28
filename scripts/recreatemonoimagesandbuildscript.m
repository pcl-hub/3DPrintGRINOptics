function calculation =  recreatemonoimagesandbuildscript...
         (calculation, qty, spacings, offcenter, interlayerparameters, savetopath)

%%
predictions = calculation.predictions;
images = calculation.baseImages;

printparameters = [predictions.printParameters];

layerthicknesses= [printparameters.layerThicknesses];
exposuretimes = [printparameters.exposureTimes];
grayscalevalues = [printparameters.grayscaleValues];

deltas= layerthicknesses';
taus = exposuretimes';
gsvs = grayscalevalues';

matsize = size(taus); % size(taus) = size(deltas) = size(gsvs)

zcompensateforglassslide = 0.05;
deltas1 = deltas(:,1)*ones(1,size(deltas,2));
bppositions = (cumsum(deltas,2)- deltas1)/1000;
bppositions(:,end) = bppositions(:,end)+zcompensateforglassslide;
buildplatepositions = bppositions';

[taus1,inds] = sort(taus,1);
stagetaus = diff([zeros(1,matsize(2)); taus1],1);
actualexposuretimes = stagetaus';

imsz = size(images); if length(imsz) == 2; imsz(3) = 1;end
images = reshape(images,[imsz(1)*imsz(2), imsz(3)])';

liftdistance = interlayerparameters(1);
liftspeed = interlayerparameters(2);
retractspeed = interlayerparameters(3);
holdtime = interlayerparameters(4);

maxtau = 600; % assumes maximum layer exposure time possible is 600s
totallayercount = sum(ceil(stagetaus(:)/maxtau));

layercounts = sum(ceil(stagetaus/maxtau),1);
liftdistances = zeros(totallayercount,1);
liftdistances(cumsum(layercounts,2),1) = liftdistance;

%%
fpath = [savetopath 'files\']; mkdir(fpath);
monobuildscriptfileheader(fpath,totallayercount);

fileID = fopen([fpath 'buildscript.ini'],'a');

layercount = 0;

for j = 1:matsize(2)

    gsvs1 = gsvs(inds(:,j),j)*ones(1,imsz(1)*imsz(2));
    
    images1 = images(inds(:,j),:).*gsvs1;
    images1 = cumsum(images1,1,'reverse');

    for i = 1:matsize(1)
        m = ceil(stagetaus(i,j)/maxtau);
        if ~m==0
            bpposition1 = bppositions(i,j);

            image1 = images1(i,:)';
            image1 = reshape(image1,[imsz(1),imsz(2)]);

            fullimage = placetheimages([],image1,qty,spacings,offcenter,[]);
            
            exposuretime1 = ones(m,1)*maxtau;
            exposuretime1(m,1) = stagetaus(i,j)-(m-1)*maxtau;
            
            for i2 = 1:m
                layercount = layercount+1;
                
                filename1 = monoimagename(layercount);
                fullfilename1 = [fpath filename1 '.png'];
                imwrite(uint8(fullimage),fullfilename1,'png');
            
                fprintf(fileID,'%1$.4f,%2$ 10s,%3$ .1f,%4$ .1f,%5$ .1g,%6$ .1f\n',...
                 bpposition1, filename1,exposuretime1(i2,1),liftdistances(layercount,1),0,0);
            end
        end
    end
end

fclose(fileID);

%%
interlayertime = (liftdistance/liftspeed+...
    liftdistance/retractspeed+holdtime)*(size(taus,2)-1);
totalprinttime = (sum(stagetaus(:))+interlayertime)/3600;
lastlayerprinttime = sum(stagetaus(:,end))/60;

disp(['total number of layers is ' num2str(layercount) ';']);
disp(['total print time will be ' num2str(totalprinttime) ' hours; ']);
disp(['last layer will take ' num2str(lastlayerprinttime) ' minutes.']);

%%
printconfiguration = struct('buildscriptInfo',struct(...
                            'buildPlatePositions',buildplatepositions,...
                            'actualExposureTimes',actualexposuretimes,...
                            'units',struct('position','mm','time','second')),...
                            'locationInfo',struct('quantity',qty,...
                            'spacings',spacings,'offcenter',offcenter,...
                            'units',struct('spacing','mm','offcenter','mm')),...
                            'interlayerParameters',struct(...
                            'liftDistance', liftdistance,'liftSpeed', liftspeed,...
                            'retractSpeed',retractspeed,'holdTime', holdtime,...
                            'units', struct('liftDistance', 'mm', ...
                            'speed', 'mm/s', 'holdTime', 'second')));

calculation.printConfiguration = printconfiguration;


function monobuildscriptfileheader(fpath,totallayercount)

fileID = fopen([fpath 'buildscript.ini'],'w');
fprintf(fileID, '%1$20s\n','Machine = MONO3_64UM');
fclose(fileID);

strnumber = floor(log10(totallayercount)+20);
str1 = ['%1$' num2str(strnumber) 's\n'];
str2 = num2str(totallayercount);

fileID = fopen([fpath 'buildscript.ini'],'a');
fprintf(fileID, '%1$21s\n','Slice thickness = 100');
fprintf(fileID, str1,['number of slices = ' str2]);
fprintf(fileID, '%1$22s\n','illumination time = 20');
fprintf(fileID, '%1$29s\n','number of override slices = 0');
fprintf(fileID, '%1$31s\n','override illumination time = 20');
fprintf(fileID, '%1$25s\n','support burn in time = 20');
fclose(fileID);


function filename = monoimagename(layernumber)
numdigits = floor(log10(layernumber)+1);
monoprefix = {'S00000' 'S0000' 'S000' 'S00' 'S0'};
filename = [monoprefix{numdigits} num2str(layernumber) '_P1'];