function conversionprofileprediction = zconversionprofile_matrixcalculation...
                        (formulation,printer,...
                         numberofLayers,layerThicknesses,...
                         exposureTimes,grayscaleValues,...
                         refresh,zinterval,particularlayers,...
                         inclusionofformulationandprinterinfo)
%%
% this function calculates z conversion profile for one pixel
% Instead of a layer-by-layer loop, matrix calculation is used

%% initialize

f = formulation;

ac = f.resinProperties.criticalConversion;
criteffdose = f.resinProperties.criticalEffectiveDose;

Dp = f.workingCurve.depthofPenetration; % resin depth of penetration
Ec = f.workingCurve.criticalEnergy; % resin critical energy

Eoi = f.cureKinetics.oxygenInhibition.oxygenInhibitionEnergy; % oxygen inhibition energy
ck = f.cureKinetics.cureKineticsModel.fitting; % cure kinetics cfit
w = f.cureKinetics.cureKineticsModel.doseRateDependency.average; % dose rate factor of cure kinetics;

gsf = printer.grayscalefunction; 

% Input printing parameters
nol = numberofLayers; % number of layers
deltas = layerThicknesses; % layer thicknesses
taus = exposureTimes; % exposure times
gsvs = grayscaleValues; % grayscale values;

% deltas, taus, and gsvs be n-by-1 vectors
if ~(length(deltas)==nol); deltas = deltas(1)*ones(nol,1); end
if ~(length(taus)==nol); taus = taus(1)*ones(nol,1); end
if ~(length(gsvs)==nol); gsvs = gsvs(1)*ones(nol,1); end


%% calculate conversion profile
dz = zinterval; % sampling points interval along z
z = flip(dz:dz:sum(deltas(:)))'; % sampling z positions for conversion prediction

deltamat = deltas*ones(1,nol);
zsminus = sum(deltamat-triu(deltamat),1);
zmat = z*ones(1,nol)-ones(size(z))*zsminus; % sampling z positions of all stages

I0 = feval(gsf,gsvs); % incident light intensities of all printing stages
I0(I0<0) = 0; % incident light intensities of all printing stages

I0mat = ones(size(z))*I0'; % light intensity matrice for all sampled z positions and all stages
taumat = ones(size(z))*taus'; % exposure time matrice for all sampled z positions and all stages

dosemat = zeros(size(zmat));
dosemat(zmat>0) = I0mat(zmat>0).*taumat(zmat>0).*...
        exp(-zmat(zmat>0)/Dp); % dose matrice for all sampled points of all stages

if logical(refresh) % if the effective dose falls short and resin refreshes, the exposure history is lost

    effdosemat = (dosemat-Eoi).*(I0mat.^(w-1)); % when effective dose is larger than
    printmat = effdosemat>criteffdose;          % the critical effective dose, it prints

    nonrefresh = logical(cumsum(printmat,1,'reverse')); % status of print cumulation in the z direction
    nonrefresh = logical(cumsum(nonrefresh, 2));        % status of print cumulation in the stage direction

    dosemat(~nonrefresh) = 0; % set the dose of all non-print points to 0 due to refreshing setting
end

if ~isempty(particularlayers) % if only the conversions of particular layers are of interest
    
    pts = false(size(z));
    pls = sort(particularlayers);
    sumdeltas = cumsum(deltas,1);
    pInds = floor([sumdeltas-deltas+dz sumdeltas]/dz);
    for i = 1:length(pls)
        pts(pInds(pls(i),1):pInds(pls(i),2),:) = true;
    end
    z = z(pts,:);
    I0mat = I0mat(pts,:);
    dosemat = dosemat(pts,:);
    
end

effdosemat = dosemat;
effdosemat(I0mat>0) = dosemat(I0mat>0).*(I0mat(I0mat>0).^(w-1)); % effective doses

% Do not modify the part between dash lines
% It has to do with determining oxygen inhibition
% ------------------------------------------------
npts = size(z,1);
inds1 = (1:npts)';

cumdosemat = cumsum(dosemat,2);
effdosemat(cumdosemat<=Eoi) = 0;

allois = all(cumdosemat<=Eoi,2);

[~,inds2] = max(cumdosemat>Eoi,[],2);

pts = inds1(~allois)+(inds2(~allois)-1)*npts;

I1 = I0(inds2(~allois));

cumdosevec = cumdosemat(:);
effdosevec = effdosemat(:);

cumresdosevec = cumdosevec(pts)-Eoi;

effdosevec(pts) = cumresdosevec.*(I1.^(w-1));
effdosemat = reshape(effdosevec,npts,nol);

cumeffdosemat = cumsum(effdosemat,2);
% ------------------------------------------------

conversions = feval(ck,cumeffdosemat); % converts accumulated effective doses to conversions using ck (cure kinetics)

stageConversionProfiles = reshape(conversions,npts,nol);
stageConversionProfiles(cumeffdosemat==0) = 0;
stageConversionProfiles(stageConversionProfiles<0) = 0;

%%
z = z/1000;
conversionProfile  = stageConversionProfiles(:,nol);

%%
if isempty(inclusionofformulationandprinterinfo)
    inclusionofformulationandprinterinfo = 1;
end

inclusion = logical(inclusionofformulationandprinterinfo);

if ~inclusion
    formulation = [];
    printer = [];
end

conversionprofileprediction = struct('formulation',formulation, ...
    'printer',printer,'targetConversionProfile',struct('z',[],'conversion',[]),...
    'printParameters', struct('numberofLayers',nol,'layerThicknesses',deltas,...
    'exposureTimes',taus,'grayscaleValues',gsvs,'refresh',refresh),...
    'prediction', struct('particularLayers',particularlayers,...
    'z',z,'effectiveDoseProfile',cumeffdosemat,...
    'stageConversionProfiles',stageConversionProfiles,...
    'conversionProfile',conversionProfile),...
    'units',struct('layerThickness','\{mu}m','exposureTime','s',...
    'z','mm','effectiveDoseProfile',['(mW/cm^{2})^{' num2str(w) '} s']));
