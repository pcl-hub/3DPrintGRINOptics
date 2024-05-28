clear
% Units in mm and seconds

[filename,savetopath] = uiputfile('*.mat',...
    'Save calculation results to:');

load('BisEMA.mat','BisEMA');
load('startingparameters.mat',...
    'deltas','taus','gsvs'); % deltas = input layer thicknesses, taus = input exposure times, gsvs = input grayscale values

n1 = 1.5650; % maximum index
n2 = 1.5550; % minimum index
discretization = 12; % the number of discrete steps in change of index
grinprofile = 'spherical'; % design of the grin profile
shapeparameters = [4 6.35 3]; % [Effective Diameter, Outter Diameter, Thickness]

qty = 4; % Number of Lenses
spacings = [25 30]; % Spacing between Lenses on the build plate
offcenter = [0 0]; % Center point value for the spacing refrence, [0 0] is the center of the build platform

interlayerparameters = [5 3 5 10]; % [liftdistance, liftspeed, retractspeed, holdtime]

refresh = 1;

calculation = grincalculation(BisEMA,'MonoPrinter1-246',... % (material data, 'Printer name-light intensity')
                n1,n2,discretization,grinprofile,shapeparameters,...
                qty,spacings,offcenter,interlayerparameters,...
                deltas,taus,gsvs,refresh,20,savetopath);
% 20 = vertical sampling rate for conv. prediction (units in microns)

save([savetopath filename],'calculation');

% calculation =  recreateimagesandbuildscript... % location adjustment
%          (calculation, qty, spacings, offcenter,...
%          interlayerparameters, savetopath);

clear filename savetopath BisEMA deltas taus gsvs
clear n1 n2 discretization grinprofile shapeparameters
clear qty spacings offcenter interlayerparameters
