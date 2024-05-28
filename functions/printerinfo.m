function printer = printerinfo(printercode)

switch printercode
    case 'MonoPrinter1-200'        
        f = fittype('a*x+b');
        c = cfit(f,0.001509,-0.063451);
        printer = struct('code','MonoPrinter1-200',...
             'printer','MonoPrinter1','powerLevel',200,...
             'grayscalefunction',c,'pixelSize',64.8);
        
    case 'MonoPrinter1-246'
        f = fittype('a*x+b');
        c = cfit(f,0.004007,-0.177702);
        printer = struct('code','MonoPrinter1-246',...
             'printer','MonoPrinter1','powerLevel',246,...
             'grayscalefunction',c,'pixelSize',64.8);
end