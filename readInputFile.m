function readInputFile(filename)
    A = readcell(filename);
    
    totalSystemCharge = NaN;
    chargeTransitingIon = NaN;
    startingZ = NaN;
    feedStartingZ = NaN;
    endMembrane = NaN;
    fe_all = [];
    chlorideDensityFilename = "";
    OPValues = [];
    z_c = NaN;
    z_c_err = NaN;
    shiftZc = "";
    shiftUncorrected = "";
    L_x = NaN;
    L_y = NaN;
    id_min = NaN;
    id_max = NaN;
    z_piston = NaN;
    boundaries = [];
    feedRegionIndices = [];
    epsilon = [];
    numImageChargeIterations = NaN;
    XImagesOfPointCharges = NaN;
    YImagesOfPointCharges = NaN;
    targetXBinSize = NaN;
    targetYBinSize = NaN;
    
    for i=1:length(A)
        if ismissing(A{i,3}) | A{i,2}~="="
            error("ERROR...INCORRECT FORMAT FOR INPUT FILE AT LINE NUMBER: " + num2str(i));
        end
    
        d = A{i,1};
        val = A{i,3};
        if(d=="totalSystemCharge")
            totalSystemCharge = val;
        elseif(d=="chargeTransitingIon")
            chargeTransitingIon = val;
        elseif(d=="startingZ")
            startingZ = val;
        elseif(d=="feedStartingZ")
            feedStartingZ = val;
        elseif(d=="endMembrane")
            endMembrane = val;
        elseif(d=="free-energy-profile")
            fe_all = [fe_all;convertCharsToStrings(val)];
        elseif(d=="chlorideDensityFilename")
            chlorideDensityFilename = val;
        elseif(d=="OPValues")
            OPValues = str2num(val);
        elseif(d=="z_c")
            z_c = val;
        elseif(d=="z_c_err")
            z_c_err = val;
        elseif(d=="shiftZc")
            shiftZc = val;
        elseif(d=="shiftUncorrected")
            shiftUncorrected = val;
        elseif(d=="L_x")
            L_x = val;
        elseif(d=="L_y")
            L_y = val;
        elseif(d=="id_min")
            id_min = val;
        elseif(d=="id_max")
            id_max = val;
        elseif(d=="z_piston")
            z_piston = val;
        elseif(d=="boundaries")
            boundaries = str2num(val);
        elseif(d=="feedRegionIndices")
            feedRegionIndices = str2num(val);
        elseif(d=="epsilon")
            epsilon = str2num(val);
        elseif(d=="numImageChargeIterations")
            numImageChargeIterations = val;
        elseif(d=="XImagesOfPointCharges")
            XImagesOfPointCharges = val;
        elseif(d=="YImagesOfPointCharges")
            YImagesOfPointCharges = val;
        elseif(d=="targetXBinSize")
            targetXBinSize = val;
        elseif(d=="targetYBinSize")
            targetYBinSize = val;
        end
        
        assignin('caller','totalSystemCharge',totalSystemCharge);
        assignin('caller',"chargeTransitingIon",chargeTransitingIon);
        assignin('caller',"startingZ",startingZ);
        assignin('caller',"feedStartingZ",feedStartingZ);
        assignin('caller',"endMembrane",endMembrane);
        assignin('caller',"fe_all",fe_all);
        assignin('caller',"chlorideDensityFilename",chlorideDensityFilename);
        assignin('caller',"OPValues",OPValues);
        assignin('caller',"z_c",z_c);
        assignin('caller',"z_c_err",z_c_err);
        assignin('caller',"shiftZc",shiftZc);
        assignin('caller',"shiftUncorrected",shiftUncorrected);
        assignin('caller',"L_x",L_x);
        assignin('caller',"L_y",L_y);
        assignin('caller',"id_min",id_min);
        assignin('caller',"id_max",id_max);
        assignin('caller',"z_piston",z_piston);
        assignin('caller',"boundaries",boundaries);
        assignin('caller',"feedRegionIndices",feedRegionIndices);
        assignin('caller',"epsilon",epsilon);
        assignin('caller',"numImageChargeIterations",numImageChargeIterations);
        assignin('caller',"XImagesOfPointCharges",XImagesOfPointCharges);
        assignin('caller',"YImagesOfPointCharges",YImagesOfPointCharges);
        assignin('caller',"targetXBinSize",targetXBinSize);
        assignin('caller',"targetYBinSize",targetYBinSize);
    end
end
