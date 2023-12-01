clear 
clc

runCorrection("sampleInputs/inputSampleSodium.txt","sodium-corrected-free-energy.txt");
runCorrection("sampleInputs/inputSampleChloride.txt","chloride-corrected-free-energy");


function runCorrection(inputFname,outputFilename)

    %Read correction parameters and uncorrected free energy profile from
    %files
    readInputFile(inputFname);
    [z,uncorrectedFreeEnergy,uncorrectedFreeEnergyEbars] = readUncorrectedValues(fe_all,id_min,id_max);
 
    %Shift the free energy profile so that it begins at zero
    shiftInFreeEnergy = uncorrectedFreeEnergy(1);
    for i=1:length(uncorrectedFreeEnergy)
        uncorrectedFreeEnergy(i) = uncorrectedFreeEnergy(i) - shiftInFreeEnergy;
    end

    %Optionally shift the value of z_c up or down by the supplied error bar
    if(shiftZc=="down")
        z_c = z_c - z_c_err;
    elseif(shiftZc=="up")
        z_c = z_c + z_c_err;
    end

    %Optionally shift the uncorrected free energy up or down by the
    %supplied error bar
    if(shiftUncorrected=="down")
        uncorrectedFreeEnergy = uncorrectedFreeEnergy - uncorrectedFreeEnergyEbars;
    elseif(shiftUncorrected=="up")
        uncorrectedFreeEnergy = uncorrectedFreeEnergy + uncorrectedFreeEnergyEbars;
    end

    correctedZ = [z(1)];
    correctedFe = [uncorrectedFreeEnergy(1)];
    cumulativeCorrection = [0];  
    netCharges = [0];

    z(1) = z(1) + .0001; %Optional offset to avoid z values at an interface
    for i=2:length(z)
       z(i) = z(i) + .0001; %Optional offset to avoid z values at an interface
       z(i) %Output just to track progress of calculation
      
       %Read in empirical charge density for non-transiting ion
       [pointFractionCharge,pointX,pointY,pointZ,slabFractionCharge,slabZ] = getCharges(OPValues,L_x,L_y,chlorideDensityFilename,z(i),endMembrane);
       

       %Determine which region the transiting ion is in
       regionInd = 'Not Found';
       if(z(i) <= boundaries(1))
            regionInd = 1;
       elseif(z(i) >= boundaries(end))
         regionInd = length(boundaries)+1;
       else
           for j=1:length(boundaries)-1
              if(z(i) > boundaries(j) && z(i) <= boundaries(j+1))
                    regionInd = j+1;
              end
           end
       end
       if(regionInd == 'NotFound')
            "ERROR: regionInd NOT FOUND"
       end

       correction = 0;
       netCharges(end+1) = 0;

       %Correction for charge induced directly by the leading ion
       if(z(i) > feedStartingZ)
           [qim,xim,yim,zim,netCharge] = getImageCharges(epsilon,boundaries,feedRegionIndices,[chargeTransitingIon],[z(i)],numImageChargeIterations,true);
           netCharges(end) = netCharges(end) + netCharge;
           Vcurr = 0;
           Vprev = 0;
           for j=1:length(qim{regionInd})
                Vprev = Vprev + potentialDueToImagesOfPointCharge(chargeTransitingIon,0,0,z(i-1),-L_x/2.0,L_x/2.0,-L_y/2.0,L_y/2.0,xim{regionInd}(j),yim{regionInd}(j),zim{regionInd}(j),XImagesOfPointCharges,YImagesOfPointCharges,qim{regionInd}(j));
                Vcurr = Vcurr + potentialDueToImagesOfPointCharge(chargeTransitingIon,0,0,z(i),-L_x/2.0,L_x/2.0,-L_y/2.0,L_y/2.0,xim{regionInd}(j),yim{regionInd}(j),zim{regionInd}(j),XImagesOfPointCharges,YImagesOfPointCharges,qim{regionInd}(j));
           end
           correction = correction + ((Vcurr - Vprev)/epsilon(regionInd));
       end

        %Correction for non-transiting point charges in the membrane
        Vcurr = 0;
        Vprev = 0;
        for j=1:length(pointX)
            if(z(i) > feedStartingZ)
                [qim,xim,yim,zim,netCharge] = getImageCharges(epsilon,boundaries,feedRegionIndices,[pointFractionCharge(j)],[pointZ(j)],numImageChargeIterations,false);
                netCharges(end) = netCharges(end) + netCharge;
                for k=1:length(qim{regionInd})
                    Vprev = Vprev + potentialDueToImagesOfPointCharge(chargeTransitingIon,0,0,z(i-1),-L_x/2.0,L_x/2.0,-L_y/2.0,L_y/2.0,xim{regionInd}(k),yim{regionInd}(k),zim{regionInd}(k),XImagesOfPointCharges,YImagesOfPointCharges,qim{regionInd}(k));
                    Vcurr = Vcurr + potentialDueToImagesOfPointCharge(chargeTransitingIon,0,0,z(i),-L_x/2.0,L_x/2.0,-L_y/2.0,L_y/2.0,xim{regionInd}(k),yim{regionInd}(k),zim{regionInd}(k),XImagesOfPointCharges,YImagesOfPointCharges,qim{regionInd}(k));
                end
            else
                 Vprev = Vprev + potentialDueToImagesOfPointCharge(chargeTransitingIon,0,0,z(i-1),-L_x/2.0,L_x/2.0,-L_y/2.0,L_y/2.0,pointX(j),pointY(j),pointZ(j),XImagesOfPointCharges,YImagesOfPointCharges,pointFractionCharge(j));
                 Vcurr = Vcurr + potentialDueToImagesOfPointCharge(chargeTransitingIon,0,0,z(i),-L_x/2.0,L_x/2.0,-L_y/2.0,L_y/2.0,pointX(j),pointY(j),pointZ(j),XImagesOfPointCharges,YImagesOfPointCharges,pointFractionCharge(j));     
            end
            correction = correction + ((Vcurr - Vprev)/epsilon(regionInd));
        end


        
        %Correction for uniform slabs from non-transiting ion
        for j=1:length(slabZ)
            Vnet = potentialDueToImagesOfChargedFinitePlane(chargeTransitingIon,z(i-1),0,0,z(i),-L_x/2.0,L_x/2.0,-L_y/2.0,L_y/2.0,slabZ(j),targetXBinSize,targetYBinSize,slabFractionCharge(j),slabZ(j)<z(i-1));
            correction = correction + (Vnet/epsilon(regionInd));
        end

       %Accounting for the "leftover" charge in the feed
       if(z(i) > feedStartingZ)
        excessCharge = totalSystemCharge - chargeTransitingIon - sum(pointFractionCharge) - sum(slabFractionCharge) - netCharges(end);
        Vnet1 = potentialDueToImagesOfChargedFinitePlane(chargeTransitingIon,z(i-1),0,0,z(i),-L_x/2.0,L_x/2.0,-L_y/2.0,L_y/2.0,z_piston,.1,.1,(excessCharge)/2.0,z_piston<z(i-1));
        Vnet2 = potentialDueToImagesOfChargedFinitePlane(chargeTransitingIon,z(i-1),0,0,z(i),-L_x/2.0,L_x/2.0,-L_y/2.0,L_y/2.0,z_c,.1,.1,(excessCharge)/2.0,z_c<z(i-1));
        correction = correction + ((Vnet1 + Vnet2)/epsilon(regionInd));
       end

       correctedZ(end+1) = z(i);
       cumulativeCorrection(end+1) = cumulativeCorrection(end) - correction;
       correctedFe(end+1) = uncorrectedFreeEnergy(i) + cumulativeCorrection(end);
    end


    
    s = length(correctedFe);
    c = {z';uncorrectedFreeEnergy';correctedZ';correctedFe'};
    outputFilename
    writecell(c,outputFilename,'FileType','text');
end

%Read in empirically-determined values for the distribution of charge
%(excluding the transiting ion) outside the feed
function [pointFractionCharge,pointX,pointY,pointZ,slabFractionCharge,slabZ] = getCharges(OPValues,xbox,ybox,chlorideDensityFilename,OP,endMembrane)   
    pointFractionCharge = [];
    pointX = [];
    pointY = [];
    pointZ = [];
    slabFractionCharge = [];
    slabZ = [];
    
    if(length(OPValues)==0 | chlorideDensityFilename=="")
        return
    end

    OPerror = 999999999999;
    ClosestOPIndex = -1;
    for i=1:length(OPValues)
        newError = abs(OP-OPValues(i));
        if(newError < OPerror)
           OPerror = newError;
           ClosestOPIndex = i;
        end
    end
    indexDensityColumn = ClosestOPIndex + 3;
    A = readmatrix(chlorideDensityFilename);
    
    dimensions = size(A);
    for i=1:dimensions(1)
        if(A(i,1) < endMembrane || A(i,2) <= endMembrane)
          pointFractionCharge(end+1) = A(i,indexDensityColumn)*xbox*ybox*(A(i,2) - A(i,1));
          pointX(end+1) = 0;
          pointY(end+1) = 0;
          pointZ(end+1) = mean([A(i,2),A(i,1)]);
       else
          slabFractionCharge(end+1) = A(i,indexDensityColumn)*xbox*ybox*(A(i,2) - A(i,1));
          slabZ(end+1) = A(i,3);
       end
    end
end

%Read in the uncorrected free energy
function [z,uncorrectedFreeEnergy,uncorrectedFreeEnergyEbars] = readUncorrectedValues(fe_all,id_min,id_max)
    ufeAll = [];
    for k=1:length(fe_all)
        A = readmatrix(fe_all(k));
        z = [];
        ufe = [];
        for i=id_min:id_max
           z(end+1) = A(i,1);
           ufe(end+1) = A(i,2);
        end
        ufeAll = [ufeAll;ufe];
    end

    uncorrectedFreeEnergy = [];
    uncorrectedFreeEnergyEbars = [];
    for i=1:length(ufeAll)
        uncorrectedFreeEnergy(end+1) = mean(ufeAll(:,i));
        uncorrectedFreeEnergyEbars(end+1) = std(ufeAll(:,i))/sqrt(length(fe_all));
    end 
end


function V = potentialDueToImagesOfPointCharge(chargeTransitingIon,x,y,z,xp_min,xp_max,yp_min,yp_max,xp,yp,zp,XImages,YImages,q)
    k = 557.005496739; 
    V = 0;

    Lx = xp_max - xp_min;
    Ly = yp_max - yp_min;
    for i=-XImages:XImages
       for j=-YImages:YImages
          if(i==0 & j==0)
              continue
          end
          xt = xp + Lx*i;
          yt = yp + Ly*j;
          zt = zp;
          V = V + k*q*chargeTransitingIon/sqrt((x-xt)^2 + (y-yt)^2 + (z-zt)^2);
       end
    end

end

function V = potentialChangeDueToInfinitePlane(chargeTransitingIon,startingZ,z,xp_min,xp_max,yp_min,yp_max,fractionalCharge,leftOfLeading)
    k = 557.005496739;
    if(leftOfLeading)
        V = -k*2*pi*chargeTransitingIon*(z-startingZ)*fractionalCharge/((xp_max - xp_min)*(yp_max - yp_min));
    else
        V = k*2*pi*chargeTransitingIon*(z-startingZ)*fractionalCharge/((xp_max - xp_min)*(yp_max - yp_min));
    end
end

function V = potentialDueToFinitePlane(chargeTransitingIon,x,y,z,xp_min,xp_max,yp_min,yp_max,zp,targetXBinSize,targetYBinSize,q)
    k = 557.005496739; 
    V = 0;
    
    actualNumXBins = ceil((xp_max - xp_min)/targetXBinSize);
    actualNumYBins = ceil((yp_max - yp_min)/targetYBinSize);
    actualXBinSize = (xp_max - xp_min) / actualNumXBins;
    actualYBinSize = (yp_max - yp_min) / actualNumYBins;
    
    for i=1:actualNumXBins
       for j=1:actualNumYBins
          xval = xp_min + (i-1)*actualXBinSize + (actualXBinSize/2.0);
          yval = yp_min + (j-1)*actualYBinSize + (actualYBinSize/2.0);
          V = V + (actualXBinSize*actualYBinSize)*k*chargeTransitingIon*(q/((yp_max-yp_min)*(xp_max-xp_min)))/sqrt((xval-x)^2 + (yval-y)^2 + (zp-z)^2);
       end
    end
end

function V = potentialDueToImagesOfChargedFinitePlane(chargeTransitingIon,startingZ,x,y,z,xp_min,xp_max,yp_min,yp_max,zp,targetXBinSize,targetYBinSize,fractionalCharge,leftOfLeading)
    deltaVInfinitePlane = potentialChangeDueToInfinitePlane(chargeTransitingIon,startingZ,z,xp_min,xp_max,yp_min,yp_max,fractionalCharge,leftOfLeading);
    VStartingFinitePlane = potentialDueToFinitePlane(chargeTransitingIon,x,y,startingZ,xp_min,xp_max,yp_min,yp_max,zp,targetXBinSize,targetYBinSize,fractionalCharge);
    VEndingFinitePlane = potentialDueToFinitePlane(chargeTransitingIon,x,y,z,xp_min,xp_max,yp_min,yp_max,zp,targetXBinSize,targetYBinSize,fractionalCharge);
    V = deltaVInfinitePlane - (VEndingFinitePlane - VStartingFinitePlane);
end