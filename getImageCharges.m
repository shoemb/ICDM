%Used to interface between MATLAB and the Python script that determines the
%image charges
function [q,x,y,z,netCharge] = getImageCharges(epsilonIn,boundariesIn,regionIndicesIn,qinitIn,zinitIn,iterationsIn,isLeadingIn)
    [q,x,y,z,netCharge] = pyrunfile("dielectricInterfacesGeneralizedGetImageCharges.py",["q","x","y","z","netCharge"],epsilon=epsilonIn,boundaries=boundariesIn,regionIndices=regionIndicesIn,qinit=qinitIn,zinit=zinitIn,iterations=iterationsIn,isLeading=isLeadingIn);
    
    q = convertTo2dArray(q);
    x = convertTo2dArray(x);
    y = convertTo2dArray(y);
    z = convertTo2dArray(z);
end

function [arrayOut] = convertTo2dArray(pythonList)
    c = cell(pythonList);
    arrayOut = {};
    for i=1:length(c)
        tempList = [];
        for j=1:length(c{i})
            tempList(end+1) = double(c{i}{j});
        end
        arrayOut{end+1} = tempList;
    end
end
