mkdir('Release');

x = input('Compile bm_polar_transform?', 's');
if strcmp(x, 'y')
    mex -R2018a -DON_LINUX -silent bm_polar_transform/mexmain.cpp bm_polar_transform/pch.cpp ...
        -output Release/bm_polar_transform
    disp('BM-Polar-Transform MEX successfully generated at ./Release/');
end

x = input('Compile Q-ary SC Decoder?', 's');
if strcmp(x, 'y')
    mex -R2018a -DON_LINUX -silent Qary_SC_Decoder/dllmain.cpp Qary_SC_Decoder/pch.cpp PolarCpp/GF.cpp PolarCpp/Qary_dist.cpp PolarCpp/SC.cpp PolarCpp/SCFrame.cpp PolarCpp/SCList.cpp ...
        -output Release/Qary_SC_Decoder
    disp('Q-ary SC Decoder MEX successfully generated at ./Release/');
end

